#!/bin/bash
set -euo pipefail

ROOT=$(cd "$(dirname "$0")/.." && pwd)
cd "$ROOT"

python3 - <<'PY'
from __future__ import annotations

from pathlib import Path
import re
import subprocess
import sys

SKIP_SUFFIXES = {
    ".png", ".jpg", ".jpeg", ".gif", ".webp", ".pdf",
    ".pack", ".idx", ".rev", ".pyc", ".zip", ".gz", ".xz",
}
SKIP_PREFIXES = (Path("smoke_tests/tmp"),)


def candidate_paths() -> list[Path]:
    """Return tracked and non-ignored files without scanning .git internals."""
    try:
        output = subprocess.check_output(
            ["git", "ls-files", "-z", "--cached", "--others", "--exclude-standard"]
        )
        return [
            Path(item)
            for item in output.decode("utf-8", errors="surrogateescape").split("\0")
            if item
        ]
    except (OSError, subprocess.CalledProcessError):
        return [path for path in Path(".").rglob("*") if path.is_file()]


# Build sensitive path roots from fragments so this scanner can safely scan itself.
private_roots = [
    "/" + "home/",
    "/" + "Users/",
    "/" + "mnt/",
    "/" + "share/",
    "/" + "opt/",
]

rules: list[tuple[str, re.Pattern[str]]] = [
    (
        "email address",
        re.compile(r"(?i)\b[A-Z0-9._%+-]+@[A-Z0-9.-]+\.[A-Z]{2,}\b"),
    ),
    (
        "IPv4 address",
        re.compile(
            r"(?<![0-9.])(?:25[0-5]|2[0-4][0-9]|1?[0-9]{1,2})"
            r"(?:\.(?:25[0-5]|2[0-4][0-9]|1?[0-9]{1,2})){3}(?![0-9.])"
        ),
    ),
    ("GitHub token", re.compile(r"\b(?:gh[pousr]_[A-Za-z0-9]{20,}|github_pat_[A-Za-z0-9_]{40,})\b")),
    ("OpenAI-style API key", re.compile(r"\bsk-(?:proj-)?[A-Za-z0-9_-]{20,}\b")),
    ("AWS access key", re.compile(r"\b(?:AKIA|ASIA)[0-9A-Z]{16}\b")),
    ("Google API key", re.compile(r"\bAIza[0-9A-Za-z_-]{35}\b")),
    ("Slack token", re.compile(r"\bxox[baprs]-[A-Za-z0-9-]{20,}\b")),
    (
        "private key material",
        re.compile(r"-----BEGIN (?:RSA |OPENSSH |EC |DSA )?PRIVATE KEY-----"),
    ),
    (
        "SSH public key material",
        re.compile(r"\bssh-(?:rsa|ed25519)\s+[A-Za-z0-9+/]{80,}={0,3}\b"),
    ),
    (
        "credential embedded in URL",
        re.compile(r"(?i)https?://[^\s/:@]+:[^\s/@]+@"),
    ),
    (
        "literal secret assignment",
        re.compile(
            r"(?i)\b(?:password|passwd|api[_-]?key|access[_-]?token|secret)\b"
            r"\s*[:=]\s*['\"](?!<)[^'\"\s]{8,}['\"]"
        ),
    ),
    (
        "Windows user profile path",
        re.compile(r"(?i)\b[A-Z]:\\Users\\[A-Za-z0-9._-]+(?:\\[^\s'\"`]+)?"),
    ),
    (
        "literal SLURM account",
        re.compile(r"(?m)^\s*#SBATCH\s+(?:-A\s+|--account(?:=|\s+))(?![<{])[A-Za-z0-9_.-]+\s*$"),
    ),
    (
        "literal SLURM partition",
        re.compile(r"(?m)^\s*#SBATCH\s+(?:-p\s+|--partition(?:=|\s+))(?![<{])[A-Za-z0-9_.-]+\s*$"),
    ),
]

for root in private_roots:
    rules.append(
        (
            "site-specific absolute path",
            re.compile(re.escape(root) + r"(?!<)[^\s'\"`]+"),
        )
    )

safe_email_domains = {
    "example.com",
    "example.org",
    "example.net",
    "users.noreply.github.com",
}

bad: list[tuple[str, str]] = []
for path in candidate_paths():
    if not path.is_file():
        continue
    if any(prefix == path or prefix in path.parents for prefix in SKIP_PREFIXES):
        continue
    if path.suffix.lower() in SKIP_SUFFIXES:
        continue
    try:
        raw = path.read_bytes()
    except OSError:
        continue
    if b"\0" in raw:
        continue
    try:
        text = raw.decode("utf-8")
    except UnicodeDecodeError:
        continue

    for rule_name, pattern in rules:
        matches = list(pattern.finditer(text))
        if rule_name == "email address":
            matches = [
                match
                for match in matches
                if match.group(0).rsplit("@", 1)[-1].lower() not in safe_email_domains
            ]
        if matches:
            bad.append((str(path), rule_name))

if bad:
    for path, rule_name in sorted(set(bad)):
        print(f"privacy rule {rule_name!r} matched in {path}")
    sys.exit(1)

print("No private identifiers or credential material found.")
PY
