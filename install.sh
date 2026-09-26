#!/bin/bash
# install.sh — Symlink vibedft skills to your agent's skill directory.
#
# Usage:
#   bash install.sh                    # OpenCode (default)
#   bash install.sh opencode           # OpenCode
#   bash install.sh claude-code        # Claude Code
#   bash install.sh claude             # alias
#   bash install.sh --prefix <dir>     # custom directory
#   bash install.sh uninstall          # remove symlinks
#
# Set SKILLS_DIR to override the target directory regardless of agent type.

set -euo pipefail

REPO_DIR="$(cd "$(dirname "$0")" && pwd)"
SKILLS_SRC="$REPO_DIR/skills"

# ── Target directory ──────────────────────────────────────────────
AGENT="${1:-opencode}"

if [ -n "${SKILLS_DIR:-}" ]; then
    DEST="$SKILLS_DIR"
elif [ "$AGENT" = "uninstall" ]; then
    # uninstall: detect from existing symlinks
    installdir=()
    for guess in \
        "${XDG_CONFIG_HOME:-$HOME/.config}/opencode/skills" \
        "$HOME/.claude/skills"; do
        if [ -d "$guess" ]; then
            for f in "$SKILLS_SRC"/*/; do
                name=$(basename "$f")
                if [ -L "$guess/$name" ] || [ -L "$guess/vibedft-$name" ]; then
                    installdir+=("$guess")
                    break
                fi
            done
        fi
    done
    if [ ${#installdir[@]} -eq 0 ]; then
        echo "No vibedft skills found installed."
        exit 1
    fi
    DEST="${installdir[0]}"
elif [ "$AGENT" = "--prefix" ]; then
    DEST="$2"
elif [ "$AGENT" = "opencode" ]; then
    DEST="${SKILLS_DIR:-${XDG_CONFIG_HOME:-$HOME/.config}/opencode/skills}"
elif [ "$AGENT" = "claude-code" ] || [ "$AGENT" = "claude" ]; then
    DEST="$HOME/.claude/skills"
else
    echo "Unknown agent: $AGENT"
    echo "Usage: $0 [opencode|claude-code|--prefix <dir>|uninstall]"
    exit 1
fi

# ── Uninstall ─────────────────────────────────────────────────────
if [ "$AGENT" = "uninstall" ]; then
    echo "Removing vibedft symlinks from $DEST ..."
    for skill_dir in "$SKILLS_SRC"/*/; do
        name=$(basename "$skill_dir")
        rm -f "$DEST/$name" "$DEST/vibedft-$name"
    done
    echo "Done."
    exit 0
fi

# ── Install ───────────────────────────────────────────────────────
mkdir -p "$DEST"
count=0
for skill_dir in "$SKILLS_SRC"/*/; do
    name=$(basename "$skill_dir")
    # Prefer short name; use vibedft- prefix only as fallback to avoid collision
    if [ -e "$DEST/$name" ] && [ ! -L "$DEST/$name" ]; then
        ln -sfn "$skill_dir" "$DEST/vibedft-$name"
        echo "  symlink $DEST/vibedft-$name  (short name $name exists)"
    else
        ln -sfn "$skill_dir" "$DEST/$name"
        echo "  symlink $DEST/$name"
    fi
    count=$((count + 1))
done

# ── PYTHONPATH hint ───────────────────────────────────────────────
echo ""
echo "Installed $count skills to $DEST"
echo ""
echo "Add these lines to your shell profile:"
echo ""
echo "  export VIBEDFT_HOME=\"$REPO_DIR\""
echo "  export PYTHONPATH=\"\$PYTHONPATH:\$VIBEDFT_HOME/skills\""
echo ""
echo "Then set your runtime variables (see .env.example)."
