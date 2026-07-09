---
name: runtime-config
description: Generic runtime configuration guide for vibedft. Use when preparing remote execution, scheduler settings, VASP command templates, or POTCAR access without storing site-specific details in the repository.
---

# Runtime Config Guide

## Purpose

This repository intentionally omits tracked site-specific runtime facts.

Before generating remote submit scripts, the agent should look for local
runtime notes that are not tracked by git. Recommended locations:

1. `.local/runtime-profile.md`
2. `.local/runtime-profile.json`
3. shell environment variables such as `VIBEDFT_*`

If those files do not exist, the agent should ask for only the missing values
needed for the current task.

## Tracked/Local Boundary

Tracked files in this repo should stay generic.

Local untracked notes may contain the runtime details required to adapt the
workflow to a specific machine, scheduler, or software stack.

## Never Commit

Do not commit any of the following:

- hostnames or IP addresses
- usernames
- SSH key paths
- scheduler partition or QoS values tied to a specific site
- core counts tied to a specific site
- software installation paths
- POTCAR library paths
- site-specific workdir conventions

## Private Runtime Notes Template

Recommended sections for `.local/runtime-profile.md`:

```md
# Runtime Profile

## Connection
- Remote access method
- Login command pattern
- Working directory rule

## Scheduler
- Required scheduler directives
- Default walltime policy
- Queue naming rules

## VASP Command
- Command template with `{vasp_bin}` placeholder

## POTCAR
- Library location
- Any wrapper or module requirement

## Validation
- One known-good dry-run command
```

## Environment Variables

Public repo code reads only generic variables:

- `VIBEDFT_VASP_CMD`
- `VIBEDFT_POTCAR_DIR`

Optional scheduler values should live in each task config JSON, not in tracked
repository defaults.
