---
name: mxl-cli
description: This skill should be used when users need to initialize, clean, or configure MaxwellLink workspaces using the mxl command family.
---

# MaxwellLink CLI workflow

## When to use
- Use this skill for questions or tasks about `mxl` commands:
  - `mxl init` / `mxl-init`
  - `mxl clean` / `mxl-clean`
  - `mxl hpc set <file>`
- Use this skill before running simulations when the workspace is not prepared yet.

## Command behaviors
- `mxl init`:
  - Copies `AGENTS.md` into the current directory.
  - Creates symlinks for `src/`, `tests/`, `skills/`, `docs/`, `media/`, `tutorials/`, and `README.md`.
  - Creates `CLAUDE.md` and `GEMINI.md` symlinks pointing to local `AGENTS.md`.
  - If `sbatch` is available, ensures `~/.maxwelllink/HPC_PROFILE.json` exists (seeded from package default if needed).
  - If `sbatch` is available, creates local `HPC_PROFILE.json` symlink to `~/.maxwelllink/HPC_PROFILE.json`.
- `mxl clean`:
  - Removes local artifacts created by `mxl init` from the current directory.
  - Leaves `~/.maxwelllink/HPC_PROFILE.json` in place.
- `mxl hpc set <file>`:
  - Validates and installs a user profile JSON to `~/.maxwelllink/HPC_PROFILE.json`.
  - Expected keys: `slurm_default_partition`, `slurm_defaults`, `slurm_resource_policy`.

## Conflict handling
- If a target path already exists and differs from the managed target/content:
  - `mxl init` and `mxl clean` fail with a conflict error.
  - Re-run with `--force` to overwrite/remove conflicting paths.
- Re-running `mxl init` or `mxl clean` without conflicts is idempotent.

## References
- CLI implementation: `src/maxwelllink/cli/mxl.py`
- Init implementation: `src/maxwelllink/cli/mxl_init.py`
- Clean implementation: `src/maxwelllink/cli/mxl_clean.py`
- HPC profile command: `src/maxwelllink/cli/mxl_hpc.py`
