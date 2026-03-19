# MaxwellLink skills

This `skills/` directory is for AI-agent workflows (not end-user runtime code).

Start here:
- `skills/mxl-index/SKILL.md` (router: choose solver/driver/HPC/post-processing workflows)
- `skills/mxl-project-scaffold/SKILL.md` (create `projects/YYYY-MM-DD-NAME/` inputs from templates)
- `skills/mxl-cli/SKILL.md` (use `mxl init`, `mxl clean`, and `mxl hpc set` commands for the light-weighted agent workflow)

Each skill is a self-contained folder with:
- `SKILL.md` (required instructions + YAML frontmatter)
- optional `scripts/`, `references/`, `assets/`

**For custom HPC setting, run `mxl hpc set path/to/HPC_PROFILE.json` (persistent user profile at `~/.maxwelllink/HPC_PROFILE.json`).**
