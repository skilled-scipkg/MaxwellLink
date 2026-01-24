# Guidelines for using MaxwellLink 

MaxwellLink is a modular, open-source Python framework for self-consistent light-matter simulations. The package utilizes a socket interface to couple classical EM solvers with a wide range of external molecular drivers.

## Preparing input files 

- Once being asked to prepare input files for using MaxwellLink, go to `projects/` and create a subfolder `YEAR-MM-DD-NAME/` with the date as today and an appropriate `NAME` matching the simulation goal. Then, add simulation input files in this subfolder. 

- Always read `skills/` first to examine whether the proposed simulation by the user is supported by the existing skills in MaxwellLink. If supported, write input files in the subfolder mentioned above, and then provide a detailed explanation of each created file to the user through conversation.

- If you feel confused, also read `docs/source/` for the documentation of MaxwellLink as well as the source code at `src/maxwelllink/`.

## Performing simulations on local machines

- Once being asked to directly perform MaxwellLink simulations, first generate the proper input files with your maximal efforts following ## Preparing input files. Then, perform an independent code review of the generated input files using knowledge from `skills/`. Finally, if `sbatch` command is available, add a SLRUM bash script to submit this simulation within the local subfolder of the input files; otherwise directly perform simulations locally.

## Performing simulations on HPC systems

If `sbatch` command is available, add a SLRUM bash script to submit this simulation within the local subfolder of the input files, according to the HPC setting given at `skills/mxl-hpc-slurm/resources/hpc_setting.md`. Then, submit the SLURM jobs.

## Postprocessing 

If you are asked to post-process the simulation data, wait the simulation if finished either on local machines or HPC SLURM, and then create Python plotting scripts accordingly to provide the visuallization file directly.

## Debugging on failed or unsuccessful simulations

If the simulation cannot be finished due to any bug or is finished but generating undesired (apparently wrong) results, write a summary file ``debugging.md`` to briefly conclude what you have done in this simulation as well as the possible causes of the bug. Then, create a new, different subfolder `YEAR-MM-DD-NAME/` within `projects/`, learning what you have done incorrectly, and redo the entire simulation procedure, from ## Preparing input files, ## Performing simulations, to ## Postprocessing.

## Other requests from the users

- For other requests from the users, follow the default guidelines and serve the user with your best efforts.

## LLM-to-simulation research workflow

- When asked to evaluate or publish on NL-to-MaxwellLink generation, load `skills/mxl-nl2sim-eval/SKILL.md`.
- Keep benchmark intents under `benchmarks/nl2sim/<task>/` with `intent.md`, `gold/config.json`, `gold/notes.md`, and `checks.yaml` (physics rules).
- Run loop: scaffold with `skills/mxl-project-scaffold/scripts/mxl_scaffold_project.py`, fill `config.json`, run `skills/mxl-project-scaffold/scripts/mxl_validate_project.py`, execute a short sanity simulation, apply physics validators, then self-critique and retry.
- Highlight research levers: harder benchmark coverage, tighter physics validators (stability/units/geometry/charge), and metrics (success rate, retries, runtime, spectral/dipole agreement) compared to human baselines.
