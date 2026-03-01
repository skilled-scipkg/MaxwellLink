# Source Map: paper_tutorial_vsc_water

## Simulation Drivers
- `maxwelllink.SingleModeSimulation`: used by `skills/paper_tutorial_vsc_water/assets/implementation_2025/vsc_water_singlemode/template/vsc_water_singlemode.py`.
- `maxwelllink.MeepSimulation`: used by `skills/paper_tutorial_vsc_water/assets/implementation_2025/vsc_water_maxwellmd/template/vsc_water_maxwellmd.py`.
- `maxwelllink.SocketHub`: socket transport between MaxwellLink and LAMMPS in both branches.

## LAMMPS Coupling Inputs
- `skills/paper_tutorial_vsc_water/assets/implementation_2025/vsc_water_singlemode/lmp_input/in_mxl.lmp`
- `skills/paper_tutorial_vsc_water/assets/implementation_2025/vsc_water_maxwellmd/lmp_input/in_mxl.lmp`

Both templates expose `HOST`/`PORT` replacement points used by `submit_driver.sh` scripts.

## Figure Assembly
- `skills/paper_tutorial_vsc_water/assets/implementation_2025/plotting/plot_fig3_vsc.py` combines:
  - single-mode spectra,
  - Meep-Bragg spectra + transmission inset,
  - timing scaling curves.
- `skills/paper_tutorial_vsc_water/assets/implementation_2025/plotting/columnplots.py` provides plotting utilities used by `plot_fig3_vsc.py`.
