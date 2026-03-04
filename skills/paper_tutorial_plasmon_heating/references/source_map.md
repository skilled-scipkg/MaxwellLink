# Source Map: paper_tutorial_plasmon_heating

## Simulation Drivers
- `maxwelllink.MeepSimulation`: coupled EM + molecule propagation in TLS/BOMD branches.
- `maxwelllink.Molecule`: per-site molecular proxy objects placed on the 16x16 lattice.
- `maxwelllink.SocketHub`: socket transport between MEEP main process and molecular drivers.
- `meep.Simulation`: EM-only baseline propagation for vacuum and dielectric flux references.

## Figure-Scoped Input Templates
- EM-only baseline:
  - `assets/implementation_2025/meep_plasmon_empty/template/emitter.py`
  - `assets/implementation_2025/meep_plasmon_empty/template/submit_main.sh`
  - `assets/implementation_2025/meep_plasmon_empty/template/submit_vac.sh`
- TLS map:
  - `assets/implementation_2025/meep_plasmon_HCN_excitation_tls_strong/template/emitter.py`
  - `assets/implementation_2025/meep_plasmon_HCN_excitation_tls_strong/template/submit_main.sh`
  - `assets/implementation_2025/meep_plasmon_HCN_excitation_tls_strong/template/submit_driver.sh`
- BOMD map:
  - `assets/implementation_2025/meep_plasmon_HCN_excitation_bomd_strong/template/emitter.py`
  - `assets/implementation_2025/meep_plasmon_HCN_excitation_bomd_strong/template/submit_main.sh`
  - `assets/implementation_2025/meep_plasmon_HCN_excitation_bomd_strong/template/submit_driver.sh`
  - `assets/implementation_2025/meep_plasmon_HCN_excitation_bomd_strong/HCN_benchmark/hcn.xyz`

## Scoped Postprocessing
- `assets/implementation_2025/plotting/fig5bde_postprocess.py` computes:
  - Figure 5b absorption from `flux0_a2.79.dat` and `flux_a2.79_r1.11.dat`.
  - Figure 5d/5e energy-gain maps from `mol_*.npz` families.
  - anisotropy checks, hotspot-position checks, and TLS-vs-BOMD cross-panel metrics.
