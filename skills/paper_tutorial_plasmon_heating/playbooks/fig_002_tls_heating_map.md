# fig_002 Playbook: HCN Heating Map with TLS Drivers (Figure 5d)

## Scientific Aim
Reproduce the spatially inhomogeneous HCN energy-gain map using MaxwellLink-coupled TLS molecular drivers on a 16x16 lattice (`nmol=256`) in the plasmonic unit cell.

## Figure-Data Mapping (Expected Artifacts)
- `implementation_2025/meep_plasmon_HCN_excitation_tls_strong/nmol_256_with_dielectric/mol_0_data.npz`
- `implementation_2025/meep_plasmon_HCN_excitation_tls_strong/nmol_256_with_dielectric/mol_255_data.npz`
- `implementation_2025/meep_plasmon_HCN_excitation_tls_strong/nmol_256_with_dielectric/flux_a2.79_r1.11.out`
- `implementation_2025/plotting/fig5d_tls_metrics.json`
- `implementation_2025/plotting/fig5d_tls_scope.pdf`

## Inputs from Skill Assets
- `assets/implementation_2025/meep_plasmon_HCN_excitation_tls_strong/template/emitter.py`
- `assets/implementation_2025/meep_plasmon_HCN_excitation_tls_strong/template/submit_all.sh`
- `assets/implementation_2025/meep_plasmon_HCN_excitation_tls_strong/template/submit_main.sh`
- `assets/implementation_2025/meep_plasmon_HCN_excitation_tls_strong/template/submit_driver.sh`
- `assets/implementation_2025/plotting/fig5bde_postprocess.py`

## Runtime Procedure
1. Stage runtime copy under `projects/YYYY-MM-DD-<scope>/`:
```bash
RUN_DATE="${RUN_DATE:-$(date +%F)}"
RUN_SCOPE="${RUN_SCOPE:-plasmon-heating}"
RUN_DIR="projects/${RUN_DATE}-${RUN_SCOPE}"
mkdir -p "$RUN_DIR"
cp -R skills/paper_tutorial_plasmon_heating/assets/implementation_2025 "$RUN_DIR/"
```

2. Prepare and submit TLS coupled runs (Slurm):
```bash
cd "$RUN_DIR/implementation_2025/meep_plasmon_HCN_excitation_tls_strong"
mkdir -p nmol_256_with_dielectric
cp -R template/* nmol_256_with_dielectric/
cd nmol_256_with_dielectric
sbatch submit_all.sh
```

3. After all jobs complete, evaluate map metrics and generate scoped panel:
```bash
cd "$RUN_DIR/implementation_2025/plotting"
python fig5bde_postprocess.py \
  --mode fig5d \
  --base-dir .. \
  --nmol 256 \
  --report-out fig5d_tls_metrics.json \
  --figure-out fig5d_tls_scope.pdf \
  --strict
```

## Acceptance Checks
- Full molecular output family exists (`mol_0_data.npz` through `mol_255_data.npz`).
- Per-molecule baseline-shifted mean gains are nonnegative.
- Directional anisotropy passes: `mean(y-gap) / mean(x-gap) >= 3.0`.
- Hotspot location check passes: dominant hotspot is y-edge aligned (`|y| >= 0.8 um`) and near x-center (`|x| <= 0.35 um`).
