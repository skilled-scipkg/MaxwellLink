# fig_003 Playbook: HCN Heating Map with ASE/BOMD Drivers (Figure 5e)

## Scientific Aim
Reproduce the HCN energy-gain map using MaxwellLink-coupled ASE/Psi4 BOMD drivers and verify stronger mean heating than TLS under identical EM geometry.

## Figure-Data Mapping (Expected Artifacts)
- `implementation_2025/meep_plasmon_HCN_excitation_bomd_strong/nmol_256_with_dielectric/mol_0_data.npz`
- `implementation_2025/meep_plasmon_HCN_excitation_bomd_strong/nmol_256_with_dielectric/mol_255_data.npz`
- `implementation_2025/meep_plasmon_HCN_excitation_bomd_strong/nmol_256_with_dielectric/flux_a2.79_r1.11.out`
- `implementation_2025/plotting/fig5e_bomd_metrics.json`
- `implementation_2025/plotting/fig5e_bomd_scope.pdf`
- `implementation_2025/plotting/fig5bde_metrics.json` (cross-panel strict report)

## Inputs from Skill Assets
- `assets/implementation_2025/meep_plasmon_HCN_excitation_bomd_strong/template/emitter.py`
- `assets/implementation_2025/meep_plasmon_HCN_excitation_bomd_strong/template/submit_all.sh`
- `assets/implementation_2025/meep_plasmon_HCN_excitation_bomd_strong/template/submit_main.sh`
- `assets/implementation_2025/meep_plasmon_HCN_excitation_bomd_strong/template/submit_driver.sh`
- `assets/implementation_2025/meep_plasmon_HCN_excitation_bomd_strong/HCN_benchmark/hcn.xyz`
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

2. Prepare and submit BOMD coupled runs (Slurm):
```bash
cd "$RUN_DIR/implementation_2025/meep_plasmon_HCN_excitation_bomd_strong"
mkdir -p nmol_256_with_dielectric
cp -R template/* nmol_256_with_dielectric/
cd nmol_256_with_dielectric
sbatch submit_all.sh
```

3. Evaluate BOMD map metrics and scoped panel:
```bash
cd "$RUN_DIR/implementation_2025/plotting"
python fig5bde_postprocess.py \
  --mode fig5e \
  --base-dir .. \
  --nmol 256 \
  --report-out fig5e_bomd_metrics.json \
  --figure-out fig5e_bomd_scope.pdf \
  --strict
```

4. Run cross-panel strict checks (requires TLS + BOMD branches complete):
```bash
python fig5bde_postprocess.py \
  --mode all \
  --base-dir .. \
  --nmol 256 \
  --report-out fig5bde_metrics.json \
  --figure-out fig5bde_scope.pdf \
  --strict
```

## Acceptance Checks
- Full molecular output family exists (`mol_0_data.npz` through `mol_255_data.npz`).
- Directional anisotropy passes: `mean(y-gap) / mean(x-gap) >= 4.0`.
- Hotspot location check passes: dominant hotspot is y-edge aligned (`|y| >= 0.8 um`) and near x-center (`|x| <= 0.35 um`).
- Cross-panel consistency passes: BOMD mean gain exceeds TLS mean gain for matched settings.
