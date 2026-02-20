# fig_001 Playbook: 3D Plasmonic Flux Spectrum (Figure 5b)

## Scientific Aim
Reproduce the no-molecule Pt/Si plasmonic absorption spectrum of one periodic unit cell and verify strong absorption near the HCN stretch frequency.

## Figure-Data Mapping (Expected Artifacts)
- `implementation_2025/meep_plasmon_empty/vac/flux0_a2.79.dat`
- `implementation_2025/meep_plasmon_empty/no_mol_with_dielectric/flux_a2.79_r1.11.dat`
- `implementation_2025/plotting/fig5b_metrics.json`
- `implementation_2025/plotting/fig5b_scope.pdf`

## Inputs from Skill Assets
- `assets/implementation_2025/meep_plasmon_empty/template/emitter.py`
- `assets/implementation_2025/meep_plasmon_empty/template/submit_main.sh`
- `assets/implementation_2025/meep_plasmon_empty/template/submit_vac.sh`
- `assets/implementation_2025/meep_plasmon_empty/vac/flux0_a2.79.dat` (archived small reference)
- `assets/implementation_2025/meep_plasmon_empty/no_mol_with_dielectric/flux_a2.79_r1.11.dat` (archived small reference)
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

2. Run vacuum and dielectric cases:
```bash
cd "$RUN_DIR/implementation_2025/meep_plasmon_empty"
mkdir -p vac no_mol_with_dielectric
cp template/* vac/
cp template/* no_mol_with_dielectric/

cd vac
mpirun -np 128 python -u emitter.py -empty -aa 2.79 > flux0_a2.79.out
grep flux1: flux0_a2.79.out | cut -d , -f2- > flux0_a2.79.dat

cd ../no_mol_with_dielectric
mpirun -np 128 python -u emitter.py -dielectric -aa 2.79 -rr 1.11 -nmol 1 > flux_a2.79_r1.11.out
grep flux1: flux_a2.79_r1.11.out | cut -d , -f2- > flux_a2.79_r1.11.dat
```

3. Evaluate acceptance gates and build scoped panel:
```bash
cd "$RUN_DIR/implementation_2025/plotting"
python fig5bde_postprocess.py \
  --mode fig5b \
  --base-dir .. \
  --report-out fig5b_metrics.json \
  --figure-out fig5b_scope.pdf \
  --strict
```

## Acceptance Checks
- Both flux files exist and share a 200-point frequency grid.
- Absorption peak lies in `3450-3535 cm^-1`.
- Maximum absorption is `> 0.90`.
- Absorption at `3465.5930 cm^-1` is `> 0.70`.
