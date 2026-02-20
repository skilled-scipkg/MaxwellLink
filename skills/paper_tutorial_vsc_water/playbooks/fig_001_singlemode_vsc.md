# fig_001 Playbook: Single-Mode Cavity VSC IR Spectra

## Scientific Aim
Reproduce liquid-water vibrational strong coupling in a lossless single-mode cavity and verify that LP/UP splitting stays approximately size-invariant when coupling is scaled as `1/size^1.5`.

## Figure-Data Mapping (Expected Artifacts)
- `implementation_2025/vsc_water_singlemode/size_1_cube/vsc_water_singlemode.py`
- `implementation_2025/vsc_water_singlemode/size_2_cube/vsc_water_singlemode.py`
- `implementation_2025/vsc_water_singlemode/size_4_cube/vsc_water_singlemode.py`
- `implementation_2025/vsc_water_singlemode/size_8_cube/vsc_water_singlemode.py`
- `implementation_2025/vsc_water_singlemode/size_1_cube/vsc_water_maxwellmd_data.npz`
- `implementation_2025/vsc_water_singlemode/size_2_cube/vsc_water_maxwellmd_data.npz`
- `implementation_2025/vsc_water_singlemode/size_4_cube/vsc_water_maxwellmd_data.npz`
- `implementation_2025/vsc_water_singlemode/size_8_cube/vsc_water_maxwellmd_data.npz`
- `implementation_2025/plotting/fig3_vsc.pdf` (generated after plotting step)

## Inputs from Skill Assets
- `assets/implementation_2025/vsc_water_singlemode/template/vsc_water_singlemode.py`
- `assets/implementation_2025/vsc_water_singlemode/lmp_input/in_mxl.lmp`
- `assets/implementation_2025/vsc_water_singlemode/lmp_input/data.lmp`
- `assets/implementation_2025/vsc_water_singlemode/run_vsc_scan_size.sh`
- `assets/implementation_2025/vsc_water_singlemode/get_time_statistics.sh`
- `assets/implementation_2025/vsc_water_singlemode/get_lammps_statistics.sh`

## Runtime Procedure
1. Stage runtime copy under `projects/`:
```bash
RUN_DIR="projects/YYYY-MM-DD-vsc-water"  # projects/YYYY-MM-DD-<scope>/
mkdir -p "$RUN_DIR"
cp -R skills/paper_tutorial_vsc_water/assets/implementation_2025 "$RUN_DIR/"
```
2. Submit the size scan:
```bash
cd "$RUN_DIR/implementation_2025/vsc_water_singlemode"
bash run_vsc_scan_size.sh
```
3. After jobs finish, collect stepping/lammps statistics:
```bash
bash get_time_statistics.sh
bash get_lammps_statistics.sh
```
4. Confirm expected single-mode artifacts exist:
```bash
for size in 1 2 4 8; do
  test -f "$RUN_DIR/implementation_2025/vsc_water_singlemode/size_${size}_cube/vsc_water_singlemode.py"
  test -f "$RUN_DIR/implementation_2025/vsc_water_singlemode/size_${size}_cube/vsc_water_maxwellmd_data.npz"
done
```

## Acceptance Checks
- Cavity-on spectra show LP/UP splitting near the `3550 cm^-1` marker.
- Baseline size (`size=1`, equivalent to `N_H2O=216`) gives O-H splitting near `~800 cm^-1` trend target.
- `size=2,4,8` runs remain approximately split-invariant under `1/size^1.5` scaling.
- Cavity-off trace remains qualitatively distinct from cavity-coupled traces.
- Quantitative helper (optional): run `implementation_2025/plotting/plot_fig3_vsc.py` and inspect panel (a) for clear LP/UP peaks and size-invariant gap ordering.
