# fig_002 Playbook: Maxwell-MD VSC in 1D Bragg Resonator

## Scientific Aim
Reproduce strong coupling in a realistic Meep 1D Bragg resonator and capture both O-H polaritonic splitting and librational-band splitting trends.

## Figure-Data Mapping (Expected Artifacts)
- `implementation_2025/vsc_water_maxwellmd/size_1_cube/vsc_water_maxwellmd.py`
- `implementation_2025/vsc_water_maxwellmd/size_2_cube/vsc_water_maxwellmd.py`
- `implementation_2025/vsc_water_maxwellmd/size_4_cube/vsc_water_maxwellmd.py`
- `implementation_2025/vsc_water_maxwellmd/size_8_cube/vsc_water_maxwellmd.py`
- `implementation_2025/vsc_water_maxwellmd/size_1_cube/vsc_water_maxwellmd_data.npz`
- `implementation_2025/vsc_water_maxwellmd/size_2_cube/vsc_water_maxwellmd_data.npz`
- `implementation_2025/vsc_water_maxwellmd/size_4_cube/vsc_water_maxwellmd_data.npz`
- `implementation_2025/vsc_water_maxwellmd/size_8_cube/vsc_water_maxwellmd_data.npz`
- `implementation_2025/vsc_water_maxwellmd/bragg_resonator_spectrum/bragg_resonator_transmission.txt`
- `implementation_2025/plotting/fig3_vsc.pdf` (generated after plotting step)

## Inputs from Skill Assets
- `assets/implementation_2025/vsc_water_maxwellmd/template/vsc_water_maxwellmd.py`
- `assets/implementation_2025/vsc_water_maxwellmd/lmp_input/in_mxl.lmp`
- `assets/implementation_2025/vsc_water_maxwellmd/lmp_input/data.lmp`
- `assets/implementation_2025/vsc_water_maxwellmd/run_vsc_scan_size.sh`
- `assets/implementation_2025/vsc_water_maxwellmd/bragg_resonator_spectrum/bragg_resonator.py`
- `assets/implementation_2025/vsc_water_maxwellmd/bragg_resonator_spectrum/bragg_resonator_transmission.txt`

## Runtime Procedure
1. Stage runtime copy under `projects/`:
```bash
RUN_DIR="projects/YYYY-MM-DD-vsc-water"  # projects/YYYY-MM-DD-<scope>/
mkdir -p "$RUN_DIR"
cp -R skills/paper_tutorial_vsc_water/assets/implementation_2025 "$RUN_DIR/"
```
2. Submit the MaxwellMD size scan:
```bash
cd "$RUN_DIR/implementation_2025/vsc_water_maxwellmd"
bash run_vsc_scan_size.sh
```
3. Extract timing metadata after completion:
```bash
bash get_time_statistics.sh
bash get_lammps_statistics.sh
```
4. Optional Bragg baseline regeneration:
```bash
cd bragg_resonator_spectrum
python bragg_resonator.py
```
5. Confirm expected Maxwell-MD artifacts exist:
```bash
for size in 1 2 4 8; do
  test -f "$RUN_DIR/implementation_2025/vsc_water_maxwellmd/size_${size}_cube/vsc_water_maxwellmd.py"
  test -f "$RUN_DIR/implementation_2025/vsc_water_maxwellmd/size_${size}_cube/vsc_water_maxwellmd_data.npz"
done
test -f "$RUN_DIR/implementation_2025/vsc_water_maxwellmd/bragg_resonator_spectrum/bragg_resonator_transmission.txt"
```

## Acceptance Checks
- Bragg transmission context includes cavity mode near `3550 cm^-1` and water-band markers near `700/1650/3550 cm^-1`.
- Baseline `size=1` with `gamma=1e5` trend yields O-H splitting near `~600 cm^-1`.
- Librational splitting around `~700 cm^-1` appears in coupled spectra.
- `size=2,4,8` runs preserve approximate splitting under `gamma/size^1.5` scaling.
- Quantitative helper (optional): run `implementation_2025/plotting/plot_fig3_vsc.py` and inspect panel (b) + inset for Bragg resonance placement and librational splitting.
