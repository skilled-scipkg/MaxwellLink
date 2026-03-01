# fig_003 Playbook: Stepping-Time Scaling Versus Driver CPU Cores

## Scientific Aim
Reproduce log-log scaling of stepping time against driver cores (`1,8,64,512`) and compare MaxwellLink+LAMMPS overhead to standalone LAMMPS.

## Figure-Data Mapping (Expected Artifacts)
- `implementation_2025/vsc_water_singlemode/size_1_cube/time_count.txt`
- `implementation_2025/vsc_water_singlemode/size_2_cube/time_count.txt`
- `implementation_2025/vsc_water_singlemode/size_4_cube/time_count.txt`
- `implementation_2025/vsc_water_singlemode/size_8_cube/time_count.txt`
- `implementation_2025/vsc_water_singlemode/size_1_cube/lammps_count.txt`
- `implementation_2025/vsc_water_singlemode/size_2_cube/lammps_count.txt`
- `implementation_2025/vsc_water_singlemode/size_4_cube/lammps_count.txt`
- `implementation_2025/vsc_water_singlemode/size_8_cube/lammps_count.txt`
- `implementation_2025/vsc_water_maxwellmd/size_1_cube/time_count.txt`
- `implementation_2025/vsc_water_maxwellmd/size_2_cube/time_count.txt`
- `implementation_2025/vsc_water_maxwellmd/size_4_cube/time_count.txt`
- `implementation_2025/vsc_water_maxwellmd/size_8_cube/time_count.txt`
- `implementation_2025/vsc_water_maxwellmd/size_1_cube/lammps_count.txt`
- `implementation_2025/vsc_water_maxwellmd/size_2_cube/lammps_count.txt`
- `implementation_2025/vsc_water_maxwellmd/size_4_cube/lammps_count.txt`
- `implementation_2025/vsc_water_maxwellmd/size_8_cube/lammps_count.txt`
- `implementation_2025/plotting/fig3_vsc.pdf`

## Inputs from Skill Assets
- `assets/implementation_2025/vsc_water_singlemode/get_time_statistics.sh`
- `assets/implementation_2025/vsc_water_singlemode/get_lammps_statistics.sh`
- `assets/implementation_2025/vsc_water_maxwellmd/get_time_statistics.sh`
- `assets/implementation_2025/vsc_water_maxwellmd/get_lammps_statistics.sh`
- `assets/implementation_2025/plotting/plot_fig3_vsc.py`
- `assets/implementation_2025/plotting/columnplots.py`

## Runtime Procedure
1. Ensure both `vsc_water_singlemode/size_*` and `vsc_water_maxwellmd/size_*` runs are completed (from `fig_001` and `fig_002` workflows).
2. Refresh timing summaries in both branches:
```bash
RUN_DIR="projects/YYYY-MM-DD-vsc-water"  # projects/YYYY-MM-DD-<scope>/

cd "$RUN_DIR/implementation_2025/vsc_water_singlemode"
bash get_time_statistics.sh
bash get_lammps_statistics.sh

cd ../vsc_water_maxwellmd
bash get_time_statistics.sh
bash get_lammps_statistics.sh
```
3. Build the panel figure (including scaling panel):
```bash
cd ../plotting
python plot_fig3_vsc.py
```
4. Confirm scaling artifacts are present:
```bash
for size in 1 2 4 8; do
  test -f "$RUN_DIR/implementation_2025/vsc_water_singlemode/size_${size}_cube/time_count.txt"
  test -f "$RUN_DIR/implementation_2025/vsc_water_singlemode/size_${size}_cube/lammps_count.txt"
  test -f "$RUN_DIR/implementation_2025/vsc_water_maxwellmd/size_${size}_cube/time_count.txt"
  test -f "$RUN_DIR/implementation_2025/vsc_water_maxwellmd/size_${size}_cube/lammps_count.txt"
done
test -f "$RUN_DIR/implementation_2025/plotting/fig3_vsc.pdf"
```

## Acceptance Checks
- X-axis core ladder appears as `1, 8, 64, 512` (log scale).
- MaxwellLink+LAMMPS curve stays within the same order of magnitude as standalone LAMMPS across ladder points.
- No pathological superlinear jump is observed in MaxwellLink timings.
- Trend supports manuscript statement that coupling introduces only modest additional stepping cost.
