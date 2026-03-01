---
name: paper_tutorial_vsc_water
description: Dual-purpose tutorial to reproduce manuscript liquid-water vibrational strong-coupling spectra and CPU scaling (single-mode cavity plus Meep Bragg resonator), and to adapt the same procedure to related VSC systems by retuning cavity/resonator parameters, size ladders, and execution settings.
---

# Paper Tutorial: Water VSC Reproduction (Figure 3a-3c)

Use this skill for the scoped paper-tutorial reproduction request targeting vibrational strong coupling in water (request alias: Fig.4 scope, manuscript mapping: Figure 3 VSC family).

## Direct Simulation Strategy
- Runtime work must happen in `projects/YYYY-MM-DD-<scope>/`.
- Recommended scope slug for this tutorial: `vsc-water`.
- Copy this skill's assets into the run folder and execute two coupled branches:
  - `implementation_2025/vsc_water_singlemode`: Figure 3a / `fig_001` single-mode cavity spectra.
  - `implementation_2025/vsc_water_maxwellmd`: Figure 3b / `fig_002` Meep-Bragg spectra.
- Reuse timing outputs from both branches to build Figure 3c / `fig_003`.
- Keep the publication parameter anchors:
  - `dt_fs=0.4`, `steps=2.5e5`, cavity marker near `3550 cm^-1`.
  - size ladder `size in {1,2,4,8}`.
  - coupling/current rescaling with system size: `1/size^1.5`.
- Acceptance is trend-level agreement for LP/UP splitting and scaling behavior, not bitwise trajectory identity.

## Minimal Execution Recipes
1. Create a dated run workspace and copy tutorial assets:
```bash
RUN_DIR="projects/YYYY-MM-DD-vsc-water"  # projects/YYYY-MM-DD-<scope>/
mkdir -p "$RUN_DIR"
cp -R skills/paper_tutorial_vsc_water/assets/implementation_2025 "$RUN_DIR/"
```
2. Run single-mode size scan and collect timing stats:
```bash
cd "$RUN_DIR/implementation_2025/vsc_water_singlemode"
bash run_vsc_scan_size.sh
bash get_time_statistics.sh
bash get_lammps_statistics.sh
```
3. Run Meep-Bragg size scan and collect timing stats:
```bash
cd "$RUN_DIR/implementation_2025/vsc_water_maxwellmd"
bash run_vsc_scan_size.sh
bash get_time_statistics.sh
bash get_lammps_statistics.sh
```
4. Build the publication panel (`fig3_vsc.pdf`):
```bash
cd "$RUN_DIR/implementation_2025/plotting"
python plot_fig3_vsc.py
```
5. Artifact gate (all figure families):
```bash
for size in 1 2 4 8; do
  test -f "$RUN_DIR/implementation_2025/vsc_water_singlemode/size_${size}_cube/vsc_water_maxwellmd_data.npz"
  test -f "$RUN_DIR/implementation_2025/vsc_water_maxwellmd/size_${size}_cube/vsc_water_maxwellmd_data.npz"
  test -f "$RUN_DIR/implementation_2025/vsc_water_singlemode/size_${size}_cube/time_count.txt"
  test -f "$RUN_DIR/implementation_2025/vsc_water_singlemode/size_${size}_cube/lammps_count.txt"
  test -f "$RUN_DIR/implementation_2025/vsc_water_maxwellmd/size_${size}_cube/time_count.txt"
  test -f "$RUN_DIR/implementation_2025/vsc_water_maxwellmd/size_${size}_cube/lammps_count.txt"
done
test -f "$RUN_DIR/implementation_2025/plotting/fig3_vsc.pdf"
```

## Figure Routing
- `fig_001`: Reproduce liquid-water VSC in a single-mode cavity, testing LP/UP splitting around `3550 cm^-1` and size-rescaled coupling invariance; playbook: `skills/paper_tutorial_vsc_water/playbooks/fig_001_singlemode_vsc.md`.
- `fig_002`: Reproduce Meep+LAMMPS VSC in a 1D Bragg resonator, including Bragg transmission context and librational/O-H splitting behavior; playbook: `skills/paper_tutorial_vsc_water/playbooks/fig_002_meep_bragg_vsc.md`.
- `fig_003`: Reproduce stepping-time scaling versus driver cores and compare MaxwellLink+LAMMPS against standalone LAMMPS trends; playbook: `skills/paper_tutorial_vsc_water/playbooks/fig_003_scaling.md`.

## Beyond Manuscript Exploration
- Sweep single-mode damping (`damping_au`) while fixing `frequency_au` and baseline size ladder to probe linewidth-splitting tradeoffs.
- Sweep Bragg resonator discretization (`resolution`) and monitor whether LP/UP positions remain stable within tolerance.
- Test shorter production windows (`steps`) for cost reduction, then validate that qualitative splitting order is preserved before claiming equivalence.
- Keep `dt_fs` close to `0.4` unless revalidated.
- Preserve coupling normalization with system size (`1/size^1.5`) for fair scaling comparisons.
- Do not mix timing data from different queue/node policies in one scaling fit.
