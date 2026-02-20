---
name: paper_tutorial_plasmon_heating
description: Dual-purpose tutorial to reproduce manuscript Figure 5b, 5d, and 5e plasmonic flux and HCN heating maps (EM-only, TLS, and ASE/BOMD), and to transfer the workflow to related plasmon-molecule systems by adjusting geometry, molecular model, and run procedures.
---

# Paper Tutorial: Plasmonic Flux and HCN Heating (Figure 5b, 5d, 5e)

Use this skill to reproduce the scoped Pt/Si plasmonic manuscript results for:
- Figure 5b: EM-only flux/absorption spectrum.
- Figure 5d: MaxwellLink + TLS heating map (`nmol=256`).
- Figure 5e: MaxwellLink + ASE/Psi4 BOMD heating map (`nmol=256`) plus TLS-vs-BOMD consistency.

## Core Simulation Strategy
- Run all simulations and postprocessing in one dated runtime folder under `projects/YYYY-MM-DD-<scope>/`.
- Recommended scope slug: `plasmon-heating`.
- Copy `assets/implementation_2025/` into the runtime folder, then run three branches with fixed paper anchors:
  - `aa=2.79`, `rr=1.11`, `nmol=256`, y-polarized excitation, coupled propagation `t=60`.
  - EM-only branch (`meep_plasmon_empty`) for Figure 5b.
  - TLS branch (`meep_plasmon_HCN_excitation_tls_strong`) for Figure 5d.
  - BOMD branch (`meep_plasmon_HCN_excitation_bomd_strong`) for Figure 5e.
- Use trend-level agreement targets (peak positions and anisotropy structure), not bitwise trajectory identity.
- Figure 5f (RT-Ehrenfest) is intentionally out of scope.

## Minimal Execution Recipes
Run commands from repository root.

1. Create runtime workspace and stage assets:
```bash
RUN_DATE="${RUN_DATE:-$(date +%F)}"
RUN_SCOPE="${RUN_SCOPE:-plasmon-heating}"
RUN_DIR="projects/${RUN_DATE}-${RUN_SCOPE}"
mkdir -p "$RUN_DIR"
cp -R skills/paper_tutorial_plasmon_heating/assets/implementation_2025 "$RUN_DIR/"
```

2. Figure 5b (EM-only spectrum):
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

3. Figure 5d (TLS coupled map via Slurm templates):
```bash
cd "$RUN_DIR/implementation_2025/meep_plasmon_HCN_excitation_tls_strong"
mkdir -p nmol_256_with_dielectric
cp -R template/* nmol_256_with_dielectric/
cd nmol_256_with_dielectric
sbatch submit_all.sh
```

4. Figure 5e (BOMD coupled map via Slurm templates):
```bash
cd "$RUN_DIR/implementation_2025/meep_plasmon_HCN_excitation_bomd_strong"
mkdir -p nmol_256_with_dielectric
cp -R template/* nmol_256_with_dielectric/
cd nmol_256_with_dielectric
sbatch submit_all.sh
```

5. Build scoped figures and enforce acceptance checks:
```bash
cd "$RUN_DIR/implementation_2025/plotting"
python fig5bde_postprocess.py \
  --mode all \
  --base-dir .. \
  --nmol 256 \
  --figure-out fig5bde_scope.pdf \
  --report-out fig5bde_metrics.json \
  --strict
```

6. Runtime guardrails:
- Do not run simulations inside `skills/paper_tutorial_plasmon_heating/`.
- Do not write runtime output under `skills/`.

## Figure Routing
- `fig_001` (Figure 5b): EM-only Pt/Si plasmonic baseline from vacuum-normalized flux, reproducing a strong absorption peak near `3492 cm^-1` and high absorption around HCN stretch; playbook: `skills/paper_tutorial_plasmon_heating/playbooks/fig_001_plasmon_flux_spectrum.md`.
- `fig_002` (Figure 5d): MaxwellLink+TLS 16x16 molecular-lattice heating map with y-gap-dominant anisotropy and y-edge hotspot structure under the same EM geometry; playbook: `skills/paper_tutorial_plasmon_heating/playbooks/fig_002_tls_heating_map.md`.
- `fig_003` (Figure 5e): MaxwellLink+ASE/Psi4 BOMD heating map with stronger y-gap anisotropy and BOMD mean gain exceeding TLS under matched settings; playbook: `skills/paper_tutorial_plasmon_heating/playbooks/fig_003_bomd_heating_map.md`.

## Beyond Manuscript Exploration
- Sweep `rr` around `1.11` (for example `1.00-1.20`) at fixed `aa=2.79` to quantify peak shifts.
- Sweep perfect-square molecular grids (`nmol=64, 100, 144, 256`) while preserving lattice placement.
- Sweep TLS dipole strength `mu12` around `0.15` and track anisotropy ratio sensitivity.
- Benchmark BOMD settings (`basis`, `memory`, `num_threads`) only after baseline reproduction passes strict checks.
- Keep EM geometry and source window fixed when comparing TLS vs BOMD to preserve interpretability.
