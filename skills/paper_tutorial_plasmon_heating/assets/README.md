# Paper Tutorial Assets

This skill bundles lightweight, manifest-backed reproducibility assets for the scoped Figure 5b/5d/5e tutorial:
- MEEP/MaxwellLink template inputs and Slurm launch scripts for EM-only, TLS, and BOMD branches
- HCN reference geometry (`hcn.xyz`) for ASE/Psi4 driver setup
- small archived flux reference files for the no-molecule spectrum
- scoped plotting and acceptance-check script (`fig5bde_postprocess.py`)

Large raw trajectory/result bundles (for example full `mol_0..255_data.npz` families) are intentionally excluded.

Use these assets by copying `assets/implementation_2025/` into `projects/YYYY-MM-DD-<scope>/` (recommended scope: `plasmon-heating`).
