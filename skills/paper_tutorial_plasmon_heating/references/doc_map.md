# Doc Map: paper_tutorial_plasmon_heating

## Scope
- Reproduction track: 3D plasmonic flux baseline and HCN heating maps for Figure 5b, 5d, and 5e.
- Figure ids covered: `fig_001`, `fig_002`, `fig_003`.
- Runtime convention: run all calculations in `projects/YYYY-MM-DD-<scope>/` (recommended scope slug: `plasmon-heating`).

## Playbook Index
- `fig_001` (EM-only Pt/Si flux spectrum baseline with HCN-window absorption checks): `skills/paper_tutorial_plasmon_heating/playbooks/fig_001_plasmon_flux_spectrum.md`
- `fig_002` (MaxwellLink + TLS 2D heating map with anisotropy/hotspot checks): `skills/paper_tutorial_plasmon_heating/playbooks/fig_002_tls_heating_map.md`
- `fig_003` (MaxwellLink + ASE/BOMD 2D heating map with TLS cross-panel checks): `skills/paper_tutorial_plasmon_heating/playbooks/fig_003_bomd_heating_map.md`

## Asset Anchors
- EM-only branch: `skills/paper_tutorial_plasmon_heating/assets/implementation_2025/meep_plasmon_empty`
- TLS branch: `skills/paper_tutorial_plasmon_heating/assets/implementation_2025/meep_plasmon_HCN_excitation_tls_strong`
- BOMD branch: `skills/paper_tutorial_plasmon_heating/assets/implementation_2025/meep_plasmon_HCN_excitation_bomd_strong`
- Scoped plotting/checking: `skills/paper_tutorial_plasmon_heating/assets/implementation_2025/plotting`
