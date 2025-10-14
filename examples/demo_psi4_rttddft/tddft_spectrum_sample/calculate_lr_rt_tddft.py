import maxwelllink as mxl

model = mxl.RTTDDFTModel(
    engine="psi4",
    molecule_xyz="../../../tests/data/hcn.xyz",
    functional="b3lyp",
    basis="sto-3g",
    dt_rttddft_au=0.04,
    delta_kick_au=1.0e-3,
    delta_kick_direction="xyz",
    memory="2GB",
    verbose=True,
    remove_permanent_dipole=False,
    dft_grid_name="SG0",
    electron_propagation="etrs",
)

model.initialize(dt_new=0.12, molecule_id=0)

model._get_lr_tddft_spectrum(states=20, tda=False)


model = mxl.RTTDDFTModel(
    engine="psi4",
    molecule_xyz="../../../tests/data/hcn.xyz",
    functional="b3lyp",
    basis="sto-3g",
    dt_rttddft_au=0.04,
    delta_kick_au=1.0e-3,
    delta_kick_direction="xyz",
    memory="2GB",
    verbose=True,
    remove_permanent_dipole=False,
    dft_grid_name="SG0",
    electron_propagation="pc",
)

model.initialize(dt_new=0.12, molecule_id=0)

model._propagate_full_rt_tddft(nsteps=10000)
