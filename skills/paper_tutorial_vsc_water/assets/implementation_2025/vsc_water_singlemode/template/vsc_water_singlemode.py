import numpy as np
import maxwelllink as mxl
from maxwelllink import sockets as mxs
import meep as mp

host, port = mxs.get_available_host_port(localhost=False)
hub = mxl.SocketHub(host=host, port=port, timeout=10.0, latency=1e-3)

print(f"SocketHub listening on {host}:{port}")

# save host and port number to a file so mxl_driver can read it
with open("tcp_host_port_info.txt", "w") as f:
    f.write(f"{host}\n{port}\n")

# rescaling factor to take into account the TIP4P water model 
molecule = mxl.Molecule(
    hub=hub,
    rescaling_factor=0.73,
)

enlarge_factor = 1.0

au_to_cminverse = 219474.63
frequency_au = 3550 / au_to_cminverse
coupling_strength = 4e-4 / enlarge_factor**1.5
print(f"Coupling strength: {coupling_strength:.3e} au")
damping_au = 0e-4
fs_to_au = 41.3413745758
dt_fs = 0.4
dt_au = dt_fs * fs_to_au

sim = mxl.SingleModeSimulation(
    hub=hub,
    molecules=[molecule],
    frequency_au=frequency_au,
    coupling_strength=coupling_strength,
    damping_au=damping_au,
    coupling_axis="z",
    drive=0.0,
    dt_au=dt_au,
    record_history=True,
    include_dse=True,
)

# Run for 100 ps
sim.run(steps=2.5e5)

# Obtain necessary data for post-processing
if mp.am_master():
    print("Simulation completed. Collecting data...")
    mux_au = np.array([ad["mux_au"] for ad in molecule.additional_data_history])
    muy_au = np.array([ad["muy_au"] for ad in molecule.additional_data_history])
    muz_au = np.array([ad["muz_au"] for ad in molecule.additional_data_history])
    time_au = np.array([ad["time_au"] for ad in molecule.additional_data_history])
    energy_molecule_au = np.array([ad["energy_au"] for ad in molecule.additional_data_history])
    np.savez(
        "vsc_water_maxwellmd_data.npz",
        mux_au=mux_au,
        muy_au=muy_au,
        muz_au=muz_au,
        time_au=time_au,
        energy_molecule_au=energy_molecule_au
    )
