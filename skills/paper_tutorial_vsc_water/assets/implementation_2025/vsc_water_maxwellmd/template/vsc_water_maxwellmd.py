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


# 1. Define the geometry for a 1D bragg resonator
resolution = 25
time_units_fs = 20
rescaling = 0.47

pml_thickness = 2.0 * rescaling
t1 = 0.125 * rescaling
t2 = 0.25 * rescaling
n1 = 2.0
n2 = 1.0

nlayer = 5
layer_indexes = np.array([n2, n1] * nlayer + [1.0] + [n1, n2] * nlayer)
layer_thicknesses = np.array([t2, t1] * nlayer + [0.5 * rescaling] + [t1, t2] * nlayer)

layer_thicknesses[0] += pml_thickness
layer_thicknesses[-1] += pml_thickness
length = np.sum(layer_thicknesses)
layer_centers = np.cumsum(layer_thicknesses) - layer_thicknesses/2
layer_centers = layer_centers - length/2
cell_size = mp.Vector3(length, 0, 0)
pml_layers = [mp.PML(thickness=pml_thickness)]

geometry = [mp.Block(mp.Vector3(layer_thicknesses[i], mp.inf, mp.inf),
    center=mp.Vector3(layer_centers[i], 0, 0), material=mp.Medium(index=layer_indexes[i]))
    for i in range(layer_thicknesses.size)]

enlarge_factor = 1.0
# assume the water box is enlarged by enlarge_factor in each dimension
rescaling_factor = 1e5 / enlarge_factor**1.5

# 2. Perform simulation
molecule = mxl.Molecule(
    hub=hub,
    center=mp.Vector3(0, 0, 0), 
    size=mp.Vector3(0.23, 0, 0), 
    sigma=0.04, 
    dimensions=1, 
    rescaling_factor=rescaling_factor
)

sim = mxl.MeepSimulation(
    hub=hub,
    molecules=[molecule],
    cell_size=cell_size,
    resolution=resolution,
    time_units_fs=time_units_fs,
    geometry=geometry,
    boundary_layers=pml_layers
)

# 3. Run for 100 ps
sim.run(steps=2.5e5)

# 4. Obtain necessary data for post-processing
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
