import numpy as np
import meep as mp

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

sourceLocation = mp.Vector3(layer_centers[0] - layer_thicknesses[0]/4, 0, 0)
transmissionMonitorLocation = mp.Vector3(layer_centers[-1]-pml_thickness/2, 0, 0)
reflectionMonitorLocation = mp.Vector3(layer_centers[0] + layer_thicknesses[0]/4, 0, 0)

frequency = 1.0 / rescaling  # target frequency in meep units
frequencyWidth = 0.5 / rescaling
numberFrequencies = 200
sources = [mp.Source(
    mp.GaussianSource(frequency=frequency,fwidth=frequencyWidth),
    component=mp.Ez,
    center=sourceLocation,
    size=mp.Vector3(0, 0, 0)
    )]

# Perform simulation in vacuum first
sim = mp.Simulation(
    cell_size=cell_size,
    sources=sources,
    resolution=resolution,
    boundary_layers=pml_layers
)

incidentRegion = mp.FluxRegion(center=reflectionMonitorLocation, size=mp.Vector3(1.0, 1.0, 0), weight=1.0, direction=mp.X)
incidentFluxMonitor = sim.add_flux(frequency, frequencyWidth, numberFrequencies, incidentRegion)

sim.run(until_after_sources=400 * rescaling)

incidentFluxToSubtract = sim.get_flux_data(incidentFluxMonitor)

# Perform simulation with the structure
sim = mp.Simulation(
    cell_size=cell_size,
    sources=sources,
    resolution=resolution,
    geometry=geometry,
    boundary_layers=pml_layers)

transmissionRegion = mp.FluxRegion(center=transmissionMonitorLocation, size=mp.Vector3(1.0, 0, 0), weight=1.0, direction=mp.X)
transmissionFluxMonitor = sim.add_flux(frequency, frequencyWidth, numberFrequencies, transmissionRegion)
reflectionRegion = incidentRegion
reflectionFluxMonitor = sim.add_flux(frequency, frequencyWidth, numberFrequencies, reflectionRegion)
sim.load_minus_flux_data(reflectionFluxMonitor, incidentFluxToSubtract)

sim.run(until_after_sources=400 * rescaling)