import meep as mp
import math
import cmath
import argparse
import maxwelllink as mxl
from maxwelllink import sockets as mxs
import numpy as np

# set up SocketHub and save host/port info to a file
host, port = mxs.get_available_host_port(localhost=False, save_to_file="tcp_host_port_info.txt")
hub = mxl.SocketHub(host=host, port=port, timeout=200.0, latency=1e-4)

def main(args):

    # default unit length is 1 um
    um_scale = 1.0
    # with speed of light c=1, default unit time is 1 um/c = 3.33564 fs, needed for MaxwellLink simulations
    time_units_fs = 3.33564

    # conversion factor for eV to 1/um [=1/hc]
    eV_um_scale = um_scale/1.23984193

    # Pt from A.D. Rakic et al., Applied Optics, Vol. 37, No. 22, pp. 5271-83, 1998
    Pt_plasma_frq = 9.59*eV_um_scale
    Pt_f0 = 0.333
    Pt_frq0 = 1e-10
    Pt_gam0 = 0.080*eV_um_scale
    Pt_sig0 = Pt_f0*Pt_plasma_frq**2/Pt_frq0**2
    Pt_f1 = 0.191
    Pt_frq1 = 0.780*eV_um_scale      # 1.590 um
    Pt_gam1 = 0.517*eV_um_scale
    Pt_sig1 = Pt_f1*Pt_plasma_frq**2/Pt_frq1**2
    Pt_f2 = 0.659
    Pt_frq2 = 1.314*eV_um_scale      # 0.944 um
    Pt_gam2 = 1.838*eV_um_scale
    Pt_sig2 = Pt_f2*Pt_plasma_frq**2/Pt_frq2**2

    Pt_susc = [ mp.DrudeSusceptibility(frequency=Pt_frq0, gamma=Pt_gam0, sigma=Pt_sig0),
                mp.LorentzianSusceptibility(frequency=Pt_frq1, gamma=Pt_gam1, sigma=Pt_sig1),
                mp.LorentzianSusceptibility(frequency=Pt_frq2, gamma=Pt_gam2, sigma=Pt_sig2) ]

    Pt = mp.Medium(epsilon=1.0, E_susceptibilities=Pt_susc)    
    
    resolution = 40     # pixels/um

    a = args.aa         # lattice periodicity
    r = args.rr         # metal rod radius
    nmol = args.nmol    # number of molecules
    nmol_per_dim = int(round(nmol**(1/2)))
    print(f"#### a = {a}, r = {r}, nmol = {nmol}, nmol_per_dim = {nmol_per_dim}")
    h = 0.2             # metal rod height
    tmet = 0.3          # metal layer thickness
    tsub = 2.0          # substrate thickness
    tabs = 5.0          # PML thickness
    tair = 4.0          # air thickness

    sz = tabs+tair+h+tmet+tsub+tabs
    cell_size = mp.Vector3(a,a,sz)

    pml_layers = [mp.PML(thickness=tabs,direction=mp.Z,side=mp.High),
                  mp.Absorber(thickness=tabs,direction=mp.Z,side=mp.Low)]

    lmin = 2.0          # source min wavelength
    lmax = 5.0          # source max wavelength
    fmin = 1/lmax       # source min frequency
    fmax = 1/lmin       # source max frequency
    fcen = 0.5*(fmin+fmax)
    df = fmax-fmin

    nSi = 3.5
    Si = mp.Medium(index=nSi)

    geometry = []
        
    if not args.empty and args.dielectric:
        geometry = [ mp.Cylinder(material=Pt, radius=r, height=h, center=mp.Vector3(0,0,0.5*sz-tabs-tair-0.5*h)),
                     mp.Block(material=Pt, size=mp.Vector3(mp.inf,mp.inf,tmet),
                              center=mp.Vector3(0,0,0.5*sz-tabs-tair-h-0.5*tmet)),
                     mp.Block(material=Si, size=mp.Vector3(mp.inf,mp.inf,tsub+tabs),
                              center=mp.Vector3(0,0,0.5*sz-tabs-tair-h-tmet-0.5*(tsub+tabs))) ]

    # CCW rotation angle (degrees) about Y-axis of PW current source; 0 degrees along -z axis
    theta = math.radians(args.theta)

    # k with correct length (plane of incidence: XZ) 
    k = mp.Vector3(math.sin(theta),0,math.cos(theta)).scale(fcen)

    def pw_amp(k, x0):
        def _pw_amp(x):
            return 1e2 * cmath.exp(1j * 2 * math.pi * k.dot(x + x0))
        return _pw_amp

    src_pos = 0.5*sz-tabs-0.2*tair
    sources = [ mp.Source(mp.GaussianSource(fcen, fwidth=df), component=mp.Ey, center=mp.Vector3(0,0,src_pos),
                          size=mp.Vector3(a,a,0),
                          amp_func=pw_amp(k, mp.Vector3(0,0,src_pos))) ]
    
    # if this is not an empty simulation and molecule flag is set, add MaxwellLink molecule
    molecules = []
    len_molecule = 0.25
    z_top_of_cylinder = 0.5*sz - tabs - tair
    if not args.empty and args.mol:
        print("### adding MaxwellLink molecules")
        # create abstract molecules
        z_end_of_cylinder = z_top_of_cylinder - h
        x_left = -0.5*a + len_molecule
        x_right = 0.5*a - len_molecule
        y_left = -0.5*a + len_molecule
        y_right = 0.5*a - len_molecule
        print(f"### z_top_of_cylinder = {z_top_of_cylinder}, len_molecule = {len_molecule}")

        for idx_x in range(nmol_per_dim):
            for idx_y in range(nmol_per_dim):
                if len(molecules) >= nmol:
                    break
                idx = idx_x*nmol_per_dim + idx_y
                x_pos = x_left + idx_x * (x_right - x_left) / (nmol_per_dim - 1) if nmol > 1 else 0.0
                y_pos = y_left + idx_y * (y_right - y_left) / (nmol_per_dim - 1) if nmol > 1 else 0.0
                print(f"### adding molecule idx = {idx} at x = {x_pos}, y = {y_pos}, z = {z_top_of_cylinder + 0.5*len_molecule}")
                molecule = mxl.Molecule(
                            hub=hub,
                            center=mp.Vector3(x_pos, y_pos, z_top_of_cylinder + 0.5*len_molecule),
                            size=mp.Vector3(len_molecule, len_molecule, len_molecule),
                            sigma=len_molecule*0.1,
                            dimensions=3,
                            rescaling_factor=1.0)
                molecules.append(molecule)

        sim = mxl.MeepSimulation(
                        hub=hub,
                        molecules=molecules,
                        time_units_fs=time_units_fs,
                        cell_size=cell_size,
                        geometry=geometry,
                        sources=sources,
                        boundary_layers=pml_layers,
                        k_point = k,
                        resolution=resolution)
    else:
        sim = mp.Simulation(cell_size=cell_size,
                            geometry=geometry,
                            sources=sources,
                            boundary_layers=pml_layers,
                            k_point = k,
                            resolution=resolution)

    nfreq = 200
    refl = sim.add_flux(fcen, df, nfreq, mp.FluxRegion(center=mp.Vector3(0, 0, z_top_of_cylinder + 0.5*len_molecule),size=mp.Vector3(a,a,0)))

    if not args.empty:
        sim.load_minus_flux('refl-flux', refl)
    
    sim.run(until=300)

    if args.empty:
        sim.save_flux('refl-flux', refl)
    
    sim.display_fluxes(refl)

    # after simulation, also save molecular information to disk for post-processing
    if not args.empty and args.mol and mp.am_master():
        for idx, molecule in enumerate(molecules):
            time_au = np.array([ad["time_au"] for ad in molecule.additional_data_history])
            energy_au = np.array([ad["energy_au"] for ad in molecule.additional_data_history])
            mux_au = np.array([ad["mux_au"] for ad in molecule.additional_data_history])
            muy_au = np.array([ad["muy_au"] for ad in molecule.additional_data_history])
            muz_au = np.array([ad["muz_au"] for ad in molecule.additional_data_history])
            np.savez(
                f"mol_{idx}_data.npz",
                time_au=time_au,
                energy_au=energy_au,
                mux_au=mux_au,
                muy_au=muy_au,
                muz_au=muz_au,
            )

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-empty', action='store_true', default=False, help="empty? (default: False)")
    parser.add_argument('-mol', action='store_true', default=False, help="add MaxwellLink molecule? (default: False)")
    parser.add_argument('-dielectric', action='store_true', default=False, help="add classical dielectric medium? (default: False)")
    parser.add_argument('-aa', type=float, default=4.5, help='lattice periodicity (default: 4.5 um)')
    parser.add_argument('-rr', type=float, default=1.5, help='Pt rod radius (default: 1.5 um)')
    parser.add_argument('-theta', type=float, default=0, help='angle of planewave current source (default: 0 degrees)')
    parser.add_argument('-nmol', type=int, default=1, help='number of molecules (default: 1)')
    args = parser.parse_args()
    main(args)
