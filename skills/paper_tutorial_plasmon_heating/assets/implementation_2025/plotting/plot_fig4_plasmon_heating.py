import numpy as np
import matplotlib.pyplot as plt
import columnplots as clp

''''
This plotting script contains three main parts:

Fig. a: Plasmonic geometry + molecules demonstration loaded from ./plasmonic_setup.png (no xy axis needed)
Fig. b: Frequency-domain plasmonic mode after the Gaussian pulse excitation
Fig. c: Ey field intensity plot loaded from ey_field_intensity_zslice.npy
Fig. d-f: HCN final energy distribution along 2D plane following previous code
'''

def get_energy_data(filename):
    # get energy data from numpy npz file
    data = np.load(filename)
    time_au = data['time_au']
    energy_au = data['energy_au']
    energy_au -= energy_au[0]  # shift to zero at t=0
    print(f"Loaded data from {filename}: time_au shape = {time_au.shape}, energy_au shape = {energy_au.shape}")
    return time_au, energy_au

def get_energy_distribution(nmol=64, path="meep_plasmon_HCN_excitation_bomd_strong"):
    e_distr_lst = []
    nmol_per_dim = int(round(nmol**(1/2)))
    for idx in range(nmol):
        filename = f"../{path}/nmol_{nmol}_with_dielectric/mol_{idx}_data.npz"
        # filename = f"../meep_plasmon_HCN_excitation_rteh_strong/nmol_{nmol}_with_dielectric/mol_{idx}_data.npz"
        time_au, energy_au = get_energy_data(filename)
        e_distr_lst.append(np.mean(energy_au))
    e_distr = np.array(e_distr_lst).reshape((nmol_per_dim, nmol_per_dim)).transpose()
    e_distr -= np.min(np.min(e_distr))
    return e_distr

def get_spectrum():
    filename_empty = "../meep_plasmon_empty/vac/flux0_a2.79.dat"
    f0 = np.genfromtxt(filename_empty, delimiter=",")
    filename_geometry = "../meep_plasmon_empty/no_mol_with_dielectric/flux_a2.79_r1.11.dat"
    f = np.genfromtxt(filename_geometry, delimiter=",")

    refl = -f[:,1]/f0[:,1]
    abs = 1 - refl
    return f0[:,0] * 1e4, abs

# preload important data
plasmon_imag = plt.imread("./plasmon_setup.png")

freq, sp = get_spectrum()

ey_field_intensity = np.load("ey_field_intensity_zslice.npy").transpose()

energy_distribution_bomd = get_energy_distribution(nmol=256, path="meep_plasmon_HCN_excitation_bomd_strong")
energy_distribution_rteh = get_energy_distribution(nmol=256, path="meep_plasmon_HCN_excitation_rteh_strong")
energy_distribution_tls = get_energy_distribution(nmol=256, path="meep_plasmon_HCN_excitation_tls_strong")


# create the figure with three subplots
fig, axes = clp.initialize(2, 3, width=4.3*3, height=4.3*2*0.84, return_fig_args=True,
                           fontsize=12, LaTeX=True)

# Fig. a: plasmonic geometry

axes[0, 0].imshow(plasmon_imag)
axes[0, 0].axis('off')

# Fig. b: Frequency-domain plasmonic mode
xs, ys = [freq], [sp]
clp.plotone(xs, ys, axes[0,1], colors=[clp.red], linestyles=['-'], xlabel='frequency (cm$^{-1}$)', ylabel='absorption', 
            xlim=[2000, 5000], ylim=[0, 1], labels=["flux spectrum"])
axes[0, 1].axvline(x=3465.5930, color=clp.dark_green, linestyle='--')
# add text to label the vertical line as the HCN stretch frequency
axes[0, 1].text(3600, 0.04, 'C-H stretch of HCN\n(B3LYP/cc-pVDZ)', color=clp.dark_green)

# Fig. c: Ey field intensity
im1 = axes[0, 2].imshow(ey_field_intensity / np.max(np.max(ey_field_intensity)), origin='lower', extent=[-1.4, 1.4, -1.4, 1.4], cmap='inferno')
cbar1 = fig.colorbar(im1, ax=axes[0, 2], fraction=0.046, pad=0.04)
cbar1.set_label(r'$E_y$ field intensity [arb. units]')
axes[0, 2].set_xlabel(r'$x$ ($\mu$m)')
axes[0, 2].set_ylabel(r'$y$ ($\mu$m)')

# Fig. d: HCN final energy distribution using tls
im2 = axes[1, 0].imshow(energy_distribution_tls, origin='lower', extent=[-1.4, 1.4, -1.4, 1.4], cmap='viridis')
cbar2 = fig.colorbar(im2, ax=axes[1, 0], fraction=0.046, pad=0.04)
cbar2.set_label(r'molecular energy gain [a.u.]')
axes[1, 0].set_xlabel(r'$x$ ($\mu$m)')
axes[1, 0].set_ylabel(r'$y$ ($\mu$m)')


# Fig. e: HCN final energy distribution using bomd
im2 = axes[1, 1].imshow(energy_distribution_bomd, origin='lower', extent=[-1.4, 1.4, -1.4, 1.4], cmap='viridis')
cbar2 = fig.colorbar(im2, ax=axes[1, 1], fraction=0.046, pad=0.04)
cbar2.set_label(r'molecular energy gain [a.u.]')
axes[1, 1].set_xlabel(r'$x$ ($\mu$m)')
# turn off y axis number labels for the last subplot
axes[1, 1].set_yticklabels([])
# axes[1, 0].set_ylabel(r'$y$ ($\mu$m)')

# Fig. f: HCN final energy distribution using rteh
im2 = axes[1, 2].imshow(energy_distribution_rteh, origin='lower', extent=[-1.4, 1.4, -1.4, 1.4], cmap='viridis')
cbar2 = fig.colorbar(im2, ax=axes[1, 2], fraction=0.046, pad=0.04)
cbar2.set_label(r'molecular energy gain [a.u.]')
axes[1, 2].set_xlabel(r'$x$ ($\mu$m)')
# turn off y axis number labels for the last subplot
axes[1, 2].set_yticklabels([])
# axes[1, 0].set_ylabel(r'$y$ ($\mu$m)')

# add labels (a), (b), (c), (d) to subplots
axes[0, 0].text(0.11, 0.95, '(a)', transform=axes[0, 0].transAxes, fontsize=14, fontweight='bold', va='top', ha='right', color='white')
axes[0, 1].text(0.11, 0.95, '(b)', transform=axes[0, 1].transAxes, fontsize=14, fontweight='bold', va='top', ha='right', color='k')
axes[0, 2].text(0.11, 0.95, '(c)', transform=axes[0, 2].transAxes, fontsize=14, fontweight='bold', va='top', ha='right', color='white')
axes[1, 0].text(0.11, 0.95, '(d)', transform=axes[1, 0].transAxes, fontsize=14, fontweight='bold', va='top', ha='right', color='white')
axes[1, 1].text(0.11, 0.95, '(e)', transform=axes[1, 1].transAxes, fontsize=14, fontweight='bold', va='top', ha='right', color='white')
axes[1, 2].text(0.11, 0.95, '(f)', transform=axes[1, 2].transAxes, fontsize=14, fontweight='bold', va='top', ha='right', color='white')

# add text about method in molecular energy distribution plots
axes[1, 0].text(0.08, 0.85, 'TLS', transform=axes[1, 0].transAxes, fontsize=12, fontweight='bold', va='top', ha='center', color='white')
axes[1, 1].text(0.12, 0.85, 'BOMD', transform=axes[1, 1].transAxes, fontsize=12, fontweight='bold', va='top', ha='center', color='white')
axes[1, 2].text(0.13, 0.85, 'RT-Eh', transform=axes[1, 2].transAxes, fontsize=12, fontweight='bold', va='top', ha='center', color='white')
clp.adjust(tight_layout=True, savefile="fig4_plasmon_heating.pdf")


