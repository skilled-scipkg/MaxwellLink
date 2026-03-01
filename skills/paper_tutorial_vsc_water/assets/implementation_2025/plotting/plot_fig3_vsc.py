import numpy as np
import matplotlib.pyplot as plt
import maxwelllink as mxl
from maxwelllink.tools import ir_spectrum
import columnplots as clp
import matplotlib.image as mpimg
from matplotlib.offsetbox import OffsetImage, AnnotationBbox


def gather_data(mode="singlemode", size=2):
    filename=f"../vsc_water_singlemode/size_{size}_cube/vsc_water_maxwellmd_data.npz"
    if mode=="maxwellmd":
        filename=f"../vsc_water_maxwellmd/size_{size}_cube/vsc_water_maxwellmd_data.npz"
    data = np.load(filename)
    time_au = data["time_au"]
    mux_au = data["mux_au"]
    muy_au = data["muy_au"]
    muz_au = data["muz_au"]
    energy_molecule_au = data["energy_molecule_au"]
    # do a fourier transfer
    time_fs = time_au * 0.02418884  # convert atomic unit time to fs
    dt = time_fs[1] - time_fs[0]
    freq, sp_x = ir_spectrum(mux_au, dt, field_description="square", smooth_window_len=25)
    freq, sp_y = ir_spectrum(muy_au, dt, field_description="square", smooth_window_len=25)
    freq, sp_z = ir_spectrum(muz_au, dt, field_description="square", smooth_window_len=25)
    sp_uncoupled = (sp_x + sp_y) / 2.0
    sp_coupled = sp_z
    # renormalized
    dfreq = freq[1] - freq[0]
    idx_max = int(1200 // dfreq)
    sp_uncoupled /= np.max(sp_uncoupled[0:idx_max])
    sp_coupled /= np.max(sp_coupled[0:idx_max])

    # also capture the time count
    filename_time = filename.replace("vsc_water_maxwellmd_data.npz", "time_count.txt")
    data_time = np.loadtxt(filename_time)
    time_per_step = np.mean(data_time)
    print(f"Mode: {mode}, Size: {size}^3, Time per step: {time_per_step:.4f} s")
    
    # standalone lammps simulation time
    filename_lammps = filename.replace("vsc_water_maxwellmd_data.npz", "lammps_count.txt")
    data_lammps = np.loadtxt(filename_lammps)
    lammps_time_per_step = 1.0 / data_lammps
    print(f"Mode: {mode}, Size: {size}^3, LAMMPS time per step: {lammps_time_per_step:.4f} s")
    return freq, sp_uncoupled, sp_coupled, time_per_step, lammps_time_per_step

size_lst = [1, 2, 4, 8]
xs_singlemode, ys_uncoupled_singlemode, ys_coupled_singlemode, time_singlemode, lammps_singlemode = [], [], [], [], []
xs_maxwellmd, ys_uncoupled_maxwellmd, ys_coupled_maxwellmd, time_maxwellmd, lammps_maxwellmd = [], [], [], [], []
increment = 1

for idx, size in enumerate(size_lst):
    freq, sp_uncoupled, sp_coupled, time_per_step, lammps_time_per_step = gather_data(mode="singlemode", size=size)
    xs_singlemode.append(freq)
    ys_uncoupled_singlemode.append(sp_uncoupled + increment*idx)
    ys_coupled_singlemode.append(sp_coupled + increment*(idx+1))
    time_singlemode.append(time_per_step)
    lammps_singlemode.append(lammps_time_per_step)

    freq, sp_uncoupled, sp_coupled, time_per_step, lammps_time_per_step = gather_data(mode="maxwellmd", size=size)
    xs_maxwellmd.append(freq)
    ys_uncoupled_maxwellmd.append(sp_uncoupled * 0.5 + increment*idx)
    ys_coupled_maxwellmd.append(sp_coupled *0.5 + increment*(idx+1))
    time_maxwellmd.append(time_per_step)
    lammps_maxwellmd.append(lammps_time_per_step)

ncpu_lst = np.array([size**3 for size in size_lst])
time_singlemode = np.array(time_singlemode)
time_maxwellmd = np.array(time_maxwellmd)
lammps_singlemode = np.array(lammps_singlemode)
lammps_maxwellmd = np.array(lammps_maxwellmd)

xs_count = [ncpu_lst]*2
ys_count = [(time_singlemode + time_maxwellmd)/2.0, (lammps_singlemode + lammps_maxwellmd)/2.0]

print("ncpu_lst", ncpu_lst)
print("Single-mode times per step (s):", time_singlemode)
print("MaxwellMD times per step (s):", time_maxwellmd)

# plotting
labels = [r"cavity off"] + [r"$N_{\text{H}_2\text{O}}=216$"] + [r"$N_{\text{H}_2\text{O}} = 216 \times %d^3$" % (size) for size in size_lst[1:]]
fig, axes = clp.initialize(3, 1, width=4.3, height=4.3*0.618*3, return_fig_args=True,
                           fontsize=12, labelthem=True, labelthemPosition=[0.1, 0.95],
                           labelsize=14,
                           LaTeX=True)

# fig a: single-mode VSC spectrum
clp.plotone([xs_singlemode[0]] + xs_singlemode, [ys_uncoupled_singlemode[0]] + ys_coupled_singlemode, axes[0], xlim=[0,4500], labels=labels,
            xlabel='frequency [cm$^{-1}$]', ylabel='IR spectrum [arb. units]', showlegend=True,
            colormap=plt.cm.hot, legendFontSize=8, alphaspacing=0.1)

axes[0].axvline(x=3550, color=clp.red, linestyle='--', linewidth=1.5, alpha=0.5)

# fig b: MaxwellMD VSC spectrum
clp.plotone([xs_maxwellmd[0]] + xs_maxwellmd, [ys_uncoupled_maxwellmd[0]] + ys_coupled_maxwellmd, axes[1], xlim=[0,4500], 
            xlabel='frequency [cm$^{-1}$]', ylabel='IR spectrum [arb. units]', showlegend=False,
            colormap=plt.cm.hot, alphaspacing=0.1)

# add inset line for bragg resonator transmission
data_bragg = np.loadtxt("../vsc_water_maxwellmd/bragg_resonator_spectrum/bragg_resonator_transmission.txt")
freq_bragg, sp_bragg = data_bragg[:,0], data_bragg[:,1]
axes_inset = axes[1].inset_axes([0.25, 0.61, 0.42, 0.37])
axes_inset.plot(freq_bragg, sp_bragg, color=clp.navy_blue)
axes_inset.set_xlim([0, 6000])
axes_inset.set_ylim([0, 1.3])
# add vertical lines to indicate H2O vibrational frequencies
for freq_line in [700, 1650, 3550]:
    axes_inset.axvline(x=freq_line, color=clp.darkgray, linestyle='--', linewidth=1.0, alpha=0.5)
# add text to indicate these vertical lines for H2O vibrations
axes_inset.text(850, 0.2, 'H$_2$O', rotation=90, fontsize=8, color=clp.black, alpha=0.7)
axes_inset.set_xlabel('frequency [cm$^{-1}$]', fontsize=8)
axes_inset.set_ylabel('transmission', fontsize=8)
axes_inset.tick_params(axis='both', which='major', labelsize=8)

# add inset for bragg geometry
img = mpimg.imread("bragg_geometry.png")
imagebox = OffsetImage(img, zoom=0.16)
ab = AnnotationBbox(imagebox, (0.46, 0.935), xycoords='axes fraction', frameon=False, zorder=10)
axes[1].add_artist(ab)


# fig c: time per step vs # CPU cores
clp.plotone(xs_count[0:2], ys_count[0:2], axes[2], colors=[clp.red, clp.darkgray], linestyles=[':', ':'],
            markers=['o', '^'], xlabel=r'\# CPU cores for driver ($=N_{\text{H}_2\text{O}}/216$)', 
            ylabel='stepping time [s]', xlog=True, ylog=True, ylim=[0.001, 0.1],
            labels=['MaxwellLink+LAMMPS', 'standalone LAMMPS'], )

clp.adjust(tight_layout=True, savefile="fig3_vsc.pdf")
