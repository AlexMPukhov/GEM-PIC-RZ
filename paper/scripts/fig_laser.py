"""Fig. laser evolution benchmark (benchmarks/laser_evolution): laser energy, peak field on the axis and spectral
centroid vs propagation distance; GEM-PIC-RZ, FBPIC and Wake-T (envelope); on-axis E_z at the end.
usage (from paper/): python3 scripts/fig_laser.py  ->  figs/fig_laser.pdf"""
import glob
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

B = "../benchmarks/laser_evolution/"
fb = np.loadtxt(B + "fbpic_out/laser_fbpic.txt")
fb[:, 1] /= fb[0, 1]
wt = np.loadtxt(B + "waket_out/laser_waket.txt")
runs = [("laser_gem_out_a1m3.txt", "#1f4e9c", "-", r"GEM-PIC-RZ, $\Delta = 2\,\mu$m"),
        ("laser_gem_out_a1m3_dz1.txt", "#6fa8dc", ":", r"GEM-PIC-RZ, $\Delta = 1\,\mu$m")]
plt.rcParams.update({"font.size": 9, "font.family": "serif", "mathtext.fontset": "cm"})
fig, axs = plt.subplots(4, 1, figsize=(3.4, 7.0), gridspec_kw={"height_ratios": [1, 1, 1, 1.1]})
for j, (col, lab) in enumerate([(1, r"$W/W_0$"), (2, r"max $|E_x(r=0)|$ $[mc\omega_0/e]$"),
                                 (3, r"$\langle |k| \rangle / k_0$")]):
    a = axs[j]
    a.plot(fb[:, 0], fb[:, col], "-", color="k", lw=2.2, alpha=0.35, label="FBPIC")
    for fn, c, ls, lb in runs:
        g = np.loadtxt(B + fn)
        g[:, 1] /= g[0, 1]
        a.plot(g[:, 0], g[:, col], ls, color=c, lw=1.2, label=lb)
    a.plot(wt[:, 0], wt[:, col], "--", color="#d1495b", lw=1.1, label="Wake-T (envelope)")
    a.set_ylabel(lab)
    a.set_xlim(0, 300)
    a.set_xticklabels([])
axs[0].legend(fontsize=6, frameon=False, loc="lower left")


def last(pattern):
    return sorted(glob.glob(B + pattern))[-1]


z_end = 300
f = np.loadtxt(last("fbpic_out/axis_*.txt"))
g = np.loadtxt(B + "gem_out_a1m3/axis_00150.txt")
w = np.loadtxt(B + f"waket_out/axis_{z_end:07.2f}.txt")
m = (f[:, 0] >= 0) & (f[:, 0] <= 24.95)
axs[3].plot(f[m, 0], f[m, 1], "-", color="k", lw=2.2, alpha=0.35, label="FBPIC")
axs[3].plot(g[:, 0], g[:, 2], "-", color="#1f4e9c", lw=1.0, label="GEM-PIC-RZ")
if w.shape[1] > 2:
    axs[3].plot(w[:, 0], w[:, 2], "--", color="#d1495b", lw=1.0, label="Wake-T")
axs[3].set_xlabel(r"$\zeta$ [$\mu$m]")
axs[3].set_ylabel(r"$E_z(r=0)$ at 300 $\mu$m")
axs[3].set_xlim(0, 24.95)
axs[2].set_xlabel(r"propagation distance [$\mu$m]")
axs[2].set_xticks([0, 100, 200, 300])
axs[2].set_xticklabels(["0", "100", "200", "300"])
for a, t in zip(axs, ["(a)", "(b)", "(c)", "(d)"]):
    a.text(0.97, 0.85, t, transform=a.transAxes, ha="right")
fig.tight_layout(h_pad=0.3)
fig.savefig("figs/fig_laser.pdf")
