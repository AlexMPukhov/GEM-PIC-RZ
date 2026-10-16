"""Fig. LWFA benchmark: on-axis E_z and E_x of GEM-PIC-RZ and FBPIC (benchmarks/lwfa_fbpic).
(a) E_z at 50 and 100 um with plasma.smooth_length = 2 k0^-1; (b) laser field E_x at 100 um.
(The comparison behind the bubble, cold and smoothed, is in fig_closure.py.)
usage (from paper/): python3 scripts/fig_lwfa.py  ->  figs/fig_lwfa.pdf"""
import glob
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

B = "../benchmarks/lwfa_fbpic/"
E0 = 4.0126e12  # m c omega0 / e at 0.8 um [V/m]


def axis(run, D):
    for fn in sorted(glob.glob(B + run + "/axis_*.txt")):
        z = float(open(fn).readline().split("z =")[1].split("mm")[0])
        if abs(z * 1000 - D) < 1e-6:
            return np.loadtxt(fn)
    raise SystemExit(f"{run}: no axis file at {D} um")


plt.rcParams.update({"font.size": 9, "font.family": "serif", "mathtext.fontset": "cm"})
fig, axs = plt.subplots(2, 1, figsize=(3.4, 4.0), gridspec_kw={"height_ratios": [1.2, 0.8]})
for D, col, fcol in [(50, "#999999", "#e8a0a8"), (100, "#1f4e9c", "#d1495b")]:
    f = np.loadtxt(B + f"fbpic_D{D}.txt")
    g = axis("gem_out_a2", D)
    m = (f[:, 0] >= 0) & (f[:, 0] <= 38.2)
    axs[0].plot(g[:, 0], g[:, 2], "-", color=col, lw=1.2, label=f"GEM-PIC-RZ, {D} µm")
    axs[0].plot(f[m, 0], f[m, 1] / E0, "--", color=fcol, lw=1.0, label=f"FBPIC, {D} µm")
f = np.loadtxt(B + "fbpic_D100.txt")
g = axis("gem_out_a2", 100)
mm = (f[:, 0] >= 0) & (f[:, 0] <= 38.2)
axs[1].plot(g[:, 0], g[:, 1], "-", color="#1f4e9c", lw=0.7, label="GEM-PIC-RZ")
axs[1].plot(f[mm, 0], f[mm, 2] / E0, "--", color="#d1495b", lw=0.7, label="FBPIC")
axs[0].set_ylabel(r"$E_z(r=0)$ $[mc\omega_0/e]$")
axs[1].set_ylabel(r"$E_x(r=0)$ $[mc\omega_0/e]$")
for a in axs:
    a.set_xlim(0, 38.2)
axs[1].set_xlabel(r"$\zeta$ [µm]")
axs[0].legend(fontsize=6, loc="lower left", frameon=False)
axs[1].legend(fontsize=6, loc="upper left", frameon=False)
for a, t in zip(axs, ["(a)", "(b)"]):
    a.text(0.97, 0.88, t, transform=a.transAxes, ha="right")
fig.tight_layout(h_pad=0.4)
fig.savefig("figs/fig_lwfa.pdf")
