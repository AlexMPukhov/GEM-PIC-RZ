"""Fig. LWFA benchmark: on-axis E_z and E_x of GEM-PIC-RZ and FBPIC (benchmarks/lwfa_fbpic).
usage (from paper/): python3 scripts/fig_lwfa.py  ->  figs/fig_lwfa.pdf"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

B = "../benchmarks/lwfa_fbpic/"
E0 = 4.0126e12  # m c omega0 / e at 0.8 um [V/m]
plt.rcParams.update({"font.size": 9, "font.family": "serif", "mathtext.fontset": "cm"})
fig, axs = plt.subplots(2, 1, figsize=(3.4, 4.2), sharex=True)
for D, step, col in [(50, 20, "#999999"), (100, 40, "#1f4e9c")]:
    f = np.loadtxt(B + f"fbpic_D{D}.txt")
    g = np.loadtxt(B + f"gem_out/axis_{step:05d}.txt")
    m = (f[:, 0] >= 0) & (f[:, 0] <= 38.2)
    axs[0].plot(g[:, 0], g[:, 2], "-", color=col, lw=1.2, label=f"GEM-PIC-RZ, {D} µm")
    axs[0].plot(f[m, 0], f[m, 1] / E0, "--", color="#d1495b" if D == 100 else "#e8a0a8", lw=1.0,
                label=f"FBPIC, {D} µm")
    if D == 100:
        axs[1].plot(g[:, 0], g[:, 1], "-", color=col, lw=0.7, label="GEM-PIC-RZ")
        axs[1].plot(f[m, 0], f[m, 2] / E0, "--", color="#d1495b", lw=0.7, label="FBPIC")
axs[0].set_ylabel(r"$E_z(r=0)$ $[mc\omega_0/e]$")
axs[1].set_ylabel(r"$E_x(r=0)$ $[mc\omega_0/e]$")
axs[1].set_xlabel(r"$\zeta$ [µm]")
axs[0].legend(fontsize=6.5, loc="lower left", frameon=False)
axs[1].legend(fontsize=6.5, loc="upper left", frameon=False)
axs[1].set_xlim(0, 38.2)
for a, t in zip(axs, ["(a)", "(b)"]):
    a.text(0.97, 0.9, t, transform=a.transAxes, ha="right")
fig.tight_layout()
fig.savefig("figs/fig_lwfa.pdf")
