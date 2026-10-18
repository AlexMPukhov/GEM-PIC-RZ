"""Fig. 1: coherent emission of a short bunch in a planar undulator, spectrum per azimuthal mode, GEM-PIC-RZ
(benchmarks/harmonics/em_long) against the radiation integral (radint_em_modes.txt)."""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

H = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "benchmarks", "harmonics")
R = np.loadtxt(os.path.join(H, "radint_em_modes.txt"))
S0 = np.loadtxt(os.path.join(H, "spectra", "em_long_step0.txt"))
S1 = np.loadtxt(os.path.join(H, "spectra", "em_long_step420.txt"))
k = S1[:, 0]
code = S1[:, 3:7].copy()
code[:, 0] -= S0[:, 3]
plt.rcParams.update({"font.size": 8, "axes.linewidth": 0.6})
fig, ax = plt.subplots(figsize=(3.4, 2.5))
cols = ["#1b6ca8", "#d1495b", "#2a9d8f", "#e9a23b"]
for m in range(4):
    ax.semilogy(R[:, 0], R[:, 1 + m], color=cols[m], lw=2.2, alpha=0.3)
    ax.semilogy(k, np.maximum(code[:, m], 1e-12), color=cols[m], lw=0.8, label=f"$m = {m}$")
ax.set_xlim(0.3, 3.6)
ax.set_ylim(3e-10, 5e-6)
ax.set_xlabel(r"$k/k_r$")
ax.set_ylabel(r"$\mathrm{d}W/\mathrm{d}k$ (code units)")
ax.legend(frameon=False, fontsize=7, loc="lower left", ncol=2)
ax.text(0.98, 0.95, "thin: GEM-PIC-RZ\nwide: radiation integral", transform=ax.transAxes, ha="right", va="top",
        fontsize=6.5)
fig.tight_layout(pad=0.3)
fig.savefig(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "figs", "fig_emission.pdf"))
