"""Fig. 3: 1D steady-state planar FEL, period-averaged vs unaveraged electron motion (benchmarks/high_rho/scan.txt):
gain length and saturation power against rho."""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

S = np.loadtxt(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "benchmarks", "high_rho",
                            "scan.txt"))
plt.rcParams.update({"font.size": 8, "axes.linewidth": 0.6})
fig, ax = plt.subplots(2, 1, figsize=(3.4, 3.6), sharex=True)
for K, c in ((1.0, "#1b6ca8"), (2.0, "#d1495b"), (4.0, "#2a9d8f")):
    r = S[S[:, 0] == K]
    ax[0].semilogx(r[:, 1], r[:, 3], "o--", ms=2.5, lw=0.8, color=c, mfc="none")
    ax[0].semilogx(r[:, 1], r[:, 6], "o-", ms=2.5, lw=0.8, color=c, label=f"$K = {K:g}$")
    ax[1].semilogx(r[:, 1], r[:, 4], "o--", ms=2.5, lw=0.8, color=c, mfc="none")
    ax[1].semilogx(r[:, 1], r[:, 7], "o-", ms=2.5, lw=0.8, color=c)
ax[0].set_ylabel(r"$L_g / L_{g,1D}$")
ax[1].set_ylabel(r"$P_\mathrm{sat} / \rho P_b$")
ax[1].set_xlabel(r"$\rho$")
ax[0].legend(frameon=False, fontsize=7)
ax[0].text(0.03, 0.08, "(a)  solid: unaveraged, dashed: averaged", transform=ax[0].transAxes, fontsize=6.5)
ax[1].text(0.03, 0.08, "(b)", transform=ax[1].transAxes)
fig.tight_layout(pad=0.3)
fig.savefig(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "figs", "fig_averaging.pdf"))
