"""Fig. ionization injection benchmark (benchmarks/ionization_fbpic): charge and mean energy of the injected electrons
above 5 MeV vs the lab position of the bunch, and their spectrum at z = 183 um; GEM-PIC-RZ and FBPIC.
usage (from paper/): python3 scripts/fig_injection.py  ->  figs/fig_injection.pdf"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

B = "../benchmarks/ionization_fbpic/"
fb = np.loadtxt(B + "injection_fbpic.txt")
runs = [("injection_gem_out.txt", "#1f4e9c", "-", r"GEM-PIC-RZ, $\Delta = 2.5\,\mu$m"),
        ("injection_gem_out_dz1.25.txt", "#6fa8dc", ":", r"GEM-PIC-RZ, $\Delta = 1.25\,\mu$m")]


def ek(a):
    return (np.sqrt(1 + a[:, 3] ** 2 + a[:, 4] ** 2 + a[:, 5] ** 2) - 1) * 0.511


plt.rcParams.update({"font.size": 9, "font.family": "serif", "mathtext.fontset": "cm"})
fig, axs = plt.subplots(3, 1, figsize=(3.4, 5.6))
m = ~np.isnan(fb[:, 1])
for j, col in enumerate([5, 6]):
    axs[j].plot(fb[m, 1], fb[m, col], "o-", color="k", lw=2.2, alpha=0.35, ms=3, label="FBPIC")
    for fn, c, ls, lb in runs:
        g = np.loadtxt(B + fn)
        k = g[:, 5] > 0 if col == 6 else slice(None)  # mean energy only where there is charge
        axs[j].plot(g[k, 1], g[k, col], ls, color=c, lw=1.2, label=lb)
    axs[j].set_xlim(0, 200)
axs[0].set_ylabel(r"$Q\,(E > 5$ MeV) [pC]")
axs[1].set_ylabel(r"$\langle E \rangle\,(E > 5$ MeV) [MeV]")
axs[1].set_xlabel(r"lab position of the bunch [$\mu$m]")
axs[0].set_xticklabels([])
axs[0].legend(fontsize=6, frameon=False, loc="upper left")
f = np.loadtxt(B + "fbpic_out/elec_N_0200.00.txt.gz")
g = np.loadtxt(B + "gem_out_s182.5/beam_00073.txt.gz")
g = g[g[:, 7] == -30]
bins = np.linspace(2, 28, 27)
for E, q, c, ls, lb in [(ek(f), np.abs(f[:, 6]), "k", "-", r"FBPIC, $ct = 200\,\mu$m"),
                        (ek(g), np.abs(g[:, 6]), "#1f4e9c", "-", r"GEM-PIC-RZ, $s = 182.5\,\mu$m")]:
    h, _ = np.histogram(E, bins, weights=q)
    axs[2].step(bins[:-1], h / np.diff(bins), where="post", color=c, lw=2.0 if c == "k" else 1.2,
                alpha=0.35 if c == "k" else 1.0, label=lb)
axs[2].set_xlabel("kinetic energy [MeV]")
axs[2].set_ylabel("dQ/dE [pC/MeV]")
axs[2].legend(fontsize=6, frameon=False)
for a, t in zip(axs, ["(a)", "(b)", "(c)"]):
    a.text(0.97, 0.08, t, transform=a.transAxes, ha="right")
fig.tight_layout(h_pad=0.4)
fig.savefig("figs/fig_injection.pdf")
