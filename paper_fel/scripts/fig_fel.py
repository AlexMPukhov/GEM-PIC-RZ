"""Fig. 2: seeded planar FEL (benchmarks/harmonics): power at harmonics 1-3 and bunching, GEM-PIC-RZ (fel_tables)
against GENESIS 1.3 v4 (genesis/harm.out.h5), steady-state slices."""
import os
import h5py
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

H = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "benchmarks", "harmonics")
zlo, zhi = 2.7, 3.45
f = os.path.join(H, "fel_tables", "fel_harmonics.txt")
hdr = [l for l in open(f) if l.startswith("# step")][0]
zj = np.array([float(v) for v in hdr.split(":")[1].split()])
sel = (zj >= zlo - 1e-9) & (zj + (zj[1] - zj[0]) <= zhi + 1e-9)
d = np.loadtxt(f)
G = {(h, m): (d[(d[:, 2] == h) & (d[:, 3] == m), 1] * 1e-3, d[(d[:, 2] == h) & (d[:, 3] == m)][:, 4:][:, sel].mean(1))
     for h in (1, 2, 3) for m in (-1, 0, 1, 2)}
g = h5py.File(os.path.join(H, "genesis", "harm.out.h5"), "r")
zs = g["Lattice/zplot"][()]
Ps = [g[k + "/power"][()].ravel() for k in ("Field", "Field2", "Field3")]
bs = [g["Beam/bunching" + s][()].ravel() for s in ("", "2", "3")]

plt.rcParams.update({"font.size": 8, "axes.linewidth": 0.6})
cols = {1: "#1b6ca8", 2: "#d1495b", 3: "#2a9d8f"}
fig, ax = plt.subplots(2, 1, figsize=(3.4, 4.4), sharex=True)
ax[0].semilogy(zs, Ps[0], "--", color=cols[1], lw=1)
ax[0].semilogy(zs, Ps[2], "--", color=cols[3], lw=1)
ax[0].semilogy(G[1, 1][0], G[1, 1][1], color=cols[1], lw=1.2, label="$h = 1$")
ax[0].semilogy(G[2, 2][0], 2 * G[2, 2][1], color=cols[2], lw=1.2, label=r"$h = 2$ ($2P_{m=2}$)")
ax[0].semilogy(G[2, 0][0], G[2, 0][1], ":", color=cols[2], lw=0.9, label=r"$h = 2$, $m = 0$")
ax[0].semilogy(G[3, 1][0], G[3, 1][1], color=cols[3], lw=1.2, label="$h = 3$")
ax[0].set_ylim(1e2, 3e10)
ax[0].set_ylabel("power (W)")
ax[0].legend(frameon=False, fontsize=6.5, loc="lower right")
ax[0].text(0.03, 0.95, "(a)  solid: GEM-PIC-RZ\n       dashed: GENESIS (no $h = 2$)", transform=ax[0].transAxes,
           va="top", fontsize=6.5)
for h in (1, 2, 3):
    ax[1].plot(zs, bs[h - 1], "--", color=cols[h], lw=1)
    ax[1].plot(G[h, -1][0], G[h, -1][1], color=cols[h], lw=1.2, label=f"$b_{h}$")
ax[1].set_xlabel("$z$ (m)")
ax[1].set_ylabel("bunching factor")
ax[1].legend(frameon=False, fontsize=7)
ax[1].text(0.03, 0.95, "(b)", transform=ax[1].transAxes, va="top")
fig.tight_layout(pad=0.3)
fig.savefig(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "figs", "fig_fel.pdf"))
