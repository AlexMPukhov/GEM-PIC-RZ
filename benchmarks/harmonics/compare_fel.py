"""Harmonics of the seeded planar FEL (fel.deck): GEM-PIC-RZ against GENESIS 1.3 v4 (genesis/fel.in).
GEM-PIC-RZ: fel_harmonics.txt, averaged over the slices in the steady-state part of the bunch (radiation that has
met only the flat top). Harmonic h, mode m: odd harmonics on the axis in m = 1; the 2nd harmonic is the TEM10-like
field E_x ~ x, which has equal energy in m = 0 and m = 2; m = 0 also holds the co-moving field of the 2nd-harmonic
bunching, so P2 = 2 P(m = 2) is used (and P(m = 0) shown for reference).
usage: python3 compare_fel.py [out_dir] [png]"""
import sys
import h5py
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

out = sys.argv[1] if len(sys.argv) > 1 else "fel_out"
png = sys.argv[2] if len(sys.argv) > 2 else "compare_fel.png"
zlo, zhi = 2.7, 3.45          # steady-state zeta range [um] (flat top 2.62 ... 4.52 um, slippage 1.06 um)

hdr = [l for l in open(f"{out}/fel_harmonics.txt") if l.startswith("# step")][0]
zeta_j = np.array([float(v) for v in hdr.split(":")[1].split()])
dsl = zeta_j[1] - zeta_j[0]
sel = (zeta_j >= zlo - 1e-9) & (zeta_j + dsl <= zhi + 1e-9)
d = np.loadtxt(f"{out}/fel_harmonics.txt")
G = {}
for h in (1, 2, 3):
    for m in (-1, 0, 1, 2):
        r = d[(d[:, 2] == h) & (d[:, 3] == m)]
        if len(r):
            G[h, m] = (r[:, 1] * 1e-3, r[:, 4:][:, sel].mean(axis=1))
z_g = G[1, 1][0]
P1g, P3g = G[1, 1][1], G[3, 1][1]
P2g, P2g0 = 2 * G[2, 2][1], G[2, 0][1]

f = h5py.File("genesis/harm.out.h5", "r")
z_s = f["Lattice/zplot"][()]
P1s, P2s, P3s = (f[k + "/power"][()].ravel() for k in ("Field", "Field2", "Field3"))
b_s = [f["Beam/bunching" + s][()].ravel() for s in ("", "2", "3")]

i_s, i_g = np.argmax(P1s), np.argmax(P1g)
print(f"steady-state slices {zeta_j[sel][0]:.4g} ... {zeta_j[sel][-1] + dsl:.4g} um")
print(f"{'':28s} {'GENESIS':>12s} {'GEM-PIC-RZ':>12s}")
print(f"{'saturation z [m]':28s} {z_s[i_s]:12.3f} {z_g[i_g]:12.3f}")
print(f"{'P1 at saturation [W]':28s} {P1s[i_s]:12.3e} {P1g[i_g]:12.3e}")
print(f"{'P2 at saturation [W]':28s} {P2s[i_s]:12.3e} {P2g[i_g]:12.3e}   (GEM m = 0: {P2g0[i_g]:.3e})")
print(f"{'P3 at saturation [W]':28s} {P3s[i_s]:12.3e} {P3g[i_g]:12.3e}")
print(f"{'P2/P1, P3/P1 at saturation':28s} {P2s[i_s] / P1s[i_s]:5.1e} {P3s[i_s] / P1s[i_s]:6.1e}"
      f"  {P2g[i_g] / P1g[i_g]:5.1e} {P3g[i_g] / P1g[i_g]:6.1e}")
for h in (1, 2, 3):
    print(f"{'bunching b%d at saturation' % h:28s} {b_s[h - 1][i_s]:12.3f} {G[h, -1][1][i_g]:12.3f}")
print(f"{'end: P1, P2, P3 [W]':28s} GENESIS {P1s[-1]:.3e} {P2s[-1]:.3e} {P3s[-1]:.3e}"
      f" | GEM {P1g[-1]:.3e} {P2g[-1]:.3e} {P3g[-1]:.3e}")
# gain lengths of the harmonics in the exponential part (P1 between 1e8 and 1e9 W)
m_s = (P1s > 1e8) & (P1s < 1e9)
m_g = (P1g > 1e8) & (P1g < 1e9)
for name, Ps, Pg in (("P1", P1s, P1g), ("P2", P2s, P2g), ("P3", P3s, P3g)):
    ls = 1 / np.polyfit(z_s[m_s], np.log(Ps[m_s]), 1)[0] if Ps[m_s].min() > 0 else np.nan
    lg = 1 / np.polyfit(z_g[m_g], np.log(Pg[m_g]), 1)[0]
    print(f"{'gain length ' + name + ' [m]':28s} {ls:12.4f} {lg:12.4f}")

fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
for P, c, l in ((P1s, "C0", "h = 1"), (P3s, "C2", "h = 3")):
    ax[0].semilogy(z_s, P, c + "--", lw=1.2)
ax[0].semilogy(z_g, P1g, "C0", label="h = 1")
ax[0].semilogy(z_g, P2g, "C1", label="h = 2 (2 P$_{m=2}$)")
ax[0].semilogy(z_g, np.maximum(P2g0, 1e-3), "C1:", lw=1, label="h = 2, m = 0 (incl. co-moving)")
ax[0].semilogy(z_g, P3g, "C2", label="h = 3")
ax[0].set_ylim(1e2, 3e10)
ax[0].set_xlabel("z [m]"); ax[0].set_ylabel("power [W]")
ax[0].set_title("solid: GEM-PIC-RZ, dashed: GENESIS (h = 2: none)", fontsize=9)
ax[0].legend(fontsize=8, loc="lower right")
for h, c in ((1, "C0"), (2, "C1"), (3, "C2")):
    ax[1].plot(z_s, b_s[h - 1], c + "--", lw=1.2)
    ax[1].plot(G[h, -1][0], G[h, -1][1], c, label=f"b$_{h}$")
ax[1].set_xlabel("z [m]"); ax[1].set_ylabel("bunching factor")
ax[1].legend(fontsize=8)
fig.suptitle("Seeded planar FEL, 10.6 nm, K = 1.5: harmonics", fontsize=10)
fig.tight_layout()
fig.savefig(png, dpi=130)
print("wrote", png)
