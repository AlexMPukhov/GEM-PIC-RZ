"""Seeded FEL benchmark: power of GEM-PIC-RZ (steady-state part of the long bunch) against GENESIS 1.3 v4.
usage: python3 compare.py 10nm | 0.8nm [gem output dirs ...]   ->  compare_<case>.png and a table on stdout"""
import sys
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

CASES = {  # GENESIS output, zeta range of the steady-state radiation [um], fit range [m], title
    "10nm": ("genesis/seeded.out.h5", (2.6, 3.6), (0.45, 0.6), "10 nm: 4 kA, 511 MeV, helical K = 1, 1 cm", ["gem_out"]),
    "0.8nm": ("genesis_x/seeded.out.h5", (0.32, 0.36), (2.0, 2.8), "0.8 nm: 4 kA, 2.56 GeV, helical K = 1, 2 cm",
              ["gemx_out"]),
}
case = sys.argv[1] if len(sys.argv) > 1 else "10nm"
gfile, zr, fit, title, dirs = CASES[case]
dirs = sys.argv[2:] or dirs
g = h5py.File(gfile, "r")
zg, Pg = g["Lattice/zplot"][...], g["Field/power"][...].ravel()
# GENESIS with the physics GEM-PIC-RZ has in addition: longitudinal space charge (&efield) and the
# widening of the beam by its transverse space charge (emulated by a beam diverging from a waist)
g2 = h5py.File(gfile.replace("seeded.out", "seeded_sigma_sc.out"), "r")
zg2, Pg2 = g2["Lattice/zplot"][...], g2["Field/power"][...].ravel()


LABELS = {
    "gem_out_nocorr": "GEM-PIC-RZ, before corrections 7 and 8",
    "gem_out": "GEM-PIC-RZ, with correction 7 (undulator kick)",
    "gem_out_v": "GEM-PIC-RZ, with corrections 7 and 8",
    "res_dz_v": "GEM-PIC-RZ, 7 + 8, dz = lambda_u / 32",
    "res_dzeta_v": "GEM-PIC-RZ, 7 + 8, dzeta = lambda / 32",
    "gemx_out_nocorr": "GEM-PIC-RZ, with correction 7 (undulator kick)",
    "gemx_out": "GEM-PIC-RZ, with corrections 7 and 8",
}


def gem(d):
    lines = open(f"{d}/fel_power.txt").read().splitlines()
    zs = np.array([float(x) for x in lines[1].split(":")[1].split()])
    a = np.loadtxt(f"{d}/fel_power.txt")
    sel = (zs >= zr[0] - 1e-9) & (zs < zr[1] - 1e-9)
    return a[:, 1] * 1e-3, a[:, 2:][:, sel].mean(axis=1)


def lgain(z, P):
    m = (z >= fit[0]) & (z <= fit[1])
    return 1.0 / np.polyfit(z[m], np.log(P[m]), 1)[0]


print(f"GENESIS: power gain length {lgain(zg, Pg) * 100:.2f} cm ({fit[0]} ... {fit[1]} m), P_end = {Pg[-1]:.3e} W")
print(f"GENESIS + LSC + beam widening: power gain length {lgain(zg2, Pg2) * 100:.2f} cm, P_end = {Pg2[-1]:.3e} W")
fig, ax = plt.subplots(figsize=(6, 4))
ax.semilogy(zg, Pg, "k-", lw=2, label="GENESIS 1.3 v4 (steady state)")
ax.semilogy(zg2, Pg2, "-", color="0.6", lw=2, label="GENESIS + long. space charge + beam widening")
for d, st in zip(dirs, ["C0--", "C1-.", "C2:", "C3--"]):
    z, P = gem(d)
    print(f"{d}: power gain length {lgain(z, P) * 100:.2f} cm, P_end = {P[-1]:.3e} W")
    for zz in np.linspace(zg[-1] / 3, zg[-1], 5):
        i, j = np.argmin(abs(z - zz)), np.argmin(abs(zg - zz))
        print(f"   z = {zz:.2f} m   GEM {P[i]:.3e} W   GENESIS {Pg[j]:.3e} W   ratio {P[i] / Pg[j]:.3f}")
    ax.semilogy(z, P, st, lw=1.5, label=LABELS.get(d, d))
ax.set_xlabel("z [m]")
ax.set_ylabel("radiation power [W]")
ax.set_title("Seeded FEL, " + title, fontsize=10)
ax.legend(fontsize=8)
ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(f"compare_{case}.png", dpi=120)
