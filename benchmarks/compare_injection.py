"""Injection benchmarks (ionization_fbpic/, downramp_fbpic/): injected charge, energy and spectrum, GEM-PIC-RZ against FBPIC.
GEM-PIC-RZ: particle files <run>/beam_<step>.txt (trapped electrons on the beam list: species -30 = electrons from the
ionization of ion block 0, -1 = self-trapped plasma electrons). FBPIC: fbpic_out/elec*_<ct>.txt (all electrons from
nitrogen, or all plasma electrons with u_z > 1). Both are compared above a kinetic-energy cut, so that the electrons
in the sheath that FBPIC holds but GEM-PIC-RZ streams through the window are left out.
Position: GEM-PIC-RZ's level s is the plane z = s (GT2), and its beam particles are at the lab position z = s. An FBPIC
output at time t holds the electrons at z = ct - (zeta_laser - zeta), i.e. a trapped bunch ~17 um behind the laser
is at z = ct - 17 um. The FBPIC rows are therefore plotted at the charge-weighted lab position of the electrons
above 5 MeV (column z_beam), not at ct.
usage (in the benchmark directory): python3 ../compare_injection.py ionization|downramp [gem_run ...]
       -> injection.txt, injection.png and a table on stdout"""
import glob
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

mode = sys.argv[1]
runs = sys.argv[2:] or ["gem_out"]
species = -30 if mode == "ionization" else -1
fpat = "fbpic_out/elec_N_*.txt*" if mode == "ionization" else "fbpic_out/elec_[0-9]*.txt*"  # (.gz in the repository)
CUTS = (2.0, 5.0)  # MeV


ZETA_LASER = 30.56  # um from the rear of the window


def stats(x, y, zeta, ux, uy, uz, q):
    E = (np.sqrt(1 + ux ** 2 + uy ** 2 + uz ** 2) - 1) * 0.511
    m = E > CUTS[-1]
    out = [np.sum(q[m] * zeta[m]) / q[m].sum() - ZETA_LASER if q[m].sum() > 0 else np.nan]  # beam behind the laser
    for cut in CUTS:
        m = E > cut
        Q = q[m].sum()
        Em = np.sum(q[m] * E[m]) / Q if Q > 0 else 0.0
        dE = np.sqrt(np.sum(q[m] * (E[m] - Em) ** 2) / Q) if Q > 0 else 0.0
        out += [Q, Em, dE]
    return out, E


def gem(run):
    rows, last = [], None
    for f in sorted(glob.glob(f"{run}/beam_*.txt*")):
        z = float(open(f.replace("beam_", "axis_").replace(".gz", "")).readline().split("z =")[1].split("mm")[0]) * 1000
        a = np.loadtxt(f, ndmin=2)
        a = a[a[:, 7] == species] if a.size else np.zeros((0, 8))
        s, E = stats(a[:, 0], a[:, 1], a[:, 2], a[:, 3], a[:, 4], a[:, 5], np.abs(a[:, 6]))
        rows.append([z, z] + s[1:])  # level s = lab position of the beam particles
        last = (E, np.abs(a[:, 6]))
    return np.array(rows), last


def fbpic():
    rows, last = [], None
    for f in sorted(glob.glob(fpat)):
        ct = float(f.split("_")[-1].replace(".gz", "")[:-4])
        a = np.loadtxt(f, ndmin=2)
        if a.size == 0:
            a = np.zeros((0, 7))
        s, E = stats(a[:, 0], a[:, 1], a[:, 2], a[:, 3], a[:, 4], a[:, 5], np.abs(a[:, 6]))
        rows.append([ct, ct + s[0]] + s[1:])  # ct, lab position of the beam (z_beam), ...
        last = (E, np.abs(a[:, 6]))
    return np.array(rows), last


fb, fbl = fbpic()
hdr = "s_or_ct[um]  z_beam[um]  " + "  ".join(f"Q(E>{c:g})[pC] <E>[MeV] dE_rms[MeV]" for c in CUTS)
np.savetxt("injection_fbpic.txt", fb, header=hdr)
fig, axs = plt.subplots(3, 1, figsize=(6, 8))
axs[0].plot(fb[:, 1], fb[:, 5], "k-", lw=2, label="FBPIC")
axs[1].plot(fb[:, 1], fb[:, 6], "k-", lw=2, label="FBPIC")
bins = np.linspace(CUTS[0], max(30, fbl[0].max() * 1.05 if fbl[0].size else 30), 57)
axs[2].hist(fbl[0], bins, weights=fbl[1] / np.diff(bins)[0], histtype="step", color="k", lw=2, label="FBPIC (last output)")
print(f"{'run':10s} {'ct|s':>6s} {'z_beam':>6s} " + " ".join(f"{'Q>' + str(c):>8s} {'<E>':>6s} {'dE':>6s}" for c in CUTS))
for r in fb:
    print(f"{'FBPIC':10s} {r[0]:6.0f} {r[1]:6.1f} " + " ".join(f"{r[2 + 3 * k]:8.3f} {r[3 + 3 * k]:6.2f} {r[4 + 3 * k]:6.2f}" for k in range(len(CUTS))))
for run, col in zip(runs, ["#1f4e9c", "#6fa8dc", "#2a9d8f"]):
    g, gl = gem(run)
    np.savetxt(f"injection_{run}.txt", g, header=hdr)
    for r in g:
        print(f"{run[:10]:10s} {r[0]:6.0f} {r[1]:6.1f} " + " ".join(f"{r[2 + 3 * k]:8.3f} {r[3 + 3 * k]:6.2f} {r[4 + 3 * k]:6.2f}" for k in range(len(CUTS))))
    axs[0].plot(g[:, 1], g[:, 5], "-", color=col, label=run)
    axs[1].plot(g[:, 1], g[:, 6], "-", color=col, label=run)
    axs[2].hist(gl[0], bins, weights=gl[1] / np.diff(bins)[0], histtype="step", color=col, label=run)
axs[0].set_ylabel(f"Q (E > {CUTS[1]:g} MeV) [pC]")
axs[1].set_ylabel(f"<E> (E > {CUTS[1]:g} MeV) [MeV]")
axs[0].set_xlabel("lab position of the bunch [µm]")
axs[1].set_xlabel("lab position of the bunch [µm]")
axs[2].set_xlabel("kinetic energy [MeV] (last output of each run)")
axs[2].set_ylabel("dQ/dE [pC/MeV]")
for a in axs:
    a.legend(fontsize=7)
    a.grid(alpha=0.3)
fig.tight_layout()
fig.savefig("injection.png", dpi=120)
