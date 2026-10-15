"""LWFA benchmark: on-axis E_z and E_x of GEM-PIC-RZ against FBPIC at the same propagation distances.
usage: python3 compare.py [run_dir]   ->  compare.png (compare_<run_dir>.png) and the differences on stdout
       run_dir: gem_out (cold plasma, default), gem_out_a2 (plasma.smooth_length = 2 k0^-1), ..."""
import glob
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

RUN = sys.argv[1] if len(sys.argv) > 1 else "gem_out"


def axis_at(D):
    """axis file of RUN at the propagation distance D [um] (header: z = ... mm)"""
    for fn in sorted(glob.glob(f"{RUN}/axis_*.txt")):
        z = float(open(fn).readline().split("z =")[1].split("mm")[0])
        if abs(z * 1000 - D) < 1e-6:
            return fn
    raise SystemExit(f"{RUN}: no axis file at {D} um")


E0 = 4.0126e12  # m c omega0 / e at 0.8 um [V/m]
fig, axs = plt.subplots(2, 1, figsize=(7, 6), sharex=True)
for D, col in [(50, "C0"), (100, "C1")]:
    f = np.loadtxt(f"fbpic_D{D}.txt")
    g = np.loadtxt(axis_at(D))  # zeta[um] Ex Ez ne/n0 ...
    zf, Ezf, Exf = f[:, 0], f[:, 1] / E0, f[:, 2] / E0
    zg, Exg, Ezg = g[:, 0], g[:, 1], g[:, 2]
    # wake region behind the laser: compare E_z (smooth); E_x (laser) over its envelope
    m = (zg > 2.0) & (zg < 26.0)
    Ezi = np.interp(zg[m], zf, Ezf)
    print(f"D = {D} um: max |Ez| GEM {abs(Ezg[m]).max():.4f}, FBPIC {abs(Ezi).max():.4f}; "
          f"rms difference {np.sqrt(np.mean((Ezg[m] - Ezi) ** 2)) / abs(Ezi).max():.3f} of the peak")
    m1 = (zg > 14.0) & (zg < 26.0)  # first period of the wake (behind it: steepened, noisy)
    E1 = np.interp(zg[m1], zf, Ezf)
    print(f"          first period (14 ... 26 um): min Ez GEM {Ezg[m1].min():.4f}, FBPIC {E1.min():.4f}; "
          f"rms difference {np.sqrt(np.mean((Ezg[m1] - E1) ** 2)) / abs(E1).max():.3f} of the peak")
    ml = (zg > 22.0) & (zg < 38.0)
    print(f"          laser: max |Ex| GEM {abs(Exg[ml]).max():.3f}, FBPIC {abs(np.interp(zg[ml], zf, Exf)).max():.3f}")
    axs[0].plot(zg, Ezg, col + "-", lw=1.5, label=f"GEM-PIC-RZ, {D} um")
    axs[0].plot(zf, Ezf, col + "--", lw=1, label=f"FBPIC, {D} um")
    axs[1].plot(zg, Exg, col + "-", lw=0.8, label=f"GEM-PIC-RZ, {D} um")
    axs[1].plot(zf, Exf, col + "--", lw=0.6, label=f"FBPIC, {D} um")
axs[0].set_ylabel("E_z on axis [m c ω0 / e]")
axs[1].set_ylabel("E_x on axis [m c ω0 / e]")
axs[1].set_xlabel("ζ [µm]")
axs[0].set_title("LWFA, a0 = 2, n0 = 3.5e18 cm⁻³: GEM-PIC-RZ vs FBPIC", fontsize=10)
for a in axs:
    a.legend(fontsize=7)
    a.grid(alpha=0.3)
fig.tight_layout()
fig.savefig("compare.png" if RUN == "gem_out" else f"compare_{RUN}.png", dpi=120)
