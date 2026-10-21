"""Figure of the flying-focus gas test: maps of E_x and n_e at the end, on-axis E_z at several z, laser amplitude."""
import glob
import re
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

d = sys.argv[1] if len(sys.argv) > 1 else "gas_out"
E0 = 4.0124e12
maps = sorted(glob.glob(d + "/map_*.txt"))
M = np.loadtxt(maps[-1])
zeta = np.unique(M[:, 0]); nz = len(zeta); nx = len(M) // nz
M = M.reshape(nz, nx, M.shape[1]); x = M[0, :, 1]
zlab = re.search(r"z = ([0-9.eE+-]+) mm", open(maps[-1]).readline()).group(1)
fig, ax = plt.subplots(2, 2, figsize=(10, 6.5))
ex = M[:, :, 2]
im = ax[0, 0].pcolormesh(zeta, x, ex.T, cmap="RdBu_r", vmin=-0.6, vmax=0.6, shading="auto")
ax[0, 0].set_title(f"E_x (a0 units, clipped +-0.6), z = {zlab} mm", fontsize=9)
plt.colorbar(im, ax=ax[0, 0])
im = ax[0, 1].pcolormesh(zeta, x, M[:, :, 4].T, cmap="viridis", vmin=0, vmax=3, shading="auto")
ax[0, 1].set_title("n_e / n0", fontsize=9)
plt.colorbar(im, ax=ax[0, 1])
for a in ax[0]:
    a.set_xlabel("zeta (um)"); a.set_ylabel("x (um)")
cols = plt.cm.viridis(np.linspace(0, 1, 6))
files = sorted(glob.glob(d + "/axis_*.txt"))
for c, f in zip(cols, files[2::2]):
    A = np.loadtxt(f)
    k = int(round(0.8 / (A[1, 0] - A[0, 0])))
    ezs = np.convolve(A[:, 2], np.ones(k) / k, mode="same")
    zz = re.search(r"z = ([0-9.eE+-]+) mm", open(f).readline()).group(1)
    ax[1, 0].plot(A[:, 0], ezs * E0 / 1e9, color=c, lw=0.9, label=f"z = {float(zz):.1f} mm")
ax[1, 0].set_xlabel("zeta (um)"); ax[1, 0].set_ylabel("E_z on axis (GV/m), averaged over 0.8 um")
ax[1, 0].legend(fontsize=7)
ax[1, 0].set_title("wake on the axis (laser centre at zeta = 90 um)", fontsize=9)
zs, amp = [], []
for f in files:
    A = np.loadtxt(f)
    zs.append(float(re.search(r"z = ([0-9.eE+-]+) mm", open(f).readline()).group(1)))
    amp.append(np.abs(A[:, 1]).max())
ax[1, 1].plot(zs, amp, "o-")
ax[1, 1].set_ylim(0, 2.5); ax[1, 1].set_xlabel("z (mm)"); ax[1, 1].set_ylabel("max |E_x| on axis (a0 units)")
ax[1, 1].set_title("laser on the axis (vacuum: 2)", fontsize=9)
fig.tight_layout()
fig.savefig(d + ".png", dpi=110)
print("wrote", d + ".png")
