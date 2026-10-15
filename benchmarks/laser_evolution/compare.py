"""Laser evolution benchmark: GEM-PIC-RZ against FBPIC and Wake-T.
Per propagation distance: laser energy W/W0 (fields of the modes m >= 1; Wake-T: from the envelope), peak |E_x| on
the axis, and the spectral centroid <|k|>/k0 (red-shift) of the laser fields.
usage: python3 compare.py [gem_dir ...]   (default: gem_out_a1m3 gem_out_a1m3_dz1 gem_out_th0.6 gem_out)
       -> laser_<gem_dir>.txt (from the HDF5 files with h5_extract), compare.png, and a table on stdout"""
import glob
import os
import subprocess
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

H5 = os.environ.get("H5_EXTRACT", "../../build/h5_extract")
runs = sys.argv[1:] or ["gem_out_a1m3", "gem_out_a1m3_dz1", "gem_out_th0.6", "gem_out"]
labels = {"gem_out": "GEM-PIC-RZ, cold (unstable wake)", "gem_out_th0.6": "GEM-PIC-RZ, cold, theta = 0.6",
          "gem_out_a1m3": "GEM-PIC-RZ, a = 1 (m = 0: 3) k0^-1", "gem_out_a1m3_dz1": "same, dz = 1 um"}


def gem_table(d):
    """z[um] W(m>=1) max|Ex|_axis <|k|>/k0 sigma_k/k0 from the HDF5 and axis files of a GEM-PIC-RZ run"""
    cache = f"laser_{d}.txt"
    files = sorted(glob.glob(f"{d}/data_*.h5"))  # (the HDF5 files are not in the repository)
    if os.path.exists(cache) and (not files or os.path.getmtime(cache) > max(os.path.getmtime(f) for f in files)):
        return np.loadtxt(cache)
    rows = []
    for f in files:
        e = subprocess.run([H5, f, "energy"], capture_output=True, text=True, check=True).stdout.split("\n")[1].split()
        z, W1 = float(e[1]) * 1000, float(e[5])
        sp = np.array([l.split() for l in subprocess.run([H5, f, "spectrum", "-1e9", "1e9", "0.3", "1.6", "261"],
                                                          capture_output=True, text=True, check=True).stdout.split("\n")
                       if l and not l.startswith("#")], float)
        k, S = sp[:, 0], sp[:, 2]
        kbar = np.sum(k * S) / np.sum(S)
        ks = np.sqrt(np.sum((k - kbar) ** 2 * S) / np.sum(S))
        step = int(f.split("_")[-1].split(".")[0])
        ax = np.loadtxt(f"{d}/axis_{step:05d}.txt")
        rows.append((z, W1, np.abs(ax[:, 1]).max(), kbar, ks))
    t = np.array(rows)
    np.savetxt(cache, t, header="z[um]  W(m>=1)[code]  max|Ex|_axis[m c w0/e]  <|k|>/k0  sigma_k/k0")
    return t


fb = np.loadtxt("fbpic_out/laser_fbpic.txt")
fb[:, 1] /= fb[0, 1]
wt = np.loadtxt("waket_out/laser_waket.txt") if os.path.exists("waket_out/laser_waket.txt") else None
plt.rcParams.update({"font.size": 9})
fig, axs = plt.subplots(3, 1, figsize=(6, 7.5), sharex=True)
for j, (lab, col) in enumerate([("W / W0 (laser)", 1), ("max |E_x| on axis [m c ω0/e]", 2), ("<|k|> / k0", 3)]):
    axs[j].plot(fb[:, 0], fb[:, col], "k-", lw=2, label="FBPIC")
    if wt is not None:
        axs[j].plot(wt[:, 0], wt[:, col], "--", color="#d1495b", lw=1.5, label="Wake-T (envelope)")
    axs[j].set_ylabel(lab)
print(f"{'run':16s} {'z':>6s} {'W/W0':>7s} {'FBPIC':>7s} {'max|Ex|':>8s} {'FBPIC':>7s} {'<k>/k0':>7s} {'FBPIC':>7s}")
for d, c in zip(runs, ["#1f4e9c", "#6fa8dc", "#2a9d8f", "#bbbbbb"]):
    if not glob.glob(f"{d}/data_*.h5") and not os.path.exists(f"laser_{d}.txt"):
        continue
    g = gem_table(d)
    g[:, 1] /= g[0, 1]
    for j, col in enumerate([1, 2, 3]):
        axs[j].plot(g[:, 0], g[:, col], "-", color=c, lw=1.2, label=labels.get(d, d))
    for z in (100, 200, 300):
        i, k = np.argmin(abs(g[:, 0] - z)), np.argmin(abs(fb[:, 0] - z))
        if abs(g[i, 0] - z) < 1 and abs(fb[k, 0] - z) < 1:
            print(f"{d:16s} {z:6.0f} {g[i, 1]:7.3f} {fb[k, 1]:7.3f} {g[i, 2]:8.3f} {fb[k, 2]:7.3f} {g[i, 3]:7.4f} {fb[k, 3]:7.4f}")
if wt is not None:
    for z in (100, 200, 300):
        i, k = np.argmin(abs(wt[:, 0] - z)), np.argmin(abs(fb[:, 0] - z))
        if abs(wt[i, 0] - z) < 1 and abs(fb[k, 0] - z) < 1:
            print(f"{'Wake-T':16s} {z:6.0f} {wt[i, 1]:7.3f} {fb[k, 1]:7.3f} {wt[i, 2]:8.3f} {fb[k, 2]:7.3f} {wt[i, 3]:7.4f} {fb[k, 3]:7.4f}")
axs[2].set_xlabel("propagation distance [µm]")
axs[0].legend(fontsize=7)
for a in axs:
    a.grid(alpha=0.3)
fig.tight_layout()
fig.savefig("compare.png", dpi=120)
