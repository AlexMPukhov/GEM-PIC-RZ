"""Flying focus in gas: per output step, the on-axis laser (amplitude, phase vs the vacuum pattern), the wake E_z
(peak amplitude behind the pulse and the zeta of its first zero crossing behind the pulse) and the density on axis.
usage: python3 check_wake.py out_dir [zeta_c=90]"""
import glob
import re
import sys
import numpy as np
from scipy.signal import hilbert

d = sys.argv[1]
zc = float(sys.argv[2]) if len(sys.argv) > 2 else 90.0
k0 = 2 * np.pi / 0.8
L = 50 * 0.299792458 / np.sqrt(2 * np.log(2))
kap = 1 - np.cos(0.05)
E0 = 4.0124e12  # V/m per unit field (m c omega0 / e at 0.8 um)
print("# z[mm]  laser amp/vac  laser phase-vac[rad]  max|Ez|[GV/m] (zeta < zc-L)  first Ez zero behind pulse [um]"
      "  ne/n0 axis (mean behind pulse)")
for f in sorted(glob.glob(d + "/axis_*.txt")):
    z = float(re.search(r"z = ([0-9.eE+-]+) mm", open(f).readline()).group(1))
    A = np.loadtxt(f)
    zeta, ex, ez, ne = A[:, 0], A[:, 1], A[:, 2], A[:, 3]
    ref = 2 * np.exp(-((zeta - zc) / L) ** 2) * np.cos(k0 * (zeta - zc) - kap * k0 * z * 1e3)
    w = np.abs(zeta - zc) < L
    ha, hr = hilbert(ex), hilbert(ref)
    ph = np.angle(np.sum(ha[w] * np.conj(hr[w])))
    amp = np.sqrt(np.sum(abs(ha[w]) ** 2) / np.sum(abs(hr[w]) ** 2))
    b = zeta < zc - L
    # smooth E_z over a laser period to remove the carrier part
    k = int(round(0.8 / (zeta[1] - zeta[0])))
    ezs = np.convolve(ez, np.ones(k) / k, mode="same")
    zb, eb = zeta[b], ezs[b]
    cross = zb[np.where(np.diff(np.sign(eb)) != 0)[0]]
    z0 = cross[-1] if len(cross) else np.nan
    print(f"{z:6.3f} {amp:8.4f} {ph:+8.4f} {np.abs(eb).max() * E0 / 1e9:9.3f} {z0:9.3f} {ne[b].mean():8.4f}")
