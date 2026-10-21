"""E_x in the plane y = 0 of a flying-focus vacuum run against the incident field, by radius bands.
usage: python3 check_map.py map_file a0 theta tau_fwhm zeta_c [lambda]"""
import re
import sys
import numpy as np
from scipy.special import j0

f, a0, th, tau, zc = sys.argv[1], *map(float, sys.argv[2:6])
lam = float(sys.argv[6]) if len(sys.argv) > 6 else 0.8
k0 = 2 * np.pi / lam
z = float(re.search(r"z = ([0-9.eE+-]+) mm", open(f).readline()).group(1)) * 1e3 * k0  # code units
A = np.loadtxt(f)
zeta, x, ex = A[:, 0] * k0, A[:, 1] * k0, A[:, 2]
L = tau * 0.299792458 / np.sqrt(2 * np.log(2)) * k0
kap = 1 - np.cos(th)
w = np.linspace(1 - 6 / L, 1 + 6 / L, 241)
dw = w[1] - w[0]
c = a0 * dw * L / (2 * np.sqrt(np.pi)) * np.exp(-0.25 * (w - 1) ** 2 * L ** 2)
kr = np.sqrt(np.maximum(0, 2 * w * kap - kap ** 2))
ref = np.zeros(len(zeta))
for j in range(len(w)):
    ref += np.real(c[j] * j0(kr[j] * np.abs(x)) * np.exp(1j * (w[j] * (zeta - zc * k0) - kap * z)))
sel = np.abs(zeta - zc * k0) < 2 * L
r = np.abs(x) / k0
print(f"z = {z / k0 * 1e-3:.3f} mm")
for lo, hi in ((0, 5), (5, 15), (15, 30), (30, 40), (40, 50), (50, 60)):
    m = sel & (r >= lo) & (r < hi)
    print(f"  r {lo:3d}-{hi:3d} um: rms(Ex)/a0 {np.sqrt(np.mean(ex[m] ** 2)) / a0:.4f}  ref {np.sqrt(np.mean(ref[m] ** 2)) / a0:.4f}"
          f"  rms diff/a0 {np.sqrt(np.mean((ex - ref)[m] ** 2)) / a0:.4f}")
