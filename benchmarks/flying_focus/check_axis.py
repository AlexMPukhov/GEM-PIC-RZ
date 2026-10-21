"""On-axis E_x of a flying-focus run against the luminal Bessel pulse:
E_x = a0 exp(-(zeta - zc)^2/L^2) cos(k0 (zeta - zc) - kappa k0 s), kappa = 1 - cos(theta).
usage: python3 check_axis.py out_dir a0 theta[rad] tau_fwhm[fs] zeta_c[um] [lambda_um]"""
import glob
import re
import sys
import numpy as np

d, a0, th, tau, zc = sys.argv[1], *map(float, sys.argv[2:6])
lam = float(sys.argv[6]) if len(sys.argv) > 6 else 0.8
k0 = 2 * np.pi / lam
L = tau * 0.299792458 / np.sqrt(2 * np.log(2))  # um, 1/e field half length
kap = 1 - np.cos(th)
print("# z[mm]  max|Ex|  max|Ex_ref|  rms(Ex - ref)/a0  envelope peak position shift [um]")
for f in sorted(glob.glob(d + "/axis_*.txt")):
    z = float(re.search(r"z = ([0-9.eE+-]+) mm", open(f).readline()).group(1))
    A = np.loadtxt(f)
    zeta, ex = A[:, 0], A[:, 1]
    ref = a0 * np.exp(-((zeta - zc) / L) ** 2) * np.cos(k0 * (zeta - zc) - kap * k0 * z * 1e3)
    w = np.abs(zeta - zc) < 2 * L
    env = np.abs(ex) 
    shift = np.sum(zeta[w] * ex[w] ** 2) / np.sum(ex[w] ** 2) - zc
    print(f"{z:7.3f} {np.abs(ex).max():8.4f} {np.abs(ref).max():8.4f} {np.sqrt(np.mean((ex - ref)[w] ** 2)) / a0:9.2e}"
          f" {shift:8.3f}")
