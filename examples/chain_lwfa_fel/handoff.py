"""Hand-off from stage 1 (FBPIC) to stage 2 (GEM-PIC-RZ): where the injected beam sits behind the laser.

Reads the openPMD particle file of stage 1 and stage1/laser_centre.txt, keeps the electrons with gamma > gamma_min
(the same cut as bunch.select of stage 2), prints their charge and moments, and prints the position of their charge
centroid in the window of stage 2: zeta_beam = zeta_laser - (z_laser - z_beam), where zeta_laser is laser.zeta_center
of stage 2 (the laser starts again there; the beam keeps its distance behind the laser centre).

usage: python3 handoff.py <openPMD file> <laser_centre.txt> <gamma_min> <zeta_laser_um>
prints (last line): bunch.zeta_center=<value>um                                        (numpy, h5py)"""
import sys
import h5py
import numpy as np


def record(sp, name, comp):
    d = sp[name][comp]
    v = d[()] if isinstance(d, h5py.Dataset) else np.full(1, d.attrs["value"])
    return v * d.attrs["unitSI"]


path, lc, gmin, zeta_laser = sys.argv[1], sys.argv[2], float(sys.argv[3]), float(sys.argv[4]) * 1e-6
with h5py.File(path, "r") as f:
    it = sorted(f["data"].keys(), key=int)[-1]
    parts = f["data"][it]["particles"]
    sp = parts[list(parts.keys())[0]]
    z = record(sp, "position", "z") + record(sp, "positionOffset", "z")
    x = record(sp, "position", "x") + record(sp, "positionOffset", "x")
    mc = 9.1093837015e-31 * 299792458.0
    ux, uy, uz = (record(sp, "momentum", c) / mc for c in "xyz")
    w = sp["weighting"][()]
z_laser = np.loadtxt(lc)[0]
g = np.sqrt(1 + ux ** 2 + uy ** 2 + uz ** 2)
k = g > gmin
if not k.any():
    sys.exit("handoff: no electron with gamma > %g in %s" % (gmin, path))
w, z, x, g, ux = w[k], z[k], x[k], g[k], ux[k]
W = w.sum()
zc = (w * z).sum() / W
gm = (w * g).sum() / W
rms = lambda v: np.sqrt((w * (v - (w * v).sum() / W) ** 2).sum() / W)
print("stage 1 beam (gamma > %g): %d macroparticles, Q = %.2f pC, E = %.1f MeV, dE/E = %.3f, sigma_z = %.3f um, "
      "sigma_x = %.3f um, eps_nx = %.3f um" % (gmin, k.sum(), W * 1.602176634e-7, (gm - 1) * 0.51099895,
                                              rms(g) / gm, rms(z) * 1e6, rms(x) * 1e6,
                                              np.sqrt(rms(x) ** 2 * rms(ux) ** 2
                                                      - ((w * (x - (w * x).sum() / W) * (ux - (w * ux).sum() / W)).sum() / W) ** 2) * 1e6))
print("laser centre at z = %.3f um, beam centroid %.3f um behind it" % (z_laser * 1e6, (z_laser - zc) * 1e6))
print("bunch.zeta_center=%.4fum" % ((zeta_laser - (z_laser - zc)) * 1e6))
