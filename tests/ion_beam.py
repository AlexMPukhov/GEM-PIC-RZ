"""Write an openPMD 1.1 particle file with ion and proton species for tests/ion_import.sh, as another code would,
and print the reference values computed here from the same particles.

usage: python3 tests/ion_beam.py <file.h5>
prints:  c6w Q[pC] E[MeV] dE/E sig_x[um] sig_y[um] eps_nx[um] eps_ny[um]   (moments of the warm C6+ beam)
         brho <B rho of the cold C6+ beam [T m]>
         ekin <kinetic energy per ion of the cold C6+ beam [MeV]>                                    (numpy, h5py)"""
import sys

import h5py
import numpy as np

QE, ME, CL, U = 1.602176634e-19, 9.1093837015e-31, 299792458.0, 1.66053906660e-27
MP = 1.67262192369e-27
MEC2 = 0.51099895000  # MeV
rng = np.random.default_rng(12345)
f = h5py.File(sys.argv[1], "w")
s = lambda t: np.bytes_(t)
f.attrs["openPMD"] = s("1.1.0")
f.attrs["openPMDextension"] = np.uint32(0)
f.attrs["basePath"] = s("/data/%T/")
f.attrs["particlesPath"] = s("particles/")
f.attrs["iterationEncoding"] = s("fileBased")
f.attrs["iterationFormat"] = s("ions_%T.h5")
it = f.create_group("data/7")
it.attrs["time"], it.attrs["dt"], it.attrs["timeUnitSI"] = 0.0, 1.0, 1.0
parts = it.create_group("particles")


def constant(g, name, value, n):
    c = g.create_group(name)
    c.attrs["value"], c.attrs["shape"], c.attrs["unitSI"] = value, np.array([n], dtype=np.uint64), 1.0
    return c


def species(name, x, y, z, px, py, pz, w, q, m, weighted=False):
    """positions [m], momenta [kg m/s] per real particle, weights; q, m per real particle (scalars or arrays);
    weighted: momentum, charge and mass stored macroWeighted (summed over the macroparticle)"""
    g = parts.create_group(name)
    n = len(x)
    pos, off, mom = g.create_group("position"), g.create_group("positionOffset"), g.create_group("momentum")
    for a, v in zip("xyz", (x, y, z - 0.5)):
        pos.create_dataset(a, data=v).attrs["unitSI"] = 1.0
    for a, v in zip("xyz", (0.0, 0.0, 0.5)):
        constant(off, a, v, n)
    for a, v in zip("xyz", (px, py, pz)):
        mom.create_dataset(a, data=v * w if weighted else v).attrs["unitSI"] = 1.0
    mom.attrs["macroWeighted"], mom.attrs["weightingPower"] = np.uint32(weighted), 1.0
    g.create_dataset("weighting", data=w).attrs["unitSI"] = 1.0
    for rec, v in (("charge", q), ("mass", m)):
        if weighted or np.ndim(v) > 0:
            d = g.create_dataset(rec, data=(np.broadcast_to(v, (n,)) * (w if weighted else 1.0)))
            d.attrs["unitSI"] = 1.0
        else:
            d = constant(g, rec, v, n)
        d.attrs["macroWeighted"], d.attrs["weightingPower"] = np.uint32(weighted), 1.0


def momenta(m, gamma, ux, uy):
    """p = m c u for u_perp given and gamma"""
    uz = np.sqrt(gamma ** 2 - 1 - ux ** 2 - uy ** 2)
    return m * CL * ux, m * CL * uy, m * CL * uz


N = 20000
# warm C6+: 50 pC, gamma 50 with 1 % spread, x-ux correlation, unequal weights; stored macroWeighted
Z, m = 6, 12 * U
x, y, z = rng.normal(0, 80e-6, N), rng.normal(0, 60e-6, N), rng.normal(0, 2e-6, N)
ux = rng.normal(0, 0.0125, N) + 0.002 * x / 80e-6
uy = rng.normal(0, 0.02, N)
gam = 50 * (1 + 0.01 * rng.normal(0, 1, N))
w = rng.uniform(0.5, 1.5, N)
w *= 50e-12 / (Z * QE) / w.sum()
species("c6w", x, y, z, *momenta(m, gam, ux, uy), w, Z * QE, m, weighted=True)
W = w.sum()
mean = lambda v: np.sum(w * v) / W
c2 = lambda a, b: mean((a - mean(a)) * (b - mean(b)))
eps = lambda a, u: np.sqrt(c2(a, a) * c2(u, u) - c2(a, u) ** 2)
g = np.sqrt(1 + ux ** 2 + uy ** 2 + (momenta(m, gam, ux, uy)[2] / (m * CL)) ** 2)
mc2 = m / ME * MEC2
print("c6w %.6e %.6e %.6e %.6e %.6e %.6e %.6e" % (W * Z * QE * 1e12, (mean(g) - 1) * mc2, np.sqrt(c2(g, g)) / mean(g),
      np.sqrt(c2(x, x)) * 1e6, np.sqrt(c2(y, y)) * 1e6, eps(x, ux) * 1e6, eps(y, uy) * 1e6))
# cold C6+: gamma 100, sigma 100 um, 1 fC (no space charge to speak of); charge and mass constant records
gc = 100.0
x, y, z = rng.normal(0, 100e-6, N), rng.normal(0, 100e-6, N), rng.normal(0, 2e-6, N)
species("c6", x, y, z, *momenta(m, np.full(N, gc), np.zeros(N), np.zeros(N)), np.full(N, 1e-15 / (Z * QE) / N),
        Z * QE, m)
print("brho %.10g" % (np.sqrt(gc ** 2 - 1) * m * CL / (Z * QE)))
print("ekin %.10g" % ((gc - 1) * mc2))
# protons of 400 GeV kinetic energy
n = 2000
gp = 1 + 400000.0 / (MP / ME * MEC2)
x, y, z = rng.normal(0, 100e-6, n), rng.normal(0, 100e-6, n), rng.normal(0, 2e-6, n)
species("p", x, y, z, *momenta(MP, np.full(n, gp), np.zeros(n), np.zeros(n)), np.full(n, 1e-15 / QE / n), QE, MP)
# refused: two charge states in one species, and neutral particles
qmix = np.where(np.arange(n) % 2 == 0, 5 * QE, 6 * QE)
species("mixed", x, y, z, *momenta(m, np.full(n, gc), np.zeros(n), np.zeros(n)), np.full(n, 100.0), qmix, m)
species("neutral", x, y, z, *momenta(m, np.full(n, gc), np.zeros(n), np.zeros(n)), np.full(n, 100.0), 0.0, m)
f.close()
