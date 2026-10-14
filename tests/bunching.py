"""Bunching noise of the beam list of a GEM-PIC-RZ data file in Hann windows.
usage: python3 bunching.py <data.h5> <wavelength_m> <window in wavelengths>
prints <|b|^2> N_e and <|b|^2> N_macro (1 for independent electrons / macroparticles)."""
import sys

import h5py
import numpy as np

f = h5py.File(sys.argv[1], "r")
L, Q0 = f.attrs["k0inv_m"], f.attrs["charge_unit_C"]
z = f["particles/beam/zeta"][...] * L
q = np.abs(f["particles/beam/charge"][...]) * Q0
lam = float(sys.argv[2])
k = 2 * np.pi / lam
W = int(sys.argv[3]) * lam
b2, ne, nm = [], [], []
for a in np.arange(z.min() + lam, z.max() - W - lam, W):
    s = (z >= a) & (z < a + W)
    w = np.sin(np.pi * (z[s] - a) / W) ** 2 * q[s]
    b2.append(abs(np.sum(w * np.exp(1j * k * z[s])) / np.sum(w)) ** 2)
    # independent electrons: sum(w^2) / sum(w)^2 = <sin^4> / <sin^2>^2 / N = 1.5 / N
    ne.append(np.sum(q[s]) / 1.602176634e-19 / 1.5)
    nm.append(s.sum() / 1.5)
print(np.mean(b2), np.mean(ne), np.mean(nm))
