#!/usr/bin/env python3
"""Reference for the undulator tests: coherent emission of a Gaussian bunch in a
planar undulator from the far-field radiation integral (Jackson 14.65),

    d2W/dw dOmega = e^2/(16 pi^3 eps0 c) |int n x ((n - b) x db/dt)/(1 - n.b)^2 e^{iw(t - n.r/c)} dt|^2,

for the trajectory in the field of GEM-PIC-RZ (B_y = B0 g(z) cos(k_u z), linear
ramps of whole periods, no natural focusing), times N_e^2 and the form factor
exp(-k^2 sigma_z^2 - k^2 theta^2 sigma_r^2). Prints the energy in J and in code
units (h5_extract energy, column W(m>=1)) and writes the spectrum dW/dk.

usage: python3 tools/undulator_radint.py [periods]   (defaults: tests/undulator.deck, 10 periods)
needs numpy; about 1 minute per 10 periods
"""
import sys
import numpy as np

lam0 = 0.8e-6                     # reference wavelength (code units)
gamma, K = 100.0, 0.5
period = 14.2222e-3               # m
periods = int(sys.argv[1]) if len(sys.argv) > 1 else 10
ramp_periods = 1
sigma_z, sigma_r = 0.1e-6, 100e-6  # m
charge = 1e-15                    # C
nu = np.linspace(0.4, 1.6, 1201)  # k / k0
theta = np.linspace(0.0, 8e-3, 161)
phi = np.linspace(0.0, 2 * np.pi, 32, endpoint=False)

c, e, eps0, me = 299792458.0, 1.602176634e-19, 8.8541878128e-12, 9.1093837015e-31
k0 = 2 * np.pi / lam0
lu = period * k0; ku = 2 * np.pi / lu; L = periods * lu; ramp = ramp_periods * lu; B0 = K * ku
sz, sr, Ne = sigma_z * k0, sigma_r * k0, charge / e
# trajectory (code units: lengths 1/k0, c = 1): du_x/dz = B_y
z = np.linspace(0.0, L, periods * 96 + 1)
g = np.minimum(1.0, np.minimum(z / ramp, (L - z) / ramp)) if ramp > 0 else np.ones_like(z)
By = B0 * g * np.cos(ku * z)
cum = lambda f: np.concatenate([[0.0], np.cumsum(0.5 * (f[1:] + f[:-1]) * np.diff(z))])
ux = cum(By); uz = np.sqrt(gamma**2 - 1 - ux**2)
bx, bz, dtdz = ux / gamma, uz / gamma, gamma / uz
t, x = cum(dtdz), cum(bx * dtdz)
dbx, dbz = np.gradient(bx, z) / dtdz, np.gradient(bz, z) / dtdz
wt = dtdz * np.gradient(z)
spec = np.zeros(len(nu))
for th in theta:
    for p in phi:
        n = np.array([np.sin(th) * np.cos(p), np.sin(th) * np.sin(p), np.cos(th)])
        a0 = n[1] * dbz; a1 = (n[2] - bz) * dbx - (n[0] - bx) * dbz; a2 = -n[1] * dbx
        A = [n[1] * a2 - n[2] * a1, n[2] * a0 - n[0] * a2, n[0] * a1 - n[1] * a0]
        w = wt / (1 - n[0] * bx - n[2] * bz) ** 2
        E = np.exp(1j * np.outer(nu, t - (n[0] * x + n[2] * z)))
        d2 = sum(np.abs(E @ (Ai * w)) ** 2 for Ai in A)
        F2 = np.exp(-nu**2 * (sz**2 + (np.sin(th) * sr) ** 2))
        spec += d2 * F2 * np.sin(th) * (theta[1] - theta[0]) * (2 * np.pi / len(phi))
unit = e**2 * (c * k0) / (16 * np.pi**3 * eps0 * c)        # J per (d nu dOmega)
E0 = me * c**2 * k0 / e; Wcode = eps0 * E0**2 / k0**3       # J per code unit of energy
W = Ne**2 * unit * spec.sum() * (nu[1] - nu[0])
print("periods %d: lambda_r = %.4f um, slippage in the undulator %.4f um" %
      (periods, period * (1 + K * K / 2) / (2 * gamma**2) * 1e6, (t[-1] - z[-1]) / k0 * 1e6))
print("coherent energy %.4e J = %.4e code units; peak of dW/dk at k = %.4f k0" %
      (W, W / Wcode, nu[np.argmax(spec)]))
np.savetxt("undulator_spectrum.txt", np.c_[nu, Ne**2 * unit * spec / Wcode], header="k/k0  dW/dk (code units)")
