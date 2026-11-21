"""Radiation reaction of the beam list (tests/radiation.sh) against an independent reference.

Helical wiggler in vacuum: |B| = B0 f(z) (f: linear ramps over one period at both ends), E = 0, v ~ z, so
w = |v x B|^2 = B^2 and chi = K gamma B with K = hbar omega0 / (m c^2). Reference (RK4 in z, code units):
   dgamma/dz = -a gamma^2,          a = (2/3) alpha K B^2                       (Landau-Lifshitz)
   dV/dz     = -4 a gamma V + (55 / (24 sqrt 3)) alpha K^2 gamma^4 B^3         (Sokolov-Ternov variance of gamma,
                                                                                 linearized drag feedback)
Checks: classical run <gamma> = reference; stochastic and Gaussian <gamma> = reference within 4 sigma statistics,
rms gamma = sqrt(V) within 4 % (sampling error of 1e5 particles ~1.5 %); no particle gains energy with the
stochastic recoil (pure magnetic field: gamma is constant without radiation), several per cent with the Gaussian
model (Fokker-Planck).
usage: python3 radiation.py <classical dir> <stochastic dir> <gaussian dir> <gamma0> <K_u> <lambda_u [m]> <L [m]> <lambda0 [m]>"""
import glob
import sys

import h5py
import numpy as np

dc, ds, dg = sys.argv[1:4]
g0, Ku, lu, L, l0 = map(float, sys.argv[4:9])
ALPHA = 7.2973525693e-3
k0 = 2 * np.pi / l0
K = 3.8615926796e-13 * k0
B0 = Ku * l0 / lu                          # K k_u / k0, units m c omega0 / e
Lc, luc = L * k0, lu * k0                  # code units


def B(z):
    return B0 * np.clip(np.minimum(z / luc, (Lc - z) / luc), 0.0, 1.0)


def rhs(z, y):
    g, V = y
    b = B(z)
    a = 2 / 3 * ALPHA * K * b * b
    return np.array([-a * g * g, -4 * a * g * V + 55 / (24 * np.sqrt(3)) * ALPHA * K * K * g ** 4 * b ** 3])


n = 20000
h = Lc / n
y = np.array([g0, 0.0])
nph = 0.0
for i in range(n):
    z = i * h
    k1 = rhs(z, y); k2 = rhs(z + h / 2, y + h / 2 * k1); k3 = rhs(z + h / 2, y + h / 2 * k2); k4 = rhs(z + h, y + h * k3)
    y = y + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    nph += 5 / (2 * np.sqrt(3)) * ALPHA * B(z + h / 2) * h
gref, Vref = y
chimax = K * g0 * B0
print(f"  chi = {chimax:.3e}, photons per particle {nph:.2f}; reference: loss {g0 - gref:.2f}, rms {np.sqrt(Vref):.2f}")


def gammas(d):
    fn = sorted(glob.glob(f"{d}/data_*.h5"))[-1]
    with h5py.File(fn) as f:
        p = f["/particles/beam"]
        return np.sqrt(1 + p["ux"][()] ** 2 + p["uy"][()] ** 2 + p["uz"][()] ** 2)


fail = 0
for name, d in (("classical", dc), ("stochastic", ds), ("gaussian", dg)):
    g = gammas(d)
    loss, rms, N = g0 - g.mean(), g.std(), len(g)
    gain = np.mean(g > g0 * (1 + 1e-8))   # 2e-3 in gamma: far below a photon (~0.46 chi gamma = 92); the bunch self-field
    # and the pusher change gamma by ~1e-5
    if name == "classical":
        ok = abs(loss / (g0 - gref) - 1) < 2e-3 and gain == 0
        print(f"  {name:10s}: loss {loss:9.3f} = {loss / (g0 - gref):.5f} x reference, rms {rms:.2e}"
              f"  {'ok' if ok else 'FAILED'}")
    else:
        se = np.sqrt(Vref / N)
        okm = abs(loss - (g0 - gref)) < 4 * se
        okv = abs(rms / np.sqrt(Vref) - 1) < 0.04
        okg = (gain == 0) if name == "stochastic" else (gain > 0.01)
        ok = okm and okv and okg
        print(f"  {name:10s}: loss {loss:9.3f} (reference {g0 - gref:.3f} +- {se:.3f}), rms / Sokolov-Ternov "
              f"{rms / np.sqrt(Vref):.4f}, particles that gained energy {gain:.4f}  {'ok' if ok else 'FAILED'}")
    fail += not ok
sys.exit(1 if fail else 0)
