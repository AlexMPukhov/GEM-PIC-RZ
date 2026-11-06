"""Radiation reaction (tests/radiation.sh), the undulator correction and the thin chicane.

undulator <classical dir> <stochastic dir> <gamma0> <K> <lambda_u [m]> <periods> <ramp periods>
  planar undulator with K ~ 1: the mean loss is the Larmor loss (1/3) r_e gamma^2 k_u^2 K(z)^2 per length, the
  variance d<dgamma^2>/dz = (7/15) r_e lambdabar_C gamma^4 k_u^3 K^2 F(K) of the radiation integral over all
  harmonics (tools/undulator_diffusion.py: F = 1.2032 K + 1/(1 + 1.34 K + 0.43 K^2), checked there to 0.8 %), i.e.
  1.29 x the local constant-field value at K = 1. Linearized drag feedback as in radiation.py (negligible here).
chicane <classical dir> <stochastic dir> <gaussian dir> <gamma0> <theta> <L [m]> <dipoles>
  B = p theta / (e L), chi = gamma B / B_S, path n L; Lambda = (5/(2 sqrt3)) alpha gamma theta n (photons),
  mean loss (2/3) r_e gamma^2 n L / rho^2 ... in gamma units: Lambda m1, variance Lambda m2 with the synchrotron
  moments m1 = (4 sqrt3/15) chi gamma, m2 = (11/12) chi^2 gamma^2; P(no loss) = exp(-Lambda) exactly. The loss of
  each particle is measured against the classical run, which removes exactly the code's Lambda m1 from the same
  state (particles matched by position): gamma - gamma_classical is a sharp spike at Lambda m1 for the particles
  without a photon, loss = spike - (gamma - gamma_classical); the self-field noise of gamma (~1e-5) drops out."""
import glob
import sys

import h5py
import numpy as np

RE, LC, ALPHA = 2.8179403262e-15, 3.8615926796e-13, 7.2973525693e-3


def gammas(d, order=False):
    fn = sorted(glob.glob(f"{d}/data_*.h5"))[-1]
    with h5py.File(fn) as f:
        p = f["/particles/beam"]
        g = np.sqrt(1 + p["ux"][()] ** 2 + p["uy"][()] ** 2 + p["uz"][()] ** 2)
        if order:  # sorted by position (the same particles in every run)
            g = g[np.lexsort((p["y"][()], p["x"][()]))]
        return g


def verdict(ok):
    return "ok" if ok else "FAILED"


fail = 0
if sys.argv[1] == "undulator":
    dc, ds = sys.argv[2:4]
    g0, K, lu, nper, nr = map(float, sys.argv[4:9])
    ku, L = 2 * np.pi / lu, nper * lu
    z = np.linspace(0, L, 200001)
    Kz = K * np.clip(np.minimum(z / (nr * lu), (L - z) / (nr * lu)), 0, 1)
    F = 1.2032 * Kz + 1 / (1 + 1.34 * Kz + 0.43 * Kz * Kz)
    g, V = g0, 0.0
    h = z[1] - z[0]
    for k in Kz[:-1]:                                     # Euler is enough: loss ~1e-6 of gamma
        a = RE * ku * ku * k * k / 3
        V += h * (-4 * a * g * V + 7 / 15 * RE * LC * g ** 4 * ku ** 3 * k * k * (1.2032 * k + 1 / (1 + 1.34 * k + 0.43 * k * k)))
        g -= h * a * g * g
    Vloc = np.trapezoid(0.5615 * RE * LC * g0 ** 4 * ku ** 3 * Kz ** 3, z)
    gc, gs = gammas(dc), gammas(ds)
    lc, ls, rms = g0 - gc.mean(), g0 - gs.mean(), (gs - gc).std()
    se = np.sqrt(V / len(gs))
    print(f"  planar undulator K = {K:g}: loss classical {lc:.5f} (Larmor {g0 - g:.5f}), stochastic {ls:.5f} +- {se:.5f};"
          f" rms(stochastic - classical) {rms:.5f} = {rms / np.sqrt(V):.4f} x radiation integral"
          f" (local-field value would be {np.sqrt(Vloc / V):.4f})")
    ok = abs(lc / (g0 - g) - 1) < 2e-3 and abs(ls - (g0 - g)) < 4 * se and abs(rms / np.sqrt(V) - 1) < 0.03
    print(f"  {verdict(ok)}")
    fail += not ok
else:
    dc, ds, dg = sys.argv[2:5]
    g0, th, L, n = map(float, sys.argv[5:9])
    rho = L / th * np.sqrt(1 - 1 / g0 ** 2)                # bending radius of the reference (= the particles)
    chi = g0 * g0 * LC / rho
    lam = 5 / (2 * np.sqrt(3)) * ALPHA * g0 * th * n
    m1, m2 = 4 * np.sqrt(3) / 15 * chi * g0, 11 / 12 * chi * chi * g0 * g0
    gc = gammas(dc, True)
    for name, d in (("stochastic", ds), ("gaussian", dg)):
        dd = gammas(d, True) - gc     # = Lambda m1 (the code's value) for the particles without a photon
        near = dd[np.abs(dd - lam * m1) < 1e-2 * lam * m1]
        c = np.median(near) if len(near) else lam * m1
        dl = c - dd
        N = len(dl)
        p0 = np.mean(np.abs(dl) < 1e-6 * m1)
        se = np.sqrt(lam * m2 / N)
        okm = abs(dl.mean() - lam * m1) < 4 * se
        okv = abs(dl.var() / (lam * m2) - 1) < 0.04
        if name == "stochastic":
            ok0 = abs(p0 - np.exp(-lam)) < 4 * np.sqrt(np.exp(-lam) / N) and dl.min() > -1e-6 * m1
        else:
            ok0 = np.mean(dl < 0) > 0.01
        ok = okm and okv and ok0
        print(f"  chicane {name:10s}: Lambda {lam:.3f}, chi {chi:.2e}; mean loss {dl.mean():.5f} (Lambda m1 {lam * m1:.5f}"
              f" +- {se:.5f}), var / (Lambda m2) {dl.var() / (lam * m2):.4f}, P(no loss) {p0:.4f} (exp(-Lambda)"
              f" {np.exp(-lam):.4f}), energy gain {np.mean(dl < -1e-6 * m1):.4f}  {verdict(ok)}")
        fail += not ok
sys.exit(1 if fail else 0)
