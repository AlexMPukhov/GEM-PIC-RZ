#!/usr/bin/env python3
"""Corrugated metal pipe as an equivalent thin dielectric layer (tests/corrugated.sh, examples/corrugated_pipe).

Small corrugations (K. Bane and G. Stupakov, PRST-AB 6, 024401 (2003)): bore radius a, depth h, period p, gap g,
h, p << a. Each slot is a short-circuited parallel-plate line: averaged over a period the wall has the surface
impedance Z_s = -i Z0 k h g / p. A dielectric layer of thickness h on the metal has Z_s = -i Z0 k h (eps - 1) / eps,
so the corrugation acts on the synchronous mode like the layer

    a < r < a + h,   eps = p / (p - g),   metal wall at a + h.

Leading order in h/a (surface-impedance boundary condition on the bore field E_z = E0 I0(tau r)):
    k0 = sqrt(2 p / (a h g)),   1 - beta_g = 4 h g / (a p),   kappa = Z0 c / (2 pi a^2)   (W(0+) = 2 kappa),
and the on-axis wake of a Gaussian bunch Q, sigma behind it: E_z = 2 kappa Q exp(-k0^2 sigma^2 / 2) cos(k0 xi).

  python3 corrugated.py theory                       # the layer (exact modal theory) vs the formulas, h/a -> 0
  python3 corrugated.py check <axis file> a h p g Q[pC] sigma[um] zeta_b[um] s[um] lambda0[um] tol_k tol_a
      on-axis wake of a GEM-PIC-RZ run of the equivalent layer (lengths in um) against the modal theory of the
      layer (the model the code solves) and against the corrugated-pipe formulas; exit 1 on failure
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dielectric import modes, group_velocity, eunit  # noqa: E402

Z0C = 376.730313668 * 299792458.0


def bane_stupakov(a, h, p, g):
    """leading-order corrugated pipe (SI lengths): k0 [1/m], beta_g, kappa [V/C/m]"""
    return np.sqrt(2 * p / (a * h * g)), 1.0 - 4 * h * g / (a * p), Z0C / (2 * np.pi * a * a)


def eps_of(p, g):
    return p / (p - g)


def theory():
    a, p, g = 1e-3, 230e-6, 115e-6
    print("equivalent layer eps = p/(p - g) = %.3f; bore a = 1 mm, p = 230 um, g = 115 um" % eps_of(p, g))
    print("   h [um]   k0 layer / formula    1-beta_g layer / formula   kappa [kV/pC/m] layer / formula")
    for h in (240e-6, 120e-6, 60e-6, 30e-6, 15e-6):
        k, bg, kap = modes(a, a + h, eps_of(p, g), 1)[0]
        k0, bg0, kap0 = bane_stupakov(a, h, p, g)
        print("   %5.0f    %8.4f / %8.4f 1/mm   %.4f / %.4f            %.3f / %.3f" %
              (h * 1e6, k * 1e-3, k0 * 1e-3, 1 - bg, 1 - bg0, kap * 1e-15, kap0 * 1e-15))
    h = 60e-6
    k, bg, kap = modes(a, a + h, eps_of(p, g), 1)[0]
    print("h = 60 um: lambda = %.1f um (%.0f GHz); GEM system beta_g = %.4f (Maxwell %.4f)" %
          (2 * np.pi / k * 1e6, 299792458.0 * k / (2 * np.pi) * 1e-9, group_velocity(k, a, a + h, eps_of(p, g),
                                                                                   "reduced"), bg))


def check(args):
    path = args[0]
    a, h, p, g, Q, sig, zb, s, lam0, tolk, tola = map(float, args[1:12])
    eps = eps_of(p, g)
    md = modes(a * 1e-6, (a + h) * 1e-6, eps, 1)
    k1, bg1, kap1 = md[0]
    k0, bg0, kap0 = bane_stupakov(a * 1e-6, h * 1e-6, p * 1e-6, g * 1e-6)
    bgr = group_velocity(k1, a * 1e-6, (a + h) * 1e-6, eps, "reduced")
    E1 = eunit(lam0)
    d = np.loadtxt(path)
    z, Ez = d[:, 0], d[:, 2]
    fill = s * (1 / bgr - 1)  # the filled part of the wake (the code solves the reduced system)
    sel = (z < zb - 5 * sig) & (z > zb - 0.8 * fill) & (z > z.min() + 5)
    if np.count_nonzero(sel) < 50:
        print("check: the filled wake is not inside the window")
        return 1
    zz, ee = (z[sel] - zb) * 1e-6, Ez[sel] * E1  # m, V/m

    def fit(k):
        A = np.column_stack([np.cos(k * zz), np.sin(k * zz)])
        c, *_ = np.linalg.lstsq(A, ee, rcond=None)
        return c, float(np.sum((A @ c - ee) ** 2))

    ks = np.linspace(0.9 * k1, 1.1 * k1, 401)  # scan, then golden-section refinement of the wavenumber
    kb = ks[np.argmin([fit(k)[1] for k in ks])]
    lo, hi = kb - 0.2 * k1 / 400 * 2, kb + 0.2 * k1 / 400 * 2
    for _ in range(80):
        m1, m2 = lo + 0.382 * (hi - lo), lo + 0.618 * (hi - lo)
        if fit(m1)[1] < fit(m2)[1]:
            hi = m2
        else:
            lo = m1
    kc = 0.5 * (lo + hi)
    c, res = fit(kc)
    amp_c = c[0]
    F = np.exp(-0.5 * (k1 * sig * 1e-6) ** 2)
    amp_l = Q * 1e-12 * 2 * kap1 * F
    amp_b = Q * 1e-12 * 2 * kap0 * np.exp(-0.5 * (k0 * sig * 1e-6) ** 2)
    rms = np.sqrt(res / len(ee)) / abs(amp_c)
    print("corrugated pipe a = %g um, h = %g um, p = %g um, g = %g um as a layer eps = %.4f, %g < r < %g um" %
          (a, h, p, g, eps, a, a + h))
    print("  wavenumber: code %.5f 1/mm, layer modal theory %.5f (%+.2e), Bane-Stupakov %.5f (%+.2e)" %
          (kc * 1e-3, k1 * 1e-3, kc / k1 - 1, k0 * 1e-3, kc / k0 - 1))
    print("  amplitude:  code %.4f MV/m, layer %.4f MV/m (%+.2e), Bane-Stupakov %.4f MV/m (%+.2e);  sin part %.1e,"
          " rms rest %.1e" % (amp_c * 1e-6, amp_l * 1e-6, amp_c / amp_l - 1, amp_b * 1e-6, amp_c / amp_b - 1,
                              c[1] / amp_c, rms))
    print("  loss factor: layer %.3f kV/pC/m, Bane-Stupakov Z0 c/(2 pi a^2) = %.3f kV/pC/m;  beta_g: layer Maxwell "
          "%.4f, GEM system %.4f, formula %.4f" % (kap1 * 1e-15, kap0 * 1e-15, bg1, bgr, bg0))
    ok = abs(kc / k1 - 1) < tolk and abs(amp_c / amp_l - 1) < tola
    print("  code vs the modal theory of the layer: %s" % ("ok" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "theory":
        theory()
    elif len(sys.argv) > 1 and sys.argv[1] == "check":
        sys.exit(check(sys.argv[2:]))
    else:
        print(__doc__)
        sys.exit(1)
