"""Checks of the thin-ring self-field test (tests/thin_ring.sh).

shell: independent reference, a long-beam shell model: every macroparticle is a charged, current-carrying cylindrical
shell; E_r and B_phi from the enclosed charge and current (Gauss, Ampere), B_z from the azimuthal current of the
rotating shells outside, the plasma lens; integrated in s with RK2 from the particles of the code's output 0. It has no
E_z (gamma constant), which limits it to distances where the longitudinal space charge of the bunch does not yet
change the matching (here 0.5 mm).
usage: python3 thin_ring.py ring <run dir> <Q pC> <length um> <G T/m> <s_end mm> <rtol> <ur_max> <Lrel_max>
       python3 thin_ring.py static <run dir> <inside_max> <max_max>
  ring:   r percentiles 1/50/99 % of the code within rtol of the shell model, u_r rms < ur_max, angular momentum
          <r u_phi> conserved to Lrel_max
  static: fields of a static or prescribed (gamma = 1e6) ring: max |E_r - B_phi| / max |E_r| inside the ring
          (r < 0.65 um) below inside_max, and over the bunch below max_max, at every output after 0"""
import glob
import os
import sys

import numpy as np

E_, ME, C, EPS0 = 1.602176634e-19, 9.1093837e-31, 2.99792458e8, 8.8541878128e-12
MU0 = 1 / (EPS0 * C * C)


def load(f):
    a = np.loadtxt(f, ndmin=2)
    return a[a[:, 7] == a[:, 7].max()]


def stats(x, y, ux, uy, uz):
    r = np.hypot(x, y)
    up = (x * uy - y * ux) / r
    ur = (x * ux + y * uy) / r
    return np.percentile(r, [1, 50, 99]), (r * up).mean(), ur.std(), np.sqrt(1 + ux ** 2 + uy ** 2 + uz ** 2).mean()


def shell(a, Q, Lb, G, send, ds=2e-6):
    x, y = a[:, 0] * 1e-6, a[:, 1] * 1e-6
    u = np.stack([a[:, 3], a[:, 4], a[:, 5]], 1)
    n = len(x)
    lam = Q / Lb / n

    def accel(x, y, u):
        r = np.hypot(x, y)
        g = np.sqrt(1 + (u * u).sum(1))
        v = u / g[:, None]
        o = np.argsort(r)
        rs = r[o]
        cq = np.arange(1, n + 1) - 0.5
        cz = np.cumsum(v[o, 2]) - 0.5 * v[o, 2]
        kz = lam * ((x * v[:, 1] - y * v[:, 0]) / r)[o] * C / (2 * np.pi * rs)
        outside = np.cumsum(kz[::-1])[::-1] - 0.5 * kz
        Er, Bp, Bz = np.empty(n), np.empty(n), np.empty(n)
        Er[o] = -lam * cq / (2 * np.pi * EPS0 * rs)
        Bp[o] = -MU0 * lam * C * cz / (2 * np.pi * rs)
        Bz[o] = -MU0 * outside
        cx, cy = x / r, y / r
        Bx, By = -Bp * cy + G * y, Bp * cx - G * x  # lens B = G (y, -x): focusing for electrons moving in +z
        vx, vy, vz = v[:, 0], v[:, 1], v[:, 2]
        k = -E_ / (ME * C * C * vz)
        return np.stack([k * (Er * cx + C * (vy * Bz - vz * By)), k * (Er * cy + C * (vz * Bx - vx * Bz)),
                         k * C * (vx * By - vy * Bx)], 1), v

    for _ in range(int(round(send / ds))):
        a1, v1 = accel(x, y, u)
        xm, ym, um = x + 0.5 * ds * v1[:, 0] / v1[:, 2], y + 0.5 * ds * v1[:, 1] / v1[:, 2], u + 0.5 * ds * a1
        a2, v2 = accel(xm, ym, um)
        x, y, u = x + ds * v2[:, 0] / v2[:, 2], y + ds * v2[:, 1] / v2[:, 2], u + ds * a2
    return x * 1e6, y * 1e6, u[:, 0], u[:, 1], u[:, 2]


def ring(d, Q, Lb, G, send, rtol, urmax, lmax):
    files = sorted(glob.glob(os.path.join(d, "beam_*.txt")))
    a0, a1 = load(files[0]), load(files[-1])
    p0, L0, _, _ = stats(a0[:, 0], a0[:, 1], a0[:, 3], a0[:, 4], a0[:, 5])
    pc, Lc, urc, gc = stats(a1[:, 0], a1[:, 1], a1[:, 3], a1[:, 4], a1[:, 5])
    ps, Ls, urs, gs = stats(*shell(a0, Q * 1e-12, Lb * 1e-6, G, send * 1e-3))
    print("r 1/50/99 %%: start %.3f %.3f %.3f | code %.3f %.3f %.3f | shell model %.3f %.3f %.3f um" %
          (*p0, *pc, *ps))
    print("u_r rms: code %.4f, shell model %.4f;  <r u_phi>: start %.4f, code %.4f, shell model %.4f" %
          (urc, urs, L0, Lc, Ls))
    rel = np.abs(pc / ps - 1)
    ok = (rel < rtol).all() and urc < urmax and (abs(L0) < 1e-9 or abs(Lc / L0 - 1) < lmax)
    print("max relative r difference %.4f (limit %.3f), u_r rms %.4f (limit %.4f)%s" %
          (rel.max(), rtol, urc, urmax, "" if abs(L0) < 1e-9 else ", <r u_phi> %.2e (limit %.1e)" % (Lc / L0 - 1, lmax)))
    return ok


def static(d, inside, worst):
    import h5py
    ok = True
    for f in sorted(glob.glob(os.path.join(d, "data_*.h5")))[1:]:
        with h5py.File(f, "r") as H:
            k = float(H.attrs["k0inv_m"]) * 1e6
            zs = H["grid/zeta_slices"][()] * k
            rf = H["grid/r_faces"][()] * k
            Er, Bp = H["fields/Er"][0].real, H["fields/Bphi"][0].real
            sel = (zs > 7.5) & (zs < 8.5)
            E0 = abs(Er[sel]).max()
            L = abs(Er[sel] - Bp[sel])
            li, lm = L[:, rf < 0.65].max() / E0, L[:, rf < 20].max() / E0
            print("step %d: max |E_r - B_phi| / max |E_r| inside the ring %.2e, over the bunch %.2e" %
                  (H.attrs["step"], li, lm))
            ok = ok and li < inside and lm < worst
    return ok


if __name__ == "__main__":
    if sys.argv[1] == "ring":
        ok = ring(sys.argv[2], *[float(v) for v in sys.argv[3:10]])
    else:
        ok = static(sys.argv[2], float(sys.argv[3]), float(sys.argv[4]))
    sys.exit(0 if ok else 1)
