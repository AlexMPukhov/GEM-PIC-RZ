"""Reference of the metal-wall test (tests/wall.sh): a laser in a perfectly conducting circular waveguide (vacuum).

Independent reference: modal theory of the guide. The initial field of the code (map_00000.txt, linear polarization
along x: E_r = f(r) cos phi, E_phi = -f(r) sin phi) is decomposed, for every Fourier component kappa in zeta, into
the m = 1 modes of a perfectly conducting tube of radius a:
  TM_1n: E_perp = grad(J1(k r) cos phi),       k = j_1n / a   (E_z = 0 at the wall),
  TE_1n: E_perp = z x grad(J1(k r) sin phi),   k = j'_1n / a  (E_phi = 0 at the wall).
A mode exp(i kappa zeta - i omega s) evolves with the dispersion relation of the scheme in vacuum (paper, Eq. disp),
  tan(omega Delta / 2) = Delta h k^2 sin(theta) / (16 sin^2(theta/2) - h^2 k^2 cos^2(theta/2)),   theta = kappa h,
so that the comparison isolates the radial boundary from the zeta and s discretization (option "continuum": the
limit omega = 2 kappa k^2 / (4 kappa^2 - k^2), to show the convergence). In a slowly tapered tube a(s) every mode
keeps its normalized amplitude (adiabatic local modes, no conversion) and collects the phase of the steps with
the radius of their new level, a(s^1). E_z of a TM mode follows from div E = 0, E_z = i div(E_perp) / kappa_e,
with kappa_e = (2/h) tan(theta/2) of the trapezoidal rule in zeta (kappa in the continuum).

usage: python3 wall.py <run dir> <final map file> <lambda0 um> <s_end um> <dz um> <a_start um> <a_end um> <tol Ex> <tol Ez>
       [continuum]
  a(s) is linear between a_start (s = 0) and a_end (s = s_end), the same law as the deck.
  The code is compared with the reference on r < a_end - 2 dr; the rms difference over the window, relative to the
  maximum of the reference, must be below tol (E_x and E_z)."""
import sys

import numpy as np

TAU = np.linspace(0.0, 2.0 * np.pi, 513)[:-1]


def jn(n, x):
    """Bessel J_n(x) = (1/2pi) int_0^2pi cos(n t - x sin t) dt (trapezoid: exponentially convergent)."""
    x = np.asarray(x, dtype=float)
    return np.cos(n * TAU[None, :] - x.reshape(-1, 1) * np.sin(TAU)[None, :]).mean(axis=1).reshape(x.shape)


def j1p(x):
    return 0.5 * (jn(0, x) - jn(2, x))


def zeros(f, n, step=0.05):
    """first n positive zeros of f (scan and bisection)"""
    out, x0 = [], step
    f0 = f(np.array([x0]))[0]
    while len(out) < n:
        x1 = x0 + step
        f1 = f(np.array([x1]))[0]
        if f0 * f1 < 0:
            a, b, fa = x0, x1, f0
            for _ in range(60):
                c = 0.5 * (a + b)
                fc = f(np.array([c]))[0]
                if fa * fc <= 0:
                    b = c
                else:
                    a, fa = c, fc
            out.append(0.5 * (a + b))
        x0, f0 = x1, f1
    return np.array(out)


def load(path):
    d = np.loadtxt(path)
    d = d[d[:, 1] > 0]  # x > 0: phi = 0
    zs = np.unique(d[:, 0])
    rs = np.unique(d[:, 1])
    ex = np.zeros((len(zs), len(rs)))
    ez = np.zeros((len(zs), len(rs)))
    iz = np.searchsorted(zs, d[:, 0])
    ir = np.searchsorted(rs, d[:, 1])
    ex[iz, ir] = d[:, 2]
    ez[iz, ir] = d[:, 3]
    return zs, rs, ex, ez


def main():
    run, final, lam, s_end, step, a0, a1, tolx, tolz = sys.argv[1:10]
    continuum = len(sys.argv) > 10 and sys.argv[10] == "continuum"
    lam, s_end, step, a0, a1, tolx, tolz = map(float, (lam, s_end, step, a0, a1, tolx, tolz))
    k0 = 2.0 * np.pi / lam  # code units: lengths times k0
    zs, rs, f, _ = load(run + "/map_00000.txt")
    zs1, rs1, ex1, ez1 = load(final)
    assert np.allclose(zs, zs1) and np.allclose(rs, rs1)
    h = (zs[1] - zs[0]) * k0  # slice spacing
    D = step * k0             # step in s
    r = rs * k0
    dr = r[1] - r[0]
    A0, A1 = a0 * k0, a1 * k0
    nm = 30
    jtm = zeros(lambda x: jn(1, x), nm)  # j_1n
    jte = zeros(j1p, nm)                 # j'_1n
    # spectrum in zeta (code units)
    F = np.fft.rfft(f, axis=0)
    kap = 2.0 * np.pi * np.fft.rfftfreq(len(zs), d=h)
    # fine radial quadrature for the projections and norms on [0, a]
    def proj(a):
        rq = (np.arange(4000) + 0.5) * a / 4000
        wq = rq * a / 4000
        fq = np.array([np.interp(rq, r, F[k].real) + 1j * np.interp(rq, r, F[k].imag) for k in range(len(kap))])
        fq[:, rq > r[-1]] = 0.0
        out = []
        for fam, js in (("tm", jtm), ("te", jte)):
            for j in js:
                k = j / a
                J0, J1, J1p = jn(0, k * rq), jn(1, k * rq), j1p(k * rq)
                norm = np.pi * np.sum((k * k * J1p * J1p + J1 * J1 / (rq * rq)) * wq)
                # <E, e_n> = pi int f k J0(k r) r dr (TM), -pi int f k J0(k r) r dr (TE)
                c = np.pi * (fq * (k * J0)[None, :] * wq[None, :]).sum(axis=1) * (1.0 if fam == "tm" else -1.0)
                out.append((fam, j, c / norm, norm))
        return out

    # the initial field is masked at a0: project on the modes of a0
    modes = proj(A0)
    # phases of the steps s^0 -> s^1 with the radius of s^1 (a linear in s)
    nstep = int(round(s_end / step))
    am = A0 + (A1 - A0) * np.arange(1, nstep + 1) / nstep
    theta = kap * h
    # kappa_e of the trapezoidal rule in zeta (E_z), kappa in the continuum
    ke = kap if continuum else 2.0 / h * np.tan(np.minimum(theta, 3.1) / 2.0)
    Ex = np.zeros((len(kap), len(r)), complex)
    Ez = np.zeros((len(kap), len(r)), complex)
    Ef = np.zeros((len(kap), len(r)), complex)  # reconstruction of the initial field (check)
    kk = kap[:, None]
    inside = r < A1
    for fam, j, c, norm in modes:
        k1 = j / A1
        kp = (j / am)[None, :]
        with np.errstate(divide="ignore", invalid="ignore"):
            if continuum:
                wd = np.nan_to_num(2.0 * kk * kp * kp / (4.0 * kk * kk - kp * kp)) * D
            else:  # omega Delta of one step (only modulo 2 pi matters)
                th = theta[:, None]
                wd = 2.0 * np.arctan(D * h * kp * kp * np.sin(th) /
                                     (16.0 * np.sin(th / 2) ** 2 - h * h * kp * kp * np.cos(th / 2) ** 2))
        ph = np.exp(-1j * np.nan_to_num(wd).sum(axis=1))
        # adiabatic: the normalized amplitude c sqrt(norm) is kept; the norm of the final mode
        rq = (np.arange(4000) + 0.5) * A1 / 4000
        J1q, J1pq = jn(1, k1 * rq), j1p(k1 * rq)
        norm1 = np.pi * np.sum((k1 * k1 * J1pq * J1pq + J1q * J1q / (rq * rq)) * rq * A1 / 4000)
        cf = c * np.sqrt(norm / norm1) * ph
        if fam == "tm":
            er = k1 * j1p(k1 * r)
            ezm = -k1 * k1 * jn(1, k1 * r)  # div E_perp of grad(J1 cos phi) = -k^2 J1 cos phi
        else:
            er = -jn(1, k1 * r) / r
            ezm = 0.0 * r
        er = er * inside
        Ex += cf[:, None] * er[None, :]
        with np.errstate(divide="ignore", invalid="ignore"):
            Ez += np.where(kk > 0, 1j * cf[:, None] * (ezm * inside)[None, :] / np.where(kk > 0, ke[:, None], 1.0), 0.0)
        if A0 == A1:
            Ef += c[:, None] * er[None, :]
    ex_ref = np.fft.irfft(Ex, n=len(zs), axis=0)
    ez_ref = np.fft.irfft(Ez, n=len(zs), axis=0)
    cmp = r < A1 - 2.0 * dr
    if A0 == A1:
        f_rec = np.fft.irfft(Ef, n=len(zs), axis=0)
        print("modal reconstruction of the initial field: rms error %.2e of max" %
              (np.sqrt(np.mean((f_rec - f)[:, cmp] ** 2)) / np.abs(f).max()))
    ok = True
    for name, code, ref, tol in (("E_x", ex1, ex_ref, tolx), ("E_z", ez1, ez_ref, tolz)):
        mx = np.abs(ref[:, cmp]).max()
        err = np.sqrt(np.mean((code - ref)[:, cmp] ** 2)) / mx
        emax = np.abs(code - ref)[:, cmp].max() / mx
        # the field without the wall would differ: the change of the field from its start, for scale
        print("%s: max |ref| %.4e, max |code| %.4e, rms diff / max %.4f, max diff / max %.4f (tol %.3g)" %
              (name, mx, np.abs(code[:, cmp]).max(), err, emax, tol))
        ok = ok and err < tol
    change = np.sqrt(np.mean((ex_ref - f)[:, cmp] ** 2)) / np.abs(f).max()
    print("E_x: rms change from the initial field / max: %.3f (the scale of the test)" % change)
    sys.exit(0 if ok else 1)


main()
