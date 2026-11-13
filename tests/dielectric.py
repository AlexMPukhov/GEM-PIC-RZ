"""Reference of the dielectric test (tests/dielectric.sh): modal theory of a dielectric-lined metal tube (dielectric
wakefield accelerator) and of Cherenkov radiation in an open dielectric.

Structure: vacuum channel r < a, dielectric eps in a < r < b, perfectly conducting wall at r = b. A bunch at v = c
excites the synchronous TM_0n modes (phase velocity c). With fields ~ exp(i k (z - c t)):
  vacuum   E_z = E0 (uniform: k_perp = 0 at v = c), E_r = B_phi = -i k E0 r/2,
  dielectric E_z = C F(r), F = J0(q r) Y0(q b) - Y0(q r) J0(q b), q = k sqrt(eps - 1),
             E_r = (i k/q^2) C F'(r), B_phi = (i k eps/q^2) C F'(r)          (B in units of E/c),
matched (E_z, B_phi continuous) at r = a:  eps F'(a)/(q^2 F(a)) = -a/2.
Per mode (SI, time averages): stored energy per length u = eps0/4 int (eps |E|^2 + |B|^2) dA, power
P = eps0 c/2 int Re(E_r B_phi*) dA, group velocity beta_g = P/(u c) (energy velocity), loss factor per length
kappa = E0^2 / (4 u (1 - beta_g)). The on-axis wake of a Gaussian bunch (charge Q, rms length sigma) behind it:
  E_z(zeta) = -q_b sum_n 2 kappa_n exp(-k_n^2 sigma^2/2) cos(k_n (zeta - zeta_b)),
filled up to zeta_b - zeta < s (1/beta_g - 1) after the length s of the structure (the energy left at the
entrance travels at beta_g).

usage:
  python3 dielectric.py modes <a um> <b um> <eps>                       print the modes
  python3 dielectric.py dwa <axis file> <a> <b> <eps> <Q pC> <sigma um> <zeta_b um> <s um> <lambda0 um> <tol k> <tol amp>
       fits the on-axis E_z of the run (axis file: zeta, Ex, Ez, ...) behind the bunch by the first modes and compares
       the wavenumbers and amplitudes; then the end of the filled wake of mode 1 against s (1/beta_g - 1)
  python3 dielectric.py slip <axis 0> <axis 1> <vacuum axis 0> <vacuum axis 1> <eps> <s um> <h um> <Delta um> <lambda0 um> <tol>
       the slippage of a laser pulse in a uniform medium against the dispersion of the scheme
  python3 dielectric.py cherenkov <map file> <eps> <zeta_b um> <r_min> <r_max> <tol>
       the field front of a v = c bunch in an open dielectric lies on zeta = zeta_b - r sqrt(eps - 1)
"""
import sys

import numpy as np

EPS0, C, QE, ME = 8.8541878128e-12, 299792458.0, 1.602176634e-19, 9.1093837015e-31
TH = np.linspace(0.0, np.pi, 2001)


def _trap(y, x):
    return float(np.sum(0.5 * (y[1:] + y[:-1]) * np.diff(x)))


def jn(n, x):
    """J_n(x) = (1/pi) int_0^pi cos(n t - x sin t) dt"""
    return _trap(np.cos(n * TH - x * np.sin(TH)), TH) / np.pi


def yn(n, x):
    """Y_n(x) = (1/pi) int_0^pi sin(x sin t - n t) dt - (1/pi) int_0^inf (e^{nt} + (-1)^n e^{-nt}) e^{-x sinh t} dt"""
    a = _trap(np.sin(x * np.sin(TH) - n * TH), TH) / np.pi
    tmax = np.arcsinh(60.0 / x) + 1.0
    t = np.linspace(0.0, tmax, 20001)
    b = _trap((np.exp(n * t) + (-1) ** n * np.exp(-n * t)) * np.exp(-x * np.sinh(t)), t) / np.pi
    return a - b


def F(q, r, b):
    return jn(0, q * r) * yn(0, q * b) - yn(0, q * r) * jn(0, q * b)


def dF(q, r, b):
    return -q * (jn(1, q * r) * yn(0, q * b) - yn(1, q * r) * jn(0, q * b))


def disp(k, a, b, eps):
    q = k * np.sqrt(eps - 1.0)
    return eps * dF(q, a, b) + 0.5 * a * q * q * F(q, a, b)


def modes(a, b, eps, n=3):
    """first n synchronous TM_0n modes: k, beta_g, kappa (V/C/m for E0 = 1 V/m) with lengths in m"""
    ks = []
    k, dk = 0.02 / (b - a), 0.02 / (b - a)
    f0 = disp(k, a, b, eps)
    while len(ks) < n:
        k1 = k + dk
        f1 = disp(k1, a, b, eps)
        if f0 * f1 < 0 and abs(F(k1 * np.sqrt(eps - 1), a, b)) > 1e-12:
            lo, hi, flo = k, k1, f0
            for _ in range(60):
                mid = 0.5 * (lo + hi)
                fm = disp(mid, a, b, eps)
                if flo * fm <= 0:
                    hi = mid
                else:
                    lo, flo = mid, fm
            ks.append(0.5 * (lo + hi))
        k, f0 = k1, f1
    out = []
    for k in ks:
        q = k * np.sqrt(eps - 1.0)
        Cc = 1.0 / F(q, a, b)  # E0 = 1
        rv = np.linspace(0.0, a, 400)
        Ezv, Erv = np.ones_like(rv), k * rv / 2  # |E_r| = |B_phi|
        rd = np.linspace(a, b, 1200)
        Ezd = Cc * np.array([F(q, r, b) for r in rd])
        dFd = Cc * np.array([dF(q, r, b) for r in rd])
        Erd, Bd = k / q**2 * dFd, k * eps / q**2 * dFd
        u = EPS0 / 4 * (_trap((Ezv**2 + 2 * Erv**2) * 2 * np.pi * rv, rv) +
                        _trap((eps * (Ezd**2 + Erd**2) + Bd**2) * 2 * np.pi * rd, rd))
        P = EPS0 * C / 2 * (_trap(Erv**2 * 2 * np.pi * rv, rv) + _trap(Erd * Bd * 2 * np.pi * rd, rd))
        bg = P / (u * C)
        out.append((k, bg, 1.0 / (4 * u * (1 - bg))))
    return out


def _j0_j1x(x2):
    """J0(x) and J1(x)/x as power series in x^2 (x^2 < 0 gives I0, I1/x)"""
    j0, t = 0.0, 1.0
    for k in range(40):
        j0 += t
        t *= -x2 / 4 / ((k + 1) ** 2)
    j1, t = 0.0, 0.5
    for k in range(40):
        j1 += t
        t *= -x2 / 4 / ((k + 1) * (k + 2))
    return j0, j1


def kperp2(kap, om, chi, model):
    """k_perp^2 of a mode exp(i kap zeta - i om s) (GT2: lab frequency kap, k_z = kap - om) in a medium chi:
    Maxwell, or the reduced GEM system (d_s L dropped): om = (k^2 - kap^2 chi) / (2 kap + kap chi/2 - k^2/(2 kap eps))"""
    if model == "exact":
        return (1 + chi) * kap**2 - (kap - om) ** 2
    return (om * (2 * kap + kap * chi / 2) + kap**2 * chi) / (1 + om / (2 * kap * (1 + chi)))


def group_velocity(k, a, b, eps, model):
    """beta_g = d Omega / d k_z = 1 / (1 - d om / d kap) at the synchronous point (kap = k, om = 0) of the tube"""
    def D(kap, om):
        j0, j1x = _j0_j1x(kperp2(kap, om, 0.0, model) * a * a)
        kd2 = kperp2(kap, om, eps - 1.0, model)
        q = np.sqrt(kd2)
        return -a * j1x / j0 - eps / kd2 * dF(q, a, b) / F(q, a, b)
    h = k * 1e-5
    dDk = (D(k + h, 0.0) - D(k - h, 0.0)) / (2 * h)
    dDo = (D(k, h) - D(k, -h)) / (2 * h)
    return 1.0 / (1.0 + dDk / dDo)


def hem_det(w, kz, a, b, eps, m=1):
    """6x6 matching matrix (real): unknowns A, B (vacuum E_z, B_z/i), C, D (E_z in dielectric: J_m, Y_m), E, F (B_z/i)"""
    k1 = np.sqrt(w * w - kz * kz)
    q2 = eps * w * w - kz * kz
    q = np.sqrt(q2)
    J = lambda x: jn(m, x)
    Jp = lambda x: 0.5 * (jn(m - 1, x) - jn(m + 1, x))
    Y = lambda x: yn(m, x)
    Yp = lambda x: 0.5 * (yn(m - 1, x) - yn(m + 1, x))
    # field rows at radius r for a column (function f, derivative fp, k^2, eps_region): E_z, Bt, E_phi, B_phi/i
    def ez_col(f, fp, kk2, e, r):  # E_z = f(r)
        return np.array([f(r), 0.0, (-kz * m * f(r) / r) / kk2, (e * w * fp(r)) / kk2])
    def bz_col(f, fp, kk2, e, r):  # B_z = i f(r)
        return np.array([0.0, f(r), (w * fp(r)) / kk2, (-kz * m * f(r) / r) / kk2])
    fv, fvp = (lambda r: J(k1 * r)), (lambda r: k1 * Jp(k1 * r))
    fj, fjp = (lambda r: J(q * r)), (lambda r: q * Jp(q * r))
    fy, fyp = (lambda r: Y(q * r)), (lambda r: q * Yp(q * r))
    Mx = np.zeros((6, 6))
    # r = a: continuity (vacuum - dielectric)
    Mx[0:4, 0] = ez_col(fv, fvp, k1 * k1, 1.0, a)
    Mx[0:4, 1] = bz_col(fv, fvp, k1 * k1, 1.0, a)
    Mx[0:4, 2] = -ez_col(fj, fjp, q2, eps, a)
    Mx[0:4, 3] = -ez_col(fy, fyp, q2, eps, a)
    Mx[0:4, 4] = -bz_col(fj, fjp, q2, eps, a)
    Mx[0:4, 5] = -bz_col(fy, fyp, q2, eps, a)
    # r = b: E_z = 0, E_phi = 0
    for j, col in ((2, ez_col(fj, fjp, q2, eps, b)), (3, ez_col(fy, fyp, q2, eps, b)), (4, bz_col(fj, fjp, q2, eps, b)), (5, bz_col(fy, fyp, q2, eps, b))):
        Mx[4, j] = col[0]
        Mx[5, j] = col[2]
    # column scaling for conditioning
    Mx /= np.max(np.abs(Mx), axis=0)
    return np.linalg.det(Mx)

def hem_modes(a, b, eps, wmax, eta=1e-5):
    """synchronous HEM_1n wavenumbers below wmax (1/m): roots of the matching determinant on k_z = omega (1 - eta)
    (the vacuum fields of v = c are the limit k_perp -> 0; eta = 1e-5 is converged to 1e-5 relative)"""
    ws = np.linspace(0.05 / (b - a), wmax, 400)
    vals = [hem_det(w, w * (1 - eta), a, b, eps) for w in ws]
    roots = []
    for i in range(len(ws) - 1):
        if vals[i] * vals[i + 1] < 0:
            lo, hi, flo = ws[i], ws[i + 1], vals[i]
            for _ in range(50):
                mid = 0.5 * (lo + hi)
                fm = hem_det(mid, mid * (1 - eta), a, b, eps)
                if flo * fm <= 0:
                    hi = mid
                else:
                    lo, flo = mid, fm
            roots.append(0.5 * (lo + hi))
    return roots


def fit_k(z, f, kv, iters=40):
    """least squares of f(z) by sum_j (c_j cos k_j z + s_j sin k_j z) with free k_j (Gauss-Newton)"""
    def fit(k):
        A = np.column_stack([fn(kk * z) for kk in k for fn in (np.cos, np.sin)])
        c, *_ = np.linalg.lstsq(A, f, rcond=None)
        return c, A @ c - f
    kv = np.array(kv, dtype=float)
    for _ in range(iters):
        c, r0 = fit(kv)
        Jm = []
        for j in range(len(kv)):
            k2 = kv.copy(); k2[j] *= 1 + 1e-7
            Jm.append((fit(k2)[1] - r0) / (kv[j] * 1e-7))
        dk = np.linalg.lstsq(np.array(Jm).T, -r0, rcond=None)[0]
        dk = np.clip(dk, -0.01 * kv, 0.01 * kv)
        kv = kv + dk
        if np.max(np.abs(dk / kv)) < 1e-8:
            break
    c, r0 = fit(kv)
    return kv, c, r0


def dipole(z, Ex, a, b, eps, zb, L, tol):
    """on-axis E_x (the m = 1 field) behind an offset bunch: wavenumbers of HEM_11, HEM_12 (5 modes fitted)"""
    th = hem_modes(a * 1e-6, b * 1e-6, eps, 1.75e5)[:5]
    sel = (z < zb - 40) & (z > zb - L)
    kv, c, r0 = fit_k((z[sel] - zb) * 1e-6, Ex[sel], th)
    bad = 0
    for j in range(len(th)):
        e = abs(kv[j] / th[j] - 1)
        ok = e < tol if j < 2 else True
        bad += not ok
        print("  HEM1%d  k = %.4f 1/mm (theory %.4f)  dk %.2e  amplitude %.3e  %s" %
              (j + 1, kv[j] * 1e-3, th[j] * 1e-3, e, abs(c[2 * j] + 1j * c[2 * j + 1]), "ok" if ok else "FAIL")
              if j < 2 else "  HEM1%d  k = %.4f 1/mm (theory %.4f)  dk %.2e  (not checked)" %
              (j + 1, kv[j] * 1e-3, th[j] * 1e-3, e))
    print("  rms residual of the 5-mode fit: %.2e of the maximum" % (np.sqrt(np.mean(r0**2)) / np.max(np.abs(Ex[sel]))))
    return bad


def eunit(lam0_um):
    """m c omega0 / e in V/m"""
    return ME * C * 2 * np.pi * C / (lam0_um * 1e-6) / QE


def dwa(args):
    path = args[0]
    a, b, eps, Q, sig, zb, s, lam0, tolk, tola = map(float, args[1:11])
    md = modes(a * 1e-6, b * 1e-6, eps, 5)
    d = np.loadtxt(path)
    z, Ez = d[:, 0], d[:, 2]
    E1 = eunit(lam0)
    # model: electron bunch (charge -Q): E_z = +Q sum 2 kappa exp(-k^2 sig^2/2) cos(k (zeta - zeta_b)) in V/m
    print("  mode   k [1/mm] (theory)   beta_g (Maxwell)   kappa [V/pC/m]   amplitude [m c w0/e]: theory   code")
    ref = np.zeros_like(z)
    fill = []
    bgr = []
    for (k, bg, kap) in md:
        amp = Q * 1e-12 * 2 * kap * np.exp(-0.5 * (k * sig * 1e-6) ** 2) / E1
        bgr.append(group_velocity(k, a * 1e-6, b * 1e-6, eps, "reduced"))
        fill.append(s * (1 / bgr[-1] - 1))  # the code solves the reduced system: its group velocity
        ref += amp * np.cos(k * 1e-6 * (z - zb))
    # fit region: behind the bunch (4 sigma) and inside the filled part of the first three modes
    L = min(fill[:3])
    sel = (z < zb - 4 * sig) & (z > zb - 0.8 * L) & (z > z.min() + 5)
    zz, ee = z[sel] - zb, Ez[sel]
    # least squares with the theory wavenumbers, then a free nonlinear refinement of k for the two main modes
    def fit(kvec):
        A = np.column_stack([f(kk * 1e-6 * zz) for kk in kvec for f in (np.cos, np.sin)])
        c, *_ = np.linalg.lstsq(A, ee, rcond=None)
        return c, float(np.sum((A @ c - ee) ** 2))
    kv = np.array([m[0] for m in md])
    for it in range(30):  # Gauss-Newton on the wavenumbers (the amplitudes are linear: variable projection)
        c0, r0 = fit(kv)
        A = np.column_stack([f(kk * 1e-6 * zz) for kk in kv for f in (np.cos, np.sin)])
        res0 = A @ c0 - ee
        Jm = []
        for j in range(len(kv)):
            k2 = kv.copy(); k2[j] *= 1 + 1e-6
            c2, _ = fit(k2)
            A2 = np.column_stack([f(kk * 1e-6 * zz) for kk in k2 for f in (np.cos, np.sin)])
            Jm.append((A2 @ c2 - ee - res0) / (kv[j] * 1e-6))
        Jm = np.array(Jm).T
        dk, *_ = np.linalg.lstsq(Jm, -res0, rcond=None)
        dk = np.clip(dk, -0.01 * kv, 0.01 * kv)
        kv = kv + dk
        if np.max(np.abs(dk / kv)) < 1e-7:
            break
    c, _ = fit(kv)
    bad = 0
    for j, (k, bg, kap) in enumerate(md):
        amp_t = Q * 1e-12 * 2 * kap * np.exp(-0.5 * (k * sig * 1e-6) ** 2) / E1
        amp_c = c[2 * j]  # cos part (phase 0 at the bunch centre)
        sn = c[2 * j + 1]
        ek, ea = abs(kv[j] / k - 1), abs(amp_c / amp_t - 1)
        if j >= 3:
            continue
        tk, ta = (tolk, tola) if j < 2 else (3 * tolk, 3 * tola)
        ok = ek < tk and ea < ta
        bad += not ok
        print("  TM0%d   %.4f (%.4f)   %.4f   %.4g   %.4e   %.4e (sin part %.1e)   dk %.2e  damp %.2e  %s" %
              (j + 1, kv[j] * 1e-3, k * 1e-3, bg, kap * 1e-12, amp_t, amp_c, sn, ek, ea, "ok" if ok else "FAIL"))
    rms = np.sqrt(np.mean((ee - np.interp(zz + zb, z, ref)) ** 2)) / np.max(np.abs(ref[sel]))
    print("  rms difference to the modal sum (3 modes) behind the bunch: %.2e of the maximum" % rms)
    # end of the filled wake of mode 1: envelope of the field filtered at k_1
    k1, bg1 = md[0][0] * 1e-6, md[0][1]
    zr = z - zb
    win = 2 * np.pi / k1
    nw = max(1, int(round(win / (z[1] - z[0]))))
    ker = np.ones(nw) / nw
    env = np.abs(np.convolve(Ez * np.exp(-1j * k1 * zr), ker, mode="same")) * 2
    amp1 = abs(c[0] + 1j * c[1])
    behind = np.where((zr < -4 * sig) & (env < 0.5 * amp1))[0]
    behind = behind[zr[behind] < -0.3 * fill[0]] if len(behind) else behind
    if len(behind):
        zt = zr[behind].max()
        et = abs(-zt / fill[0] - 1)
        ok = et < 0.03
        bad += not ok
        print("  end of the filled TM01 wake: %.1f um behind the bunch, reduced model s (1/beta_g - 1) = %.1f um "
              "(beta_g = %.4f; %.1f %%) %s" % (-zt, fill[0], bgr[0], 100 * et, "ok" if ok else "FAIL"))
        print("  (Maxwell: beta_g = %.4f, %.1f um: the GEM system drops d_s L, exact only for the synchronous field)"
              % (md[0][1], s * (1 / md[0][1] - 1)))
    else:
        print("  end of the filled wake not inside the window (theory %.1f um behind the bunch)" % fill[0])
    if np.max(np.abs(d[:, 1])) > 0:  # an offset bunch: the m = 1 field on the axis
        bad += dipole(z, d[:, 1], a, b, eps, zb, 600.0, 2 * tolk)
    return bad


def cherenkov(args):
    """peak of |E_z| at each radius r0 < r < r1 (x > 0 of the map), searched in zeta_b - 1.5 r tan - 10 um ... zeta_b + 10 um
    (ahead of the cone reflected at the wall); its position against zeta_b - r tan(theta), tan(theta) = sqrt(eps - 1)"""
    path = args[0]
    eps, zb, r0, r1, tol = map(float, args[1:6])
    d = np.loadtxt(path)
    z, x, Ez = d[:, 0], d[:, 1], d[:, 3]
    tanth = np.sqrt(eps - 1)
    rs = np.unique(x[(x > r0) & (x < r1)])
    pk = []
    for r in rs:
        sel = (x == r) & (z > zb - 1.5 * r * tanth - 10) & (z < zb + 10)
        zz, ee = z[sel], np.abs(Ez[sel])
        i = int(np.argmax(ee))
        if 0 < i < len(zz) - 1:  # parabola through the peak
            y0, y1, y2 = ee[i - 1], ee[i], ee[i + 1]
            pk.append(zz[i] + 0.5 * (y0 - y2) / (y0 - 2 * y1 + y2) * (zz[1] - zz[0]))
        else:
            pk.append(zz[i])
    p = np.polyfit(rs, np.array(pk), 1)
    err = abs(-p[0] / tanth - 1)
    print("  Cherenkov cone: |E_z| peak at zeta = %.2f - %.4f r um, theory zeta_b - r sqrt(eps - 1) = %.2f - %.4f r "
          "(slope %.2f %%) %s" % (p[1], -p[0], zb, tanth, 100 * err, "ok" if err < tol else "FAIL"))
    return 0 if err < tol else 1


def pattern_shift(f0, f1):
    """b(zeta) = a(zeta - d) of the on-axis E_x of two axis files: d from the phase slope of the cross spectrum"""
    d0, d1 = np.loadtxt(f0), np.loadtxt(f1)
    h = d0[1, 0] - d0[0, 0]
    A, B = np.fft.rfft(d0[:, 1]), np.fft.rfft(d1[:, 1])
    k = np.fft.rfftfreq(len(d0), h) * 2 * np.pi
    X = np.conj(A) * B
    w = np.abs(X)
    sel = w > 0.1 * w.max()
    return -np.polyfit(k[sel], np.unwrap(np.angle(X[sel])), 1, w=w[sel])[0]


def slip(args):
    """laser in a uniform medium: the pulse slips back by (n_g - 1) s. Reference: the scheme in the limit k_perp -> 0,
    L = -chi R/(4 + chi), J from the trapezoidal rule (d_zeta -> (2i/h) tan(kappa h/2)) and the trapezoid in s:
    omega(kappa) = (2/Delta) arctan((Delta u/h) tan(kappa h/2)), u = 2 chi/(4 + chi), n_g - 1 = d omega/d kappa;
    continuum of the GEM system 2 chi/(4 + chi), Maxwell sqrt(1 + chi) - 1"""
    d0, d1, v0, v1 = args[:4]
    eps, s, h, D, lam0, tol = map(float, args[4:10])
    chi = eps - 1
    meas = -(pattern_shift(d0, d1) - pattern_shift(v0, v1)) / s
    u, k0 = 2 * chi / (4 + chi), 2 * np.pi / lam0
    om = lambda kk: (2 / D) * np.arctan((D * u / h) * np.tan(kk * h / 2))
    sch = (om(k0 * (1 + 1e-5)) - om(k0 * (1 - 1e-5))) / (2e-5 * k0)
    err = abs(meas / sch - 1)
    print("  laser in eps = %g over %.0f um: n_g - 1 = %.6f measured, %.6f scheme (%.1e), %.6f GEM continuum, %.6f Maxwell %s"
          % (eps, s, meas, sch, err, u, np.sqrt(1 + chi) - 1, "ok" if err < tol else "FAIL"))
    return 0 if err < tol else 1


if __name__ == "__main__":
    what = sys.argv[1]
    if what == "modes":
        a, b, eps = map(float, sys.argv[2:5])
        for k, bg, kap in modes(a * 1e-6, b * 1e-6, eps):
            print("k = %.5f 1/mm  lambda = %.2f um  beta_g = %.4f  kappa = %.4g V/pC/m" %
                  (k * 1e-3, 2 * np.pi / k * 1e6, bg, kap * 1e-12))
    elif what == "dwa":
        sys.exit(1 if dwa(sys.argv[2:]) else 0)
    elif what == "slip":
        sys.exit(slip(sys.argv[2:]))
    elif what == "cherenkov":
        sys.exit(1 if cherenkov(sys.argv[2:]) else 0)
