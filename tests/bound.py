"""Reference of the bound-electron test (tests/bound.sh): Maxwell-Bloch medium (block bound_electrons).

linear: weak field. The on-axis E_x(zeta) of a laser after the length z of a uniform medium, compared per Fourier
  component kappa with the initial one: the ratio must be exp(+-i (n(kappa) - 1) kappa z) with the Lorentz index
  n = 1 + 2 chi / (4 + chi) of the GEM system (= sqrt(1 + chi) to O(chi^2)), chi = n_a sum_j f_j / (w_j^2 - kappa^2)
  in code units (n_a in n_c, frequencies in omega0), and |ratio| = 1 (no damping); both divided by the same ratio of
  a run in vacuum (diffraction). The discrete susceptibility of
  the Bloch step, chi_h = n_a sum f_j (w_j'/w_j) / (w_j'^2 - kappa'^2), x' = (2/h) tan(x h/2), is printed as well.
  usage: bound.py linear <axis 0> <axis z> <vacuum axis 0> <vacuum axis z> <z um> <lambda0 um> <n_a cm^-3> <h um> <tol> <lambda_1 um> <f_1> [<lambda_2> <f_2>]

rabi: Rabi flopping of a thin resonant medium. The pulse area is measured from the code's own field at the atoms of
  the first cell centre (map file; or the axis field),
  A = a |2 int E_x(zeta) exp(i kappa zeta) d zeta| (kappa = omega0: the envelope integral of E = env cos), a = d/(e a0
  alpha) with d^2 = f E_h / (2 hbar omega). Reference: the two-level Schroedinger equation without the rotating-wave
  approximation integrated (RK4, field refined 8x) in the code's on-axis field; sin^2(A/2) is printed for comparison.
  usage: bound.py rabi <map or axis file> <run log> <lambda0 um> <lambda_j um> <f> <tol>
"""
import sys

import numpy as np


def nc_cm3(lam0):
    return 1.1148e21 / lam0**2


def chi(kap, na, lines, h=None):
    c = np.zeros_like(kap)
    for w, f in lines:
        if h is None:
            c += na * f / (w * w - kap * kap)
        else:
            T = lambda x: (2 / h) * np.tan(0.5 * x * h)
            c += na * f * (T(w) / w) / (T(w) ** 2 - T(kap) ** 2)
    return c


def linear(a):
    f0, f1, v0, v1 = a[0], a[1], a[2], a[3]
    z, lam0, na_cm3, h_um, tol = map(float, a[4:9])
    rest = list(map(float, a[9:]))
    lines = [(lam0 / rest[i], rest[i + 1]) for i in range(0, len(rest), 2)]
    k0 = 2 * np.pi / lam0
    na = na_cm3 / nc_cm3(lam0)
    d0, d1 = np.loadtxt(f0), np.loadtxt(f1)
    hz = d0[1, 0] - d0[0, 0]
    A, B = np.fft.rfft(d0[:, 1]), np.fft.rfft(d1[:, 1])
    kap = np.fft.rfftfreq(len(d0), hz) * 2 * np.pi / k0  # in k0
    sel = np.abs(A) > 0.1 * np.abs(A).max()
    # the same run in vacuum: diffraction (Gouy phase, the amplitude on the axis) is removed
    w0_, w1_ = np.loadtxt(v0), np.loadtxt(v1)
    Av, Bv = np.fft.rfft(w0_[:, 1]), np.fft.rfft(w1_[:, 1])
    ratio = (B[sel] / A[sel]) / (Bv[sel] / Av[sel])
    ph = np.unwrap(np.angle(ratio))
    amp = np.abs(ratio)
    zc = z * k0
    c = chi(kap[sel], na, lines)
    ref = 2 * c / (4 + c) * kap[sel] * zc
    ch = chi(kap[sel], na, lines, h_um * k0)
    refh = 2 * ch / (4 + ch) * kap[sel] * zc
    # sign convention of the transform: compare with the sign that fits; the phase is defined up to 2 pi
    best = None
    for sg in (1.0, -1.0):
        d = sg * ph - ref
        d -= 2 * np.pi * np.round(np.mean(d) / (2 * np.pi))
        e = np.sqrt(np.mean(d**2)) / np.max(np.abs(ref))
        if best is None or e < best[0]:
            best = (e, sg, d)
    e, sg, d = best
    dh = sg * ph - refh
    dh -= 2 * np.pi * np.round(np.mean(dh) / (2 * np.pi))
    eh = np.sqrt(np.mean(dh**2)) / np.max(np.abs(ref))
    core = np.abs(A[sel]) > 0.3 * np.abs(A).max()  # amplitude: inside the band (the edges carry little energy)
    ea = np.max(np.abs(amp[core] - 1))
    print("  band kappa = %.4f ... %.4f k0, chi = %.3e ... %.3e; phase (n - 1) kappa z = %.4f ... %.4f rad"
          % (kap[sel].min(), kap[sel].max(), c.min(), c.max(), ref.min(), ref.max()))
    print("  rms phase difference / max phase: Lorentz %.2e, discrete chi_h %.2e; max | |ratio| - 1 | = %.1e"
          % (e, eh, ea))
    ok = e < tol and ea < tol
    print("  %s" % ("ok" if ok else "FAIL"))
    return 0 if ok else 1


def rabi(a):
    ax, log = a[0], a[1]
    lam0, lamj, f, tol = map(float, a[2:6])
    d = np.loadtxt(ax)
    if d.shape[1] >= 4 and np.any(d[:, 1] != d[0, 1]):  # a map (zeta, x, Ex, Ez): the field at the first centre x > 0
        x0 = np.min(d[d[:, 1] > 0, 1])
        d = d[d[:, 1] == x0][:, [0, 2]]
        d = d[np.argsort(d[:, 0])]
    k0 = 2 * np.pi / lam0
    z = d[:, 0] * k0  # code units
    A_env = abs(2 * np.trapezoid(d[:, 1] * np.exp(1j * z), z)) if hasattr(np, "trapezoid") else         abs(2 * np.trapz(d[:, 1] * np.exp(1j * z), z))
    dip = np.sqrt(f * 27.211386 / (2 * 1.2398420 / lamj))
    A = dip * 137.035999084 * A_env
    p = None
    for line in open(log):
        if "excited on the axis at the rear" in line:
            p = float(line.split("excited on the axis at the rear")[1].split()[0])
    rwa = np.sin(A / 2) ** 2
    # exact two-level dynamics (no rotating-wave approximation) in the code's on-axis field: H = [[0, -a E], [-a E, w]],
    # time t = -zeta from the front, RK4 on the field refined 8 times by band-limited (Fourier) interpolation
    a_ = dip * 137.035999084
    w = lam0 / lamj
    E = d[::-1, 1]  # time order: front first
    n = len(E)
    F = np.fft.rfft(E)
    ref8 = 8
    Ef = np.fft.irfft(F, n * ref8) * ref8
    dt = (z[1] - z[0]) / ref8
    c = np.array([1.0 + 0j, 0.0 + 0j])
    def rhs(c, e):
        return -1j * np.array([-a_ * e * c[1], -a_ * e * c[0] + w * c[1]])
    for k in range(len(Ef) - 1):
        e0, e1 = Ef[k], Ef[k + 1]
        em = 0.5 * (e0 + e1)
        k1 = rhs(c, e0); k2 = rhs(c + 0.5 * dt * k1, em); k3 = rhs(c + 0.5 * dt * k2, em); k4 = rhs(c + dt * k3, e1)
        c = c + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    ref = abs(c[1]) ** 2 / (abs(c[0]) ** 2 + abs(c[1]) ** 2)
    ok = p is not None and abs(p - ref) < tol
    print("  pulse area on the axis A = %.5f pi: excited population %.6f, exact two-level %.6f (diff %.1e), "
          "rotating-wave sin^2(A/2) %.6f %s" % (A / np.pi, p, ref, abs(p - ref), rwa, "ok" if ok else "FAIL"))
    return 0 if ok else 1


def static(a):
    """static limit: a line far above the laser (w_j = 10) against the dielectric with eps = 1 + chi(omega0): the
    group delays differ by the dispersion of the Lorentz line, (n_g - 1)_bound - (n_g - 1)_eps = n_a f w / (w_j^2 - w^2)^2
    (GEM: n - 1 = chi / 2 to O(chi^2)); the group delay from the pattern shift against a vacuum run"""
    sys.path.insert(0, __file__.rsplit("/", 1)[0])
    import dielectric as D
    b0, b1, d0, d1, v0, v1 = a[:6]
    z, lam0, lamj, f, na_cm3, tol = map(float, a[6:12])
    vac = D.pattern_shift(v0, v1)
    gb = -(D.pattern_shift(b0, b1) - vac) / z
    gd = -(D.pattern_shift(d0, d1) - vac) / z
    na = na_cm3 / nc_cm3(lam0)
    w = lam0 / lamj
    disp = na * f / (w * w - 1) ** 2  # (omega/2) d chi / d omega at omega = 1
    err = abs((gb - gd) / disp - 1)
    print("  static limit: n_g - 1 = %.6f (bound, w_j = %g), %.6f (eps = 1 + chi = %.6f); difference %.3e, Lorentz "
          "dispersion %.3e (%.1f %%) %s" % (gb, w, gd, 1 + na * f / (w * w - 1), gb - gd, disp, 100 * err,
                                            "ok" if err < tol else "FAIL"))
    return 0 if err < tol else 1


def adk_code(U_eV, E_code, lam0):
    """ADK rate (ionization.hpp: s-state, Z after = 1, capped at its maximum) [1/omega0] in the field E [m c omega0/e]"""
    from math import gamma
    UH, Ea, wa = 13.605693122994, 5.14220674763e11, 4.1341373335e16
    w0 = 2 * np.pi * 299792458.0 / (lam0 * 1e-6)
    Eunit = 9.1093837015e-31 * 299792458.0 * w0 / 1.602176634e-19
    E = np.abs(np.asarray(E_code, dtype=float)) * Eunit
    ns = np.sqrt(UH / U_eV)
    C2 = 2 ** (2 * ns) / (ns * gamma(ns + 1) * gamma(ns))
    x = (U_eV / UH) ** 1.5 * Ea
    Epk = (2.0 / 3.0) * x / (2 * ns - 1)
    Ee = np.minimum(np.maximum(E, 1e-30), Epk)
    w = wa * C2 * (U_eV / (2 * UH)) * (2 * x / Ee) ** (2 * ns - 1) * np.exp(-(2.0 / 3.0) * x / Ee)
    return np.where(E > 0, w / w0, 0.0)


def vsystem(E, dt, lines, ips, lam0):
    """V-type atom (ground, excited lines) with ADK losses per level in the field samples E (time order), RK4 with the
    field refined 8x (band-limited): returns the final density matrix"""
    n = len(E)
    F = np.fft.rfft(E)
    Ef = np.fft.irfft(F, n * 8) * 8
    h = dt / 8
    nl = len(lines)
    w = [lam0 / l for l, f in lines]
    a = [np.sqrt(f * 27.211386 / (2 * 1.2398420 / l)) * 137.035999084 for l, f in lines]
    def rhs(rho, e):
        H = np.zeros((nl + 1, nl + 1))
        for j in range(nl):
            H[0, j + 1] = H[j + 1, 0] = -a[j] * e
            H[j + 1, j + 1] = w[j]
        return -1j * (H @ rho - rho @ H)
    def loss(rho, e, t):  # exact decay of the levels over t (stiff rates beyond barrier suppression)
        g = np.array([adk_code(ips[k], e, lam0) for k in range(nl + 1)])
        f = np.exp(-0.5 * g * t)
        return rho * np.outer(f, f)
    rho = np.zeros((nl + 1, nl + 1), dtype=complex)
    rho[0, 0] = 1.0
    for k in range(len(Ef) - 1):
        e0, e1 = Ef[k], Ef[k + 1]
        em = 0.5 * (e0 + e1)
        rho = loss(rho, e0, 0.5 * h)
        k1 = rhs(rho, e0); k2 = rhs(rho + 0.5 * h * k1, em); k3 = rhs(rho + 0.5 * h * k2, em); k4 = rhs(rho + h * k3, e1)
        rho = rho + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        rho = loss(rho, e1, 0.5 * h)
    return rho


def ionize(a):
    """neutral fraction on the axis (grid) against the V-system with ADK per level in the code's field at the first
    centre; the particles' mean charge state leaving the window against the grid's area mean"""
    mp, log = a[0], a[1]
    lam0, ipg, tol, tolp = map(float, a[2:6])
    rest = list(map(float, a[6:]))
    lines = [(rest[i], rest[i + 1]) for i in range(0, len(rest), 2)]
    ips = [ipg] + [ipg - 1.2398420 / l for l, f in lines]
    d = np.loadtxt(mp)
    x0 = np.min(d[d[:, 1] > 0, 1])
    d = d[d[:, 1] == x0][:, [0, 2]]
    d = d[np.argsort(d[:, 0])]
    k0 = 2 * np.pi / lam0
    E = d[::-1, 1]
    rho = vsystem(E, (d[1, 0] - d[0, 0]) * k0, lines, ips, lam0)
    nu_ref = np.real(np.trace(rho))
    ex_ref = np.real(rho[1:, 1:].trace()) / nu_ref
    nu_ax = nu_mean = zbar = None
    for line in open(log):
        if "neutral on the axis" in line:
            nu_ax = float(line.split("neutral on the axis")[1].split(",")[0])
            nu_mean = float(line.split("area mean")[1].split()[0])
        if "<Z> leaving the window" in line:
            zbar = float(line.split("<Z> leaving the window =")[1].split()[0])
    e1 = abs(nu_ax - nu_ref)
    e2 = abs(zbar - (1 - nu_mean))
    ok = e1 < tol and e2 < tolp
    print("  neutral fraction on the axis: grid %.5f, V-system with ADK per level %.5f (diff %.1e); excited share of the"
          " neutrals %.4f" % (nu_ax, nu_ref, e1, ex_ref))
    print("  particles: <Z> leaving the window %.5f, grid area mean of the ionized fraction %.5f (diff %.1e) %s"
          % (zbar, 1 - nu_mean, e2, "ok" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    if sys.argv[1] == "linear":
        sys.exit(linear(sys.argv[2:]))
    if sys.argv[1] == "rabi":
        sys.exit(rabi(sys.argv[2:]))
    if sys.argv[1] == "static":
        sys.exit(static(sys.argv[2:]))
    if sys.argv[1] == "ionize":
        sys.exit(ionize(sys.argv[2:]))
