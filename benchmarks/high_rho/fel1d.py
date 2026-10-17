"""1D steady-state FEL amplifier in a planar undulator, two models (independent reference for the high-rho study):
  averaged   the standard period-averaged equations (coupling K [JJ]/2, mean longitudinal velocity), as in GENESIS
  unaveraged the electron motion in the undulator kept explicitly: v_x/c = K cos(k_u z)/gamma, the figure-eight
             oscillation of v_z, coupling to the field at every z (the field itself is a slowly varying amplitude of
             one frequency, d^2E/dz^2 neglected: error ~ rho lambda_r / lambda_u)
Electrons: one ponderomotive bucket, N macroparticles with uniform phases (quiet start), cold beam. SI units; the beam
is uniform over the area A (1D). Optional longitudinal space charge at lambda_r (the 1D field of the bunching).
usage (module): run(gamma0, current, area, K, lambda_u, z_end, P0, model, ...) -> z, P, bunching"""
import numpy as np
from scipy.special import j0, j1

e, me, c, eps0 = 1.602176634e-19, 9.1093837015e-31, 2.99792458e8, 8.8541878128e-12
Z0 = 1.0 / (eps0 * c)
IA = 4 * np.pi * eps0 * me * c**3 / e   # Alfven current 17.045 kA


def JJ(K):
    xi = K**2 / (4 + 2 * K**2)
    return j0(xi) - j1(xi)


def rho_1d(gamma0, current, area, K, lambda_u):
    """Pierce parameter of a uniform beam of area A (= the Xie value for a Gaussian of 2 pi sigma^2 = A)"""
    ku = 2 * np.pi / lambda_u
    return (K**2 * JJ(K)**2 * current * 2 * np.pi / (16 * IA * gamma0**3 * area * ku**2)) ** (1 / 3)


def run(gamma0, current, area, K, lambda_u, z_end, P0, model="unaveraged", npart=256, steps_per_period=96,
        space_charge=False, detune=0.0):
    ku = 2 * np.pi / lambda_u
    lam = lambda_u * (1 + K**2 / 2) / (2 * gamma0**2) * (1 + detune)
    k = 2 * np.pi / lam
    JJK = K * JJ(K)
    J0 = current / area                       # |current density|
    E0 = np.sqrt(2 * Z0 * P0 / area)          # P = |E|^2 A / (2 Z0)
    g = np.full(npart, float(gamma0))
    psi = 2 * np.pi * (np.arange(npart) + 0.5) / npart
    E = complex(E0)
    q = e / (me * c**2)

    def deriv(z, g, psi, E):
        if model == "averaged":
            th = psi
            bz = np.sqrt(1 - (1 + K**2 / 2) / g**2)
            coup = JJK / 2
            dg = -q * coup / g * np.real(E * np.exp(1j * th))
            dpsi = ku + k - k / bz
            src = np.mean(np.exp(-1j * th) / g)
            dE = Z0 * J0 * coup * src
        else:
            cz = np.cos(ku * z)
            bz = np.sqrt(1 - (1 + (K * cz)**2) / g**2)
            vxz = K * cz / (g * bz)           # v_x / v_z
            dg = -q * vxz * np.real(E * np.exp(1j * psi))
            dpsi = k - k / bz
            dE = Z0 * J0 * np.mean(vxz * np.exp(-1j * psi))
            th = psi + ku * z
        if space_charge:  # 1D longitudinal field of the bunching at k (ponderomotive phase), E_z = Re[Ez e^{i th}]
            bunch = np.mean(np.exp(-1j * th))
            Ez = 1j * J0 / (eps0 * c * k) * bunch * 2
            dg = dg - q * np.real(Ez * np.exp(1j * th))
        return dg, dpsi, dE

    dz = lambda_u / steps_per_period
    n = int(round(z_end / dz))
    zs = np.empty(n + 1)
    P = np.empty(n + 1)
    b = np.empty(n + 1)
    z = 0.0
    for i in range(n + 1):
        th = psi + (0 if model == "averaged" else ku * z)
        zs[i], P[i], b[i] = z, abs(E)**2 * area / (2 * Z0), abs(np.mean(np.exp(1j * th)))
        if i == n:
            break
        k1 = deriv(z, g, psi, E)
        k2 = deriv(z + dz / 2, g + dz / 2 * k1[0], psi + dz / 2 * k1[1], E + dz / 2 * k1[2])
        k3 = deriv(z + dz / 2, g + dz / 2 * k2[0], psi + dz / 2 * k2[1], E + dz / 2 * k2[2])
        k4 = deriv(z + dz, g + dz * k3[0], psi + dz * k3[1], E + dz * k3[2])
        g = g + dz / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
        psi = psi + dz / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
        E = E + dz / 6 * (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2])
        z += dz
    return zs, P, b


def gain_length(z, P, lo=1e2, hi=1e-2):
    """power gain length from the exponential part: P between lo P(0) and hi P_max"""
    m = (P > lo * P[0]) & (P < hi * P.max()) & (z < z[np.argmax(P)])
    if m.sum() < 5:
        return np.nan
    return 1 / np.polyfit(z[m], np.log(P[m]), 1)[0]
