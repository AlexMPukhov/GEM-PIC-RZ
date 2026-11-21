"""Quantum diffusion of the energy in an undulator from the classical radiation integral (all harmonics), the
reference of the undulator factor of numerics.radiation_reaction (src/radiation.hpp, radUndulatorFactor).

R(K) = [d<dgamma^2>/dz] / (r_e lambdabar_C gamma^4 k_u^3), k_u = 1, c = 1, gamma >> 1. Jackson 14.67 for one period:
  energy at harmonic n per solid angle  dI_n/dOmega = (e^2 w_n^2 / 4 pi^2) |M_n|^2 (2 pi / T_eff),
  M_n = int_period n x (n x beta) exp(i w_n (t - n.r)) dt = FFT over the retarded time tau (uniform grid, dt/dtau),
  second moment per period = sum_n int dOmega hbar w_n dI_n/dOmega  ->  R = sum int w_n^3 |M_n|^2 / (2 pi T_eff) dOmega / (2 pi gamma^4).
Limits (derived in the average rest frame and from the local constant field): helical (14/15) K^2 for K -> 0 and
1.4174 K^3 for K -> inf; planar (peak K) (7/15) K^2 and 1.2032 (7/15) K^3 = 0.5615 K^3.
Results (gamma = 400): F(K) = R / (c K^2), c = 14/15 helical, 7/15 planar:
  helical: Saldin, Schneidmiller, Yurkov (NIM A 381, 545 (1996)) F = 1.42 K + 1 / (1 + 1.5 K + 0.95 K^2) to 0.5 % (K 0.25-3)
  planar:  F = 1.2032 K + 1 / (1 + 1.34 K + 0.43 K^2) to 0.8 % (K 0.25-4)
usage: python3 undulator_diffusion.py planar|helical K [K ...]   (~1 min per K)"""
import sys

import numpy as np


def R(K, kind, g=400.0, Nt=1024, nth=160, nph=48):
    Ns = 8192
    s = np.arange(Ns) * 2 * np.pi / Ns                    # phase k_u z over one period
    bx = K / g * np.cos(s)
    by = 0 * s if kind == "planar" else K / g * np.sin(s)
    bz = np.sqrt(1 - 1 / g**2 - bx**2 - by**2)
    dt = (2 * np.pi / Ns) / bz
    t = np.concatenate([[0], np.cumsum(dt)[:-1]])
    T = dt.sum()
    x = np.concatenate([[0], np.cumsum(bx * dt)[:-1]])
    y = np.concatenate([[0], np.cumsum(by * dt)[:-1]])
    x -= x.mean(); y -= y.mean()
    Kf = 1 + K * K / 2 if kind == "planar" else 1 + K * K
    thmax = 4 * np.sqrt(Kf) / g + (2 * K / g if kind == "planar" else 0)
    th = (np.arange(nth) + 0.5) * thmax / nth
    ph = (np.arange(nph) + 0.5) * 2 * np.pi / nph if kind == "planar" else np.array([0.0])
    b = np.array([bx, by, bz])
    tot = 0.0
    for a in th:
        for p in ph:
            n = np.array([np.sin(a) * np.cos(p), np.sin(a) * np.sin(p), np.cos(a)])
            tau = t - (n[0] * x + n[1] * y + n[2] * s)
            Teff = T - n[2] * 2 * np.pi
            nb = n @ b
            v = (n[:, None] * nb - b) / (1 - nb)              # n x (n x beta) dt/dtau
            tu = tau[0] + np.arange(Nt) * Teff / Nt
            f = np.array([np.interp(tu, tau, v[c]) for c in range(3)])
            M = np.fft.fft(f, axis=1) * (Teff / Nt)
            nn = np.arange(1, Nt // 2)
            wn = 2 * np.pi * nn / Teff
            dOm = (thmax / nth) * np.sin(a) * (2 * np.pi / len(ph))
            tot += np.sum(wn**3 * np.sum(np.abs(M[:, nn]) ** 2, axis=0)) / (2 * np.pi * Teff) * dOm
    return tot / (2 * np.pi * g**4)


if __name__ == "__main__":
    kind = sys.argv[1]
    c, A, a, b = (14 / 15, 1.42, 1.5, 0.95) if kind == "helical" else (7 / 15, 1.2032, 1.34, 0.43)
    for K in map(float, sys.argv[2:]):
        r = R(K, kind)
        fit = c * K * K * (A * K + 1 / (1 + a * K + b * K * K))
        print(f"{kind} K = {K:g}: R = {r:.5g}, fit {fit:.5g} ({r / fit - 1:+.4f}), factor to the local field "
              f"{r / (c * K * K * (1.4174 if kind == 'helical' else 1.2032) * K):.4f}")
