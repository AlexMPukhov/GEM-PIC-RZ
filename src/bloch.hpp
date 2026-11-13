// ============================================================================
//  bloch.hpp -- bound electrons: Maxwell-Bloch medium of a V-type atom
//  (block 'bound_electrons'; README "Bound electrons")
//
//  Every node holds the density matrix rho of a V-type atom: ground state 0,
//  excited states 1, 2 (e.g. Rb 5S, 5P_1/2, 5P_3/2), Hamiltonian (units of
//  omega0, field E = E_x in m c omega0/e)
//      H = [[0, -a1 E, -a2 E], [-a1 E, w1, 0], [-a2 E, 0, w2]],
//  a_j = d_j / (e a0 alpha) the Rabi frequency per unit field, w_j the
//  transition frequencies. The carrier is resolved (no rotating-wave
//  approximation). Polarization per atom density P = sum_j cP_j Re rho_0j,
//  cP_j = 2 d_j/(e a0) a0 k0; with n_a in n_c the medium current is
//  J = d_t (n_a P). For rho_00 = 1 the response is that of Lorentz
//  oscillators, chi = n_a sum_j f_j / (w_j^2 - omega^2).
//
//  One step over a half-cell (time h, field E0 at its start, E1 at its end):
//  Strang splitting of the losses (ionization rates g0, g1, g2 of the levels,
//  relaxation T1, T2: exact exponentials over h/2) around the unitary Cayley
//  step rho -> U rho U^+, U = (1 + i K)^-1 (1 - i K), K = (h/2) H((E0 + E1)/2)
//  (implicit midpoint of the Schroedinger equation: unconditionally stable,
//  second order, keeps the trace without losses), with the transition
//  frequencies pre-warped (exact free precession). Weak-field susceptibility
//  of the step: chi_h = n_a sum_j f_j (w_j'/w_j) / (w_j'^2 - Omega'^2),
//  x' = (2/h) tan(x h/2).
//
//  The node functions are GEM_HD (the serial code and the engine share them).
// ============================================================================
#pragma once
#include <math.h>

#include "portable.hpp"

namespace gem {

// constants of the lines (code units)
struct BlochLines {
  int n = 0;                     // 1 or 2 excited states
  double w[2] = {0, 0};          // transition frequencies
  double a[2] = {0, 0};          // Rabi frequency per unit field
  double cP[2] = {0, 0};         // polarization per atom density per Re rho_0j
  double T1 = 0.0, T2 = 0.0;     // relaxation (0: none)
};

// density matrix of one node: populations p0, p1, p2; coherences c01, c02, c12 (9 doubles)
struct BlochNode {
  double p0, p1, p2;
  cplx c01, c02, c12;
};

GEM_HD inline void blochGround(BlochNode& s) {
  s.p0 = 1.0; s.p1 = 0.0; s.p2 = 0.0;
  s.c01 = cplx(); s.c02 = cplx(); s.c12 = cplx();
}

GEM_HD inline double blochP(const BlochLines& L, const BlochNode& s) {
  return L.cP[0] * s.c01.re + (L.n > 1 ? L.cP[1] * s.c02.re : 0.0);
}

// d P / dt at field E (the slice value of the polarization current per atom density, without the losses of the
// ionization, which remove the dipole with the atom; see BoundMedium)
GEM_HD inline double blochJ(const BlochLines& L, const BlochNode& s, double E) {
  // (d rho/dt)_0j = -i [H, rho]_0j;  H_00 = 0, H_0j = H_j0 = -a_j E, H_jj = w_j, H_12 = 0
  const double h1 = -L.a[0] * E, h2 = -L.a[1] * E;
  // [H, rho]_01 = H_00 r01 + H_01 r11 + H_02 r21 - (r00 H_01 + r01 H_11 + r02 H_21)
  //            = h1 p1 + h2 conj(c12) - p0 h1 - c01 w1
  const cplx c12c(s.c12.re, -s.c12.im);
  cplx k01 = cplx(h1 * (s.p1 - s.p0), 0.0) + h2 * c12c - L.w[0] * s.c01;
  double d01 = k01.im;  // Re(-i k) = Im k
  if (L.T2 > 0) d01 -= s.c01.re / L.T2;
  double J = L.cP[0] * d01;
  if (L.n > 1) {
    // [H, rho]_02 = h1 c12 + h2 p2 - p0 h2 - c02 w2
    cplx k02 = h1 * s.c12 + cplx(h2 * (s.p2 - s.p0), 0.0) - L.w[1] * s.c02;
    double d02 = k02.im;
    if (L.T2 > 0) d02 -= s.c02.re / L.T2;
    J += L.cP[1] * d02;
  }
  return J;
}

// losses over the time t: ionization rates g0, g1, g2 (population of level k decays as exp(-g_k t), coherence kl as
// exp(-(g_k + g_l) t / 2)), T1 (excited -> ground), T2 (pure dephasing of the coherences)
GEM_HD inline void blochLoss(const BlochLines& L, BlochNode& s, double t, double g0, double g1, double g2) {
  if (L.T1 > 0) {
    const double e = exp(-t / L.T1), d1 = s.p1 * (1.0 - e), d2 = s.p2 * (1.0 - e);
    s.p1 -= d1; s.p2 -= d2; s.p0 += d1 + d2;
    const double c = exp(-0.5 * t / L.T1);
    s.c01 = c * s.c01; s.c02 = c * s.c02; s.c12 = (c * c) * s.c12;
  }
  if (L.T2 > 0) {
    const double c = exp(-t / L.T2);
    s.c01 = c * s.c01; s.c02 = c * s.c02; s.c12 = c * s.c12;
  }
  if (g0 > 0 || g1 > 0 || g2 > 0) {
    const double e0 = exp(-0.5 * g0 * t), e1 = exp(-0.5 * g1 * t), e2 = exp(-0.5 * g2 * t);
    s.p0 *= e0 * e0; s.p1 *= e1 * e1; s.p2 *= e2 * e2;
    s.c01 = (e0 * e1) * s.c01; s.c02 = (e0 * e2) * s.c02; s.c12 = (e1 * e2) * s.c12;
  }
}

// unitary Cayley step with the field Em (midpoint) over the time h
GEM_HD inline void blochRotate(const BlochLines& L, BlochNode& s, double Em, double h) {
  // K = (h/2) H (real symmetric); A = 1 + i K; U = 2 A^-1 - 1
  const double q = 0.5 * h;
  // coupling of line j scaled by g = tan(x) / (x cos x), x = w_j h / 2: the Cayley step of the pre-warped Hamiltonian
  // turns a weak coupling by (x / tan x) of the exact one, and the midpoint field (E0 + E1)/2 of a carrier at w_j is
  // cos(x) of its amplitude; with g a resonant drive is exact to this order (Rabi flopping of the test: error
  // 0.023 -> see README at h = lambda/32)
  auto gfac = [&](double w) { const double x = q * w; return x > 1e-8 ? tan(x) / (x * cos(x)) : 1.0; };
  const double k01 = -q * L.a[0] * gfac(L.w[0]) * Em, k02 = (L.n > 1) ? -q * L.a[1] * gfac(L.w[1]) * Em : 0.0;
  // the free precession is pre-warped, w -> (2/h) tan(w h/2): the Cayley step then turns it by exactly w h (without
  // the warp the resonance of the scheme would sit at (2/h) atan(w h/2), e.g. 0.3 % low at w h = 0.2, which near a
  // resonance changes the susceptibility by factors)
  const double k11 = tan(q * L.w[0]), k22 = (L.n > 1) ? tan(q * L.w[1]) : 0.0;
  // A = [[1, i k01, i k02], [i k01, 1 + i k11, 0], [i k02, 0, 1 + i k22]]
  const cplx A00(1.0, 0.0), A01(0.0, k01), A02(0.0, k02), A11(1.0, k11), A22(1.0, k22);
  // inverse by cofactors (A symmetric, A12 = 0)
  const cplx C00 = A11 * A22;
  const cplx C01 = cplx() - A01 * A22;  // -(A01 A22 - A02 * 0)
  const cplx C02 = cplx() - A02 * A11;  // (A01 * 0 - A02 A11)
  const cplx C11 = A00 * A22 - A02 * A02;
  const cplx C12 = A01 * A02;           // -(A00 * 0 - A01 A02)
  const cplx C22 = A00 * A11 - A01 * A01;
  const cplx det = A00 * C00 + A01 * C01 + A02 * C02;
  const double dn = det.re * det.re + det.im * det.im;
  const cplx inv(det.re / dn, -det.im / dn);
  // U = 2 A^-1 - 1 (symmetric)
  const cplx U00 = 2.0 * (C00 * inv) - cplx(1.0, 0.0), U01 = 2.0 * (C01 * inv), U02 = 2.0 * (C02 * inv);
  const cplx U11 = 2.0 * (C11 * inv) - cplx(1.0, 0.0), U12 = 2.0 * (C12 * inv);
  const cplx U22 = 2.0 * (C22 * inv) - cplx(1.0, 0.0);
  const cplx U[3][3] = {{U00, U01, U02}, {U01, U11, U12}, {U02, U12, U22}};
  // rho (Hermitian)
  auto cj = [](const cplx& z) { return cplx(z.re, -z.im); };
  const cplx R[3][3] = {{cplx(s.p0, 0.0), s.c01, s.c02}, {cj(s.c01), cplx(s.p1, 0.0), s.c12}, {cj(s.c02), cj(s.c12), cplx(s.p2, 0.0)}};
  cplx T[3][3];  // T = U rho
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) T[i][j] = U[i][0] * R[0][j] + U[i][1] * R[1][j] + U[i][2] * R[2][j];
  // rho' = T U^+ ; (U^+)_kj = conj(U_jk)
  auto el = [&](int i, int j) { return T[i][0] * cj(U[j][0]) + T[i][1] * cj(U[j][1]) + T[i][2] * cj(U[j][2]); };
  s.p0 = el(0, 0).re; s.p1 = el(1, 1).re; s.p2 = el(2, 2).re;
  s.c01 = el(0, 1); s.c02 = el(0, 2); s.c12 = el(1, 2);
}

// one half-cell: losses h/2, rotation with (E0 + E1)/2, losses h/2
GEM_HD inline void blochStep(const BlochLines& L, BlochNode& s, double E0, double E1, double h, double g0 = 0.0,
                             double g1 = 0.0, double g2 = 0.0) {
  blochLoss(L, s, 0.5 * h, g0, g1, g2);
  blochRotate(L, s, 0.5 * (E0 + E1), h);
  blochLoss(L, s, 0.5 * h, g0, g1, g2);
}

}  // namespace gem
