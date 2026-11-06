// ============================================================================
//  radiation.hpp -- radiation reaction of the beam list (numerics.radiation_reaction)
//
//  Classical synchrotron emission (chi << 1), local constant crossed field, of a particle with unit charge and
//  mass M = 1/|qm| (electrons, positrons, muons, protons; heavier particles radiate ~M^-3 less). Code units: t in
//  1/omega0, fields in m c omega0 / e, energies in M c^2:
//    chi    = K qm^2 gamma sqrt(w),   w = (E + v x B)^2 - (v.E)^2,   K = hbar omega0 / (m_e c^2)
//    Lambda = (5 / (2 sqrt 3)) alpha (M / K) chi dt / gamma                 (mean number of photons of the step)
//    m1     = (4 sqrt 3 / 15) chi gamma,   m2 = (11/12) chi^2 gamma^2         (one-photon moments)
//  Lambda m1 = (2/3) alpha (M/K) chi^2 dt is the Landau-Lifshitz loss, Lambda m2 = (55 / (24 sqrt 3)) alpha (M/K)
//  chi^3 gamma dt the quantum excitation (Sokolov-Ternov).
//    classical   : S = Lambda m1 (deterministic drag, leading Landau-Lifshitz term)
//    stochastic  : S = 0 with probability exp(-Lambda), else Gamma(k, theta) matched to the conditional mean r1 and
//                  variance s2 of the compound-Poisson sum (zero-inflated Gamma, method of N. Elkina and
//                  A. Pukhov): exact mean and variance, S >= 0. The classical synchrotron spectrum is itself close
//                  to a Gamma law (k = 0.30 for one photon), so the shape is good for any Lambda.
//    gaussian    : S = Lambda m1 + sqrt(Lambda m2) N(0, 1) (Fokker-Planck; can gain energy; for comparison)
//  The loss S is removed along u after the Lorentz push of a sub-step; chi from u before it and from the fields at
//  the point (external fields not averaged over the sub-step: radiation goes with |F|^2, the kick with <F>).
//  The random numbers are a hash of the particle state at the start of the sub-step (x, y, zeta, u_z, z):
//  independent of the order, thread, backend and MPI rank.
// ============================================================================
#pragma once
#include <math.h>
#include <stdint.h>

#include "hd.hpp"
#include "particles.hpp"
#include "rng.hpp"

namespace gem {

enum { RAD_OFF = 0, RAD_CLASSICAL = 1, RAD_STOCHASTIC = 2, RAD_GAUSSIAN = 3 };

struct RadParams {
  int model = RAD_OFF;
  double K = 0.0;          // hbar omega0 / (m_e c^2)
  uint64_t seed = 0;
};

constexpr double RAD_ALPHA = 7.2973525693e-3;

struct RadStep { double lam, m1, m2; };

GEM_HD inline RadStep radStep(double gamma, double chi, double M, double K, double dt) {
  RadStep s;
  s.lam = 1.4433756729740645 * RAD_ALPHA * (M / K) * chi * dt / gamma;   // 5 / (2 sqrt 3)
  s.m1 = 0.46188021535170065 * chi * gamma;                               // 4 sqrt 3 / 15
  s.m2 = (11.0 / 12.0) * chi * chi * gamma * gamma;
  return s;
}

GEM_HD inline double radUniform(uint64_t key, uint64_t d) { return uniform01(key + d * 0x9E3779B97F4A7C15ULL); }

GEM_HD inline double radNormal(uint64_t key, uint64_t d) {
  const double u1 = fmax(radUniform(key, d), 1e-300), u2 = radUniform(key, d + 1);
  return sqrt(-2.0 * log(u1)) * cos(6.283185307179586 * u2);
}

// Gamma(k, 1), Marsaglia & Tsang (2000); k < 1 via Gamma(k + 1) U^(1/k); acceptance > 0.95 (cap 64 trials)
GEM_HD inline double radGamma(double k, uint64_t key, uint64_t d0) {
  double boost = 1.0;
  if (k < 1.0) {
    boost = exp(log(fmax(radUniform(key, d0), 1e-300)) / k);
    k += 1.0;
  }
  const double d = k - 1.0 / 3.0, c = 1.0 / sqrt(9.0 * d);
  uint64_t n = d0 + 1;
  for (int it = 0; it < 64; ++it, n += 3) {
    const double x = radNormal(key, n);
    double v = 1.0 + c * x;
    if (v <= 0.0) continue;
    v = v * v * v;
    const double u = fmax(radUniform(key, n + 2), 1e-300);
    if (log(u) < 0.5 * x * x + d - d * v + d * log(v)) return boost * d * v;
  }
  return boost * k;
}

// energy loss of one step (units M c^2)
GEM_HD inline double radLoss(const RadStep& s, int model, uint64_t key) {
  if (!(s.lam > 0.0)) return 0.0;
  if (model == RAD_CLASSICAL) return s.lam * s.m1;
  if (model == RAD_GAUSSIAN) return s.lam * s.m1 + sqrt(s.lam * s.m2) * radNormal(key, 0);
  const double e = exp(-s.lam);
  if (radUniform(key, 0) < e) return 0.0;  // no photon
  const double q = -expm1(-s.lam);
  const double r1 = s.lam * s.m1 / q;
  const double s2 = (q * s.lam * s.m2 - s.lam * s.lam * s.m1 * s.m1 * e) / (q * q);
  if (!(s2 > 0.0)) return r1;
  return (s2 / r1) * radGamma(r1 * r1 / s2, key, 1);
}

GEM_HD inline uint64_t radBits(double x) {
  union { double d; uint64_t u; } b;
  b.d = x;
  return b.u;
}

// radiation of one sub-step: u0 = momentum before the Lorentz push (chi, Lambda), p = particle after it; the key
// hashes the state at the start of the sub-step (x0, y0, zeta0, z0)
GEM_HD inline void radiate(const RadParams& R, Particle& p, double ux0, double uy0, double uz0, double x0, double y0,
                           double zeta0, double z0, const EB& f, double dt) {
  const double aq = fabs(p.qm);
  if (R.model == RAD_OFF || aq == 0.0) return;
  const double g0 = sqrt(1.0 + ux0 * ux0 + uy0 * uy0 + uz0 * uz0);
  const double vx = ux0 / g0, vy = uy0 / g0, vz = uz0 / g0;
  const double Fx = f.Ex + vy * f.Bz - vz * f.By, Fy = f.Ey + vz * f.Bx - vx * f.Bz, Fz = f.Ez + vx * f.By - vy * f.Bx;
  const double ve = vx * f.Ex + vy * f.Ey + vz * f.Ez;
  const double w = fmax(Fx * Fx + Fy * Fy + Fz * Fz - ve * ve, 0.0);
  const double chi = R.K * aq * aq * g0 * sqrt(w);
  const RadStep st = radStep(g0, chi, 1.0 / aq, R.K, fabs(dt));
  const uint64_t key =
      mix64(R.seed ^ mix64(radBits(x0) ^ mix64(radBits(y0) ^ mix64(radBits(zeta0) ^ mix64(radBits(uz0) ^ mix64(radBits(z0)))))));
  const double S = radLoss(st, R.model, key);
  if (S == 0.0) return;
  const double p2 = p.ux * p.ux + p.uy * p.uy + p.uz * p.uz;
  if (!(p2 > 0.0)) return;
  double gn = sqrt(1.0 + p2) - S;
  if (gn < 1.0) gn = 1.0;
  const double s = sqrt((gn - 1.0) * (gn + 1.0) / p2);
  p.ux *= s; p.uy *= s; p.uz *= s;
}

}  // namespace gem
