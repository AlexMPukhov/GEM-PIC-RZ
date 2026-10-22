// ============================================================================
//  simulation_fields.cpp -- beam fields, field gather, smoothing and filtering of currents, beam deposition
// ============================================================================
#include "simulation.hpp"

#include "eigen.hpp"
#include "rng.hpp"

#include <algorithm>
#if defined(_OPENMP)
#include <omp.h>
#endif
#include <cmath>

namespace gem {

static const cd I1(0.0, 1.0);

// adds sign x the co-moving self-field of the beam list to fold, per group of particles
// (bunches, and with deckGamma = false also trapped particles), with the gamma of the
// bunch block (deckGamma) or the mean gamma of the group
void Simulation::addBeamFields(double sign, bool deckGamma) {
  const int N = g.N;
  // radial eigenmodes of -Lap^c_m (computed on demand, per mode)
  std::vector<std::vector<double>> lam(M + 1), Q(M + 1);
  auto eigen = [&](int m) {
    if (!lam[m].empty()) return;
    const Tridiag& T = solver.poissonMatrix(m);
    std::vector<double> d(N), e(N, 0.0);
    for (int i = 0; i < N; ++i) d[i] = T.b[i];
    for (int i = 0; i + 1 < N; ++i) e[i] = -std::sqrt(T.cc[i] * T.a[i + 1]);  // symmetrized off-diagonal
    if (!symTridiagEigen(d, e, Q[m])) throw std::runtime_error("bunch fields: eigenvalue iteration failed");
    lam[m] = d;
  };
  Currents C;
  C.init(M, Nz, N);
  std::vector<cd> rho(size_t(Nz) * N), psi(size_t(Nz) * N), col(Nz), Er(N + 1), Ep(N);
  std::vector<double> sq(N);
  for (int i = 0; i < N; ++i) sq[i] = std::sqrt(g.Ac[i]);
  std::vector<int> groups;
  for (size_t b = 0; b < cfg.bunches.size(); ++b) groups.push_back(int(b));
  if (!deckGamma) {
    for (const auto& p : beam)
      if (p.sp < 0 && std::find(groups.begin(), groups.end(), p.sp) == groups.end()) groups.push_back(p.sp);
  }
  for (const int b : groups) {
    C.zero();
    bool any = false;
    double gsum = 0.0, wsum = 0.0;
    for (const auto& p : beam)
      if (p.sp == b) {
        depositBeamParticle(C, p);
        any = true;
        gsum += std::fabs(p.q) * gammaOf(p);
        wsum += std::fabs(p.q);
      }
    if (!any) continue;
    filterBeamCurrent(C);
    const double gam = deckGamma ? cfg.bunches[size_t(b)].gamma : gsum / wsum;
    const double beta = (gam > 1e5) ? 1.0 : std::sqrt(1.0 - 1.0 / (gam * gam));
    const double ig2 = (gam > 1e5) ? 0.0 : 1.0 / (gam * gam);
    for (int m = 0; m <= M; ++m) {
      // charge density at slice n: linear interpolation of the half-cells n-1 and n
      for (int n = 0; n < Nz; ++n)
        for (int i = 0; i < N; ++i)
          rho[size_t(n) * N + i] = 0.5 * (zv.wLo(n) * ((n > 0) ? C.jz[C.c(m, n - 1, i)] : cd(0.0)) +
                                          zv.wHi(n) * C.jz[C.c(m, n, i)]) / beta;
      if (ig2 == 0.0) {
        for (int n = 0; n < Nz; ++n) solver.poisson(m, &rho[size_t(n) * N], &psi[size_t(n) * N]);
      } else {
        eigen(m);
        const std::vector<double>& Z = Q[m];
        // rho_k(n) = sum_i Q_ik sqrt(A_i) rho_i(n)
        std::vector<cd> rk(size_t(Nz) * N);
        for (int n = 0; n < Nz; ++n) {
          const cd* r = &rho[size_t(n) * N];
          cd* out = &rk[size_t(n) * N];
          for (int i = 0; i < N; ++i) {
            const cd v = sq[i] * r[i];
            if (v == cd(0.0)) continue;
            const double* zi = &Z[size_t(i) * N];
            for (int k = 0; k < N; ++k) out[k] += zi[k] * v;
          }
        }
        // tridiagonal solve along zeta for every eigenmode
        std::vector<double> lo(Nz), di(Nz), up(Nz), cp(Nz);
        for (int n = 0; n < Nz; ++n) {
          const double hm = (n > 0) ? zv.len(n - 1) : zv.len(0), hp = zv.len(n);
          const double f = 2.0 * ig2 / (hm + hp);
          lo[n] = -f / hm; up[n] = -f / hp; di[n] = f / hm + f / hp;
        }
        for (int k = 0; k < N; ++k) {
          const double lk = lam[m][k];
          for (int n = 0; n < Nz; ++n) col[n] = rk[size_t(n) * N + k];
          double den = di[0] + lk;
          cp[0] = up[0] / den;
          col[0] /= den;
          for (int n = 1; n < Nz; ++n) {
            den = di[n] + lk - lo[n] * cp[n - 1];
            cp[n] = up[n] / den;
            col[n] = (col[n] - lo[n] * col[n - 1]) / den;
          }
          for (int n = Nz - 2; n >= 0; --n) col[n] -= cp[n] * col[n + 1];
          for (int n = 0; n < Nz; ++n) rk[size_t(n) * N + k] = col[n];
        }
        // psi_i(n) = sum_k Q_ik psi_k(n) / sqrt(A_i)
        for (int n = 0; n < Nz; ++n) {
          const cd* in = &rk[size_t(n) * N];
          cd* out = &psi[size_t(n) * N];
          for (int i = 0; i < N; ++i) {
            const double* zi = &Z[size_t(i) * N];
            cd a = 0.0;
            for (int k = 0; k < N; ++k) a += zi[k] * in[k];
            out[i] = a / sq[i];
          }
        }
      }
      for (int n = 0; n < Nz; ++n) {
        solver.gradient(m, &psi[size_t(n) * N], Er.data(), Ep.data());
        for (int i = 1; i <= N; ++i) {
          fold.Rr[fold.f(m, n, i)] += sign * ((1.0 + beta) * Er[i]);
          if (beta < 1.0) fold.Lr[fold.f(m, n, i)] += sign * ((1.0 - beta) * Er[i]);
        }
        for (int i = 0; i < N; ++i) {
          fold.Rp[fold.c(m, n, i)] += sign * ((1.0 + beta) * Ep[i]);
          if (beta < 1.0) fold.Lp[fold.c(m, n, i)] += sign * ((1.0 - beta) * Ep[i]);
        }
        if (ig2 > 0.0) {  // E_z = -gamma^-2 d_zeta psi (central difference, psi = 0 beyond the window)
          const double hm = (n > 0) ? zv.len(n - 1) : zv.len(0), hp = zv.len(n);
          for (int i = 0; i < N; ++i) {
            const cd up1 = (n + 1 < Nz) ? psi[size_t(n + 1) * N + i] : cd(0.0);
            const cd dn1 = (n > 0) ? psi[size_t(n - 1) * N + i] : cd(0.0);
            fold.Ez[fold.c(m, n, i)] -= sign * (ig2 * (up1 - dn1) / (hm + hp));
          }
        }
      }
    }
  }
}

// ============================================================================
//  field interpolation
// ============================================================================
// add w * (mode-m field components at radius r, slice n) to the accumulators
void Simulation::accumulate(const Fields& F, int m, int n, double w, double r, const RadialGrid::W& wc,
                            const RadialGrid::W& wf, cd& er, cd& ep, cd& ez, cd& br, cd& bp, cd& bz) const {
  const int N = g.N;
  const size_t cb = F.c(m, n, 0), fb = F.f(m, n, 0);
  auto EpC = [&](int i) { return 0.5 * (F.Rp[cb + i] + F.Lp[cb + i]); };
  auto BrC = [&](int i) { return 0.5 * (F.Lp[cb + i] - F.Rp[cb + i]); };
  auto EzC = [&](int i) { return F.Ez[cb + i]; };
  // linear extrapolation of a centre quantity to the axis
  const double c0 = g.c[0], c1 = g.c[1];
  auto axis = [&](cd f0, cd f1) { return (c1 * f0 - c0 * f1) / (c1 - c0); };
  cd ep0 = 0.0, br0 = 0.0;
  if (m == 1) { ep0 = axis(EpC(0), EpC(1)); br0 = axis(BrC(0), BrC(1)); }

  // centres
  cd vep, vbr, vez;
  if (wc.mode == 0) {
    vep = wc.wa * EpC(wc.a) + wc.wb * EpC(wc.a + 1);
    vbr = wc.wa * BrC(wc.a) + wc.wb * BrC(wc.a + 1);
    vez = wc.wa * EzC(wc.a) + wc.wb * EzC(wc.a + 1);
  } else if (wc.mode == 2) {
    vep = EpC(N - 1); vbr = BrC(N - 1); vez = EzC(N - 1);
  } else {  // r < c_0: regular behaviour on the axis
    double x = r / c0;
    vez = EzC(0) * ((m == 0) ? 1.0 : std::pow(x, m));
    if (m == 0) { vep = EpC(0) * x; vbr = BrC(0) * x; }
    else if (m == 1) { vep = ep0 + (EpC(0) - ep0) * x; vbr = br0 + (BrC(0) - br0) * x; }
    else { double fx = std::pow(x, m - 1); vep = EpC(0) * fx; vbr = BrC(0) * fx; }
  }
  // faces (face 0 for m = 1 from regularity: E_+ = B_+ = 0 on the axis)
  auto ErF = [&](int i) -> cd {
    if (i == 0) return (m == 1) ? -I1 * ep0 : cd(0.0);
    return 0.5 * (F.Rr[fb + i] + F.Lr[fb + i]);
  };
  auto BpF = [&](int i) -> cd {
    if (i == 0) return (m == 1) ? I1 * br0 : cd(0.0);
    return 0.5 * (F.Rr[fb + i] - F.Lr[fb + i]);
  };
  auto BzF = [&](int i) { return F.Bz[fb + i]; };
  cd ver, vbp, vbz;
  if (wf.mode == 0) {
    ver = wf.wa * ErF(wf.a) + wf.wb * ErF(wf.a + 1);
    vbp = wf.wa * BpF(wf.a) + wf.wb * BpF(wf.a + 1);
    vbz = wf.wa * BzF(wf.a) + wf.wb * BzF(wf.a + 1);
  } else {
    ver = ErF(N); vbp = BpF(N); vbz = BzF(N);
  }
  er += w * ver; ep += w * vep; ez += w * vez; br += w * vbr; bp += w * vbp; bz += w * vbz;
}

// Cartesian fields at (x,y), weighted combination of slices n1, n2
EB Simulation::gather(const Fields& F, int n1, double w1, int n2, double w2, double x, double y) const {
  double r = std::sqrt(x * x + y * y);
  double cs = 1.0, sn = 0.0;
  if (r > 0.0) { cs = x / r; sn = y / r; }
  auto wc = g.locC_lin(r), wf = g.locF_lin(r);  // gather: linear in r
  const cd eip(cs, sn);
  cd em(1.0, 0.0);
  double Er = 0, Ep = 0, Ez = 0, Br = 0, Bp = 0, Bz = 0;
  for (int m = 0; m <= M; ++m) {
    cd er = 0, ep = 0, ez = 0, br = 0, bp = 0, bz = 0;
    if (w1 != 0.0) accumulate(F, m, n1, w1, r, wc, wf, er, ep, ez, br, bp, bz);
    if (w2 != 0.0) accumulate(F, m, n2, w2, r, wc, wf, er, ep, ez, br, bp, bz);
    const double fac = (m == 0) ? 1.0 : 2.0;
    Er += fac * std::real(er * em); Ep += fac * std::real(ep * em);
    Ez += fac * std::real(ez * em); Br += fac * std::real(br * em);
    Bp += fac * std::real(bp * em); Bz += fac * std::real(bz * em);
    em *= eip;
  }
  return {Er * cs - Ep * sn, Er * sn + Ep * cs, Ez, Br * cs - Bp * sn, Br * sn + Bp * cs, Bz};
}

// ============================================================================
//  deposition
// ============================================================================
// current density j_m = sum W v S e^{-i m phi} / (2 pi V h_k) into half-cell k
void Simulation::deposit(Currents& C, int k, double wk, double x, double y, double vx, double vy, double vz,
                         double W, bool ion, int kc) {
  if (k < 0 || k >= Nz || wk == 0.0 || W == 0.0) return;
  const int kd = (kc < 0) ? k : kc;  // index in C
  const int N = g.N;
  double r = std::sqrt(x * x + y * y);
  double cs = 1.0, sn = 0.0;
  if (r > 0.0) { cs = x / r; sn = y / r; }
  const double vr = cs * vx + sn * vy, vp = -sn * vx + cs * vy;
  auto wc = g.locC(r), wf = g.locF(r);
  const double norm = W * wk / (2.0 * PI * zv.len(k));
  const cd emi(cs, -sn);
  cd em(1.0, 0.0);
  for (int m = 0; m <= M; ++m) {
    // centres: j_phi, j_z
    auto addC = [&](int i, double sp, double sz) {
      double f = norm / g.Vc[i];
      C.jp[C.c(m, kd, i)] += f * sp * vp * em;
      C.jz[C.c(m, kd, i)] += f * sz * vz * em;
    };
    if (wc.mode == 0) { addC(wc.a, wc.wa, wc.wa); addC(wc.a + 1, wc.wb, wc.wb); }
    else if (wc.mode == 2) addC(N - 1, 1.0, 1.0);
    else addC(0, 1.0, (m == 0) ? 1.0 : std::pow(r / g.c[0], m));
    // faces: j_r (face 0 is not used by the solver)
    auto addF = [&](int i, double s) {
      if (i == 0) return;
      C.jr[C.f(m, kd, i)] += norm * s / g.Vf[i] * vr * em;
    };
    if (wf.mode == 0) { addF(wf.a, wf.wa); addF(wf.a + 1, wf.wb); }
    else addF(N, 1.0);
    em *= emi;
  }
  // m = 0 charge density (diagnostic), separately for ions
  std::vector<double>& dens = ion ? C.rhoi : C.rho;
  auto addR = [&](int i, double s) { dens[size_t(kd) * N + i] += norm * s / g.Vc[i]; };
  if (wc.mode == 0) { addR(wc.a, wc.wa); addR(wc.a + 1, wc.wb); }
  else if (wc.mode == 2) addR(N - 1, 1.0);
  else addR(0, 1.0);
}

// plasma.smooth_length: the currents and densities that the streaming particles (plasma
// electrons, ions, freed electrons) deposited into half-cell n are smoothed radially,
// (1 - a^2 Lap)^-1, per mode (smoothing.hpp). The beam list is not smoothed, nor are the
// fields acting on the particles.
void Simulation::smoothStream(Currents& C, int n) {
  if (!sfilter.on()) return;
  const int N = g.N;
  std::vector<cd> jr(size_t(N) + 1, cd(0.0)), jp(size_t(N), cd(0.0)), jz(size_t(N), cd(0.0));
  for (int m = 0; m <= M; ++m) {
    for (int i = 0; i <= N; ++i) jr[i] = C.jr[C.f(m, n, i)];
    for (int i = 0; i < N; ++i) { jp[i] = C.jp[C.c(m, n, i)]; jz[i] = C.jz[C.c(m, n, i)]; }
    sfilter.vector(m, jr.data(), jp.data());
    sfilter.scalar(m, jz.data());
    for (int i = 1; i <= N; ++i) C.jr[C.f(m, n, i)] = jr[i];
    for (int i = 0; i < N; ++i) { C.jp[C.c(m, n, i)] = jp[i]; C.jz[C.c(m, n, i)] = jz[i]; }
  }
  sfilter.scalar0(&C.rho[size_t(n) * N]);
  sfilter.scalar0(&C.rhoi[size_t(n) * N]);
}

// Numerical Cherenkov: the dispersion relation of the scheme has omega = 0 at
// the zeta-Nyquist wavenumber (kappa h = pi) for every transverse k, so this
// odd-even mode is resonant with every v ~ c source. Beam particles keep their
// zeta from step to step, so the Nyquist part of their current (shot noise)
// drives it secularly. Passes of the (1,2,1)/4 filter remove that part
// (response cos^2(kappa h/2)); zero beyond the ends of the window.
void Simulation::filterBeamCurrent(Currents& C) const {
  const int passes = cfg.numerics.beam_filter;
  if (passes <= 0) return;
  auto column = [&](std::vector<cd>& v, int L) {
    for (int m = 0; m <= M; ++m)
      for (int i = 0; i < L; ++i)
        for (int p = 0; p < passes; ++p) {
          cd prev(0.0);
          for (int k = 0; k < Nz; ++k) {
            cd& x = v[(size_t(m) * Nz + k) * L + i];
            const cd cur = x;
            const cd next = (k + 1 < Nz) ? v[(size_t(m) * Nz + k + 1) * L + i] : cd(0.0);
            x = 0.25 * prev + 0.5 * cur + 0.25 * next;
            prev = cur;
          }
        }
  };
  column(C.jr, g.N + 1);
  column(C.jp, g.N);
  column(C.jz, g.N);
}

// beam-type deposit (current at level s^1): W = q / v_z, linear in zeta
// between half-cell centres
void Simulation::depositBeamParticle(Currents& C, const Particle& p) {
  const int nsub = std::max(1, int(std::ceil(Delta / cfg.numerics.beam_dz_max)));
  double ux, uy, uz;
  velocityAtPosition(ext, p, Delta / nsub, ux, uy, uz);
  const double gm = std::sqrt(1.0 + ux * ux + uy * uy + uz * uz);
  double vx = ux / gm, vy = uy / gm, vz = uz / gm;
  if (vz <= 0.0) return;
  double f;
  const int j = zv.beamCell(p.zeta, f);
  double W = p.q / vz;
  const bool ion = isIon(p.sp);
  deposit(C, j - 1, 1.0 - f, p.x, p.y, vx, vy, vz, W, ion);
  deposit(C, j, f, p.x, p.y, vx, vy, vz, W, ion);
}

}  // namespace gem
