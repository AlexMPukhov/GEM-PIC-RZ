// ============================================================================
//  bound.cpp -- bound-electron medium (see bound.hpp, bloch.hpp)
// ============================================================================
#include "bound.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <stdexcept>

namespace gem {

namespace {
constexpr double ALPHA_INV = 137.035999084;
constexpr double A0_UM = 5.29177210903e-5;  // Bohr radius [um]
}  // namespace

void BoundMedium::setup(const BoundElectronsConfig& c, const RadialGrid& grid, int M_, const Units& u, double T) {
  const double lambda0_um = u.lambda0_um;
  on = c.present;
  if (!on) return;
  cfg = c;
  g = &grid;
  N = grid.N;
  M = M_;
  nphi = c.n_phi > 0 ? c.n_phi : 4 * (M + 1);
  L.n = c.lines;
  const double a0k0 = A0_UM * 2.0 * M_PI / lambda0_um;
  for (int j = 0; j < c.lines; ++j) {
    L.w[j] = c.omega[j];
    L.a[j] = c.dipole[j] * ALPHA_INV;
    L.cP[j] = 2.0 * c.dipole[j] * a0k0;
  }
  L.T1 = c.T1;
  L.T2 = c.T2;
  ion = c.ionIndex >= 0;
  if (ion) {
    adkL[0].setup({c.ip_ground_eV}, 0, 1, u, T);
    for (int j = 0; j < c.lines; ++j) adkL[j + 1].setup({c.ip_excited_eV[j]}, 0, 1, u, T);
  }
  if (M > 63) throw std::runtime_error("bound_electrons: m_max must be < 64");
  cs.resize(nphi);
  sn.resize(nphi);
  for (int k = 0; k < nphi; ++k) {
    cs[k] = std::cos(2.0 * M_PI * k / nphi);
    sn[k] = std::sin(2.0 * M_PI * k / nphi);
  }
  const size_t nn = size_t(N) * nphi;
  st0.resize(nn); st1.resize(nn); tmp1.resize(nn);
  Ex0.assign(nn, 0.0); Ex1.assign(nn, 0.0); Exa.assign(nn, 0.0);
  Jx.assign(nn, 0.0); Hx.assign(nn, 0.0); Jr0.assign(nn, 0.0); Jr1.assign(nn, 0.0);
  for (auto* v : {&J0r, &J1r, &Hr}) v->assign(size_t(M + 1) * (N + 1), 0.0);
  for (auto* v : {&J0p, &J1p, &Hp}) v->assign(size_t(M + 1) * N, 0.0);
  na0.assign(N, 0.0); na1.assign(N, 0.0);
}

void BoundMedium::begin(const std::vector<double>& a0, const std::vector<double>& a1, int K) {
  nc = 0;
  for (int i = 0; i < N && i < K; ++i) {
    na0[i] = a0[i];
    na1[i] = a1[i];
    if (na0[i] > 0 || na1[i] > 0) nc = i + 1;
  }
  for (int i = nc; i < N; ++i) na0[i] = na1[i] = 0.0;
  for (auto& s : st0) blochGround(s);
  for (auto& s : st1) blochGround(s);
  std::fill(Ex0.begin(), Ex0.end(), 0.0);  // no field and no current ahead of the window
  std::fill(Jr0.begin(), Jr0.end(), 0.0);
  std::fill(Jr1.begin(), Jr1.end(), 0.0);
  std::fill(Ex1.begin(), Ex1.end(), 0.0);
}

// E_x at the centres: E_r interpolated from the faces, E_phi of the centre; F = F_0 + 2 Re sum F_m e^{i m phi}
void BoundMedium::fieldAt(const Fields& F, int n, std::vector<double>& Ex) const {
  for (int i = 0; i < nc; ++i) {
    const double w = (g->c[i] - g->rho[i]) / (g->rho[i + 1] - g->rho[i]);
    cd er[64], ep[64];
    for (int m = 0; m <= M; ++m) {
      const size_t f = F.f(m, n, 0), c = F.c(m, n, 0);
      ep[m] = 0.5 * (F.Rp[c + i] + F.Lp[c + i]);
      // face 0 (the axis) carries no E_r; there E_r = 0 for m != 1 and, by regularity of E_x + i E_y,
      // E_r,1 = -i E_phi,1 for m = 1 (the laser field is finite on the axis)
      const cd e0 = (i > 0) ? 0.5 * (F.Rr[f + i] + F.Lr[f + i]) : (m == 1 ? cd(0.0, -1.0) * ep[m] : cd(0.0));
      const cd e1 = 0.5 * (F.Rr[f + i + 1] + F.Lr[f + i + 1]);
      er[m] = (1.0 - w) * e0 + w * e1;
    }
    for (int k = 0; k < nphi; ++k) {
      double Er = er[0].real(), Ep = ep[0].real();
      const cd e(cs[k], sn[k]);
      cd em = e;
      for (int m = 1; m <= M; ++m) {
        Er += 2.0 * (er[m] * em).real();
        Ep += 2.0 * (ep[m] * em).real();
        em *= e;
      }
      Ex[size_t(i) * nphi + k] = Er * cs[k] - Ep * sn[k];
    }
  }
}

void BoundMedium::project(const std::vector<double>& x, cd* fr, cd* fp) const {
  // centre values of J_r = J_x cos, J_phi = -J_x sin, mode m: (1/n_phi) sum F e^{-i m phi}
  std::vector<cd> jrc(size_t(M + 1) * N, 0.0);
  for (size_t j = 0; j < size_t(M + 1) * N; ++j) fp[j] = 0.0;
  for (size_t j = 0; j < size_t(M + 1) * (N + 1); ++j) fr[j] = 0.0;
  for (int i = 0; i < nc; ++i) {
    for (int k = 0; k < nphi; ++k) {
      const double v = x[size_t(i) * nphi + k];
      if (v == 0.0) continue;
      const cd e(cs[k], -sn[k]);
      cd em(1.0, 0.0);
      for (int m = 0; m <= M; ++m) {
        jrc[size_t(m) * N + i] += (v * cs[k]) * em;
        fp[size_t(m) * N + i] += (-v * sn[k]) * em;
        em *= e;
      }
    }
  }
  const double s = 1.0 / nphi;
  for (int m = 0; m <= M; ++m) {
    for (int i = 0; i < N; ++i) {
      fp[size_t(m) * N + i] *= s;
      jrc[size_t(m) * N + i] *= s;
    }
    // faces: mean of the neighbouring centres (face 0 is not used by the solver)
    for (int i = 1; i <= N; ++i) {
      const cd a = jrc[size_t(m) * N + i - 1], b = (i < N) ? jrc[size_t(m) * N + i] : cd(0.0);
      fr[size_t(m) * (N + 1) + i] = 0.5 * (a + b);
    }
  }
}

void BoundMedium::level0(int n, double h, const Fields& fold) {
  fieldAt(fold, n, Exa);
  for (int i = 0; i < nc; ++i)
    for (int k = 0; k < nphi; ++k) {
      const size_t j = size_t(i) * nphi + k;
      if (na0[i] == 0.0) { Jx[j] = 0.0; continue; }
      const double P1 = blochP(L, st0[j]);
      double g0 = 0, g1 = 0, g2 = 0;
      if (ion) rates(Ex0[j], Exa[j], g0, g1, g2);
      blochStep(L, st0[j], Ex0[j], Exa[j], h, g0, g1, g2);
      // slice current by the trapezoidal rule over the half-cell (as the dielectric): (J_n + J_{n+1})/2 = d P / h
      Jx[j] = 2.0 * na0[i] * (blochP(L, st0[j]) - P1) / h - Jr0[j];
      Jr0[j] = Jx[j];
    }
  std::swap(Ex0, Exa);  // E of slice n = the start of the next half-cell
  project(Jx, J0r.data(), J0p.data());
}

void BoundMedium::level1(int n, double h, const Fields& fnew, bool commit) {
  fieldAt(fnew, n, Exa);
  double ex = 0.0;
  for (int i = 0; i < nc; ++i)
    for (int k = 0; k < nphi; ++k) {
      const size_t j = size_t(i) * nphi + k;
      if (na1[i] == 0.0) { Jx[j] = 0.0; Hx[j] = 0.0; continue; }
      BlochNode s = st1[j];
      const double P1 = blochP(L, s);
      double g0 = 0, g1 = 0, g2 = 0;
      if (ion) rates(Ex1[j], Exa[j], g0, g1, g2);
      blochStep(L, s, Ex1[j], Exa[j], h, g0, g1, g2);
      Hx[j] = na1[i] * (blochP(L, s) - P1) / h;
      Jx[j] = 2.0 * Hx[j] - Jr1[j];
      tmp1[j] = s;
      ex = std::max(ex, s.p1 + s.p2);
    }
  project(Jx, J1r.data(), J1p.data());
  project(Hx, Hr.data(), Hp.data());
  if (commit) {
    if (n == 0) diagRear(tmp1);
    std::swap(st1, tmp1);
    for (int i = 0; i < nc; ++i)
      for (int k = 0; k < nphi; ++k) Jr1[size_t(i) * nphi + k] = Jx[size_t(i) * nphi + k];
    for (int i = 0; i < nc; ++i)
      for (int k = 0; k < nphi; ++k) Ex1[size_t(i) * nphi + k] = Exa[size_t(i) * nphi + k];
    exMax = ex;
  }
}

void BoundMedium::rates(double E0, double E1, double& g0, double& g1, double& g2) const {
  const double a0 = std::fabs(E0), a1 = std::fabs(E1);  // |E| of the scalar model: the laser component
  g0 = 0.5 * (adkL[0].rate(0, a0) + adkL[0].rate(0, a1));
  g1 = 0.5 * (adkL[1].rate(0, a0) + adkL[1].rate(0, a1));
  g2 = L.n > 1 ? 0.5 * (adkL[2].rate(0, a0) + adkL[2].rate(0, a1)) : 0.0;
}

void BoundMedium::populations(double x, double y, double& p0, double& p1, double& p2) const {
  p0 = 1.0; p1 = 0.0; p2 = 0.0;
  const double r = std::sqrt(x * x + y * y);
  const auto it = std::upper_bound(g->rho.begin(), g->rho.end(), r);
  const int i = int(it - g->rho.begin()) - 1;
  if (i < 0 || i >= nc || na1[i] == 0.0) return;
  double ph = std::atan2(y, x);
  if (ph < 0) ph += 2.0 * M_PI;
  const int k = int(std::lround(ph / (2.0 * M_PI) * nphi)) % nphi;
  const BlochNode& q = st1[size_t(i) * nphi + k];
  const double t = q.p0 + q.p1 + q.p2;
  if (!(t > 0)) return;
  p0 = q.p0 / t; p1 = q.p1 / t; p2 = q.p2 / t;
}

// diagnostics at the rear: excited population and neutral fraction on the axis (centre 0, mean over the angles),
// area mean of the neutral fraction over the medium
void BoundMedium::diagRear(const std::vector<BlochNode>& s) {
  double e = 0.0, t = 0.0, tw = 0.0, w = 0.0;
  for (int k = 0; k < nphi; ++k) {
    e += s[size_t(k)].p1 + s[size_t(k)].p2;
    t += s[size_t(k)].p0 + s[size_t(k)].p1 + s[size_t(k)].p2;
  }
  for (int i = 0; i < nc; ++i) {
    if (na1[i] == 0.0) continue;
    const double a = na1[i] * g->Ac[i];
    for (int k = 0; k < nphi; ++k) {
      const BlochNode& q = s[size_t(i) * nphi + k];
      tw += a * (q.p0 + q.p1 + q.p2);
      w += a;
    }
  }
  exAxis = nc > 0 ? e / nphi : 0.0;
  nuAxis = nc > 0 ? t / nphi : 1.0;
  nuMean = w > 0 ? tw / w : 1.0;
}

void BoundMedium::getState(std::vector<double>& out) const {
  out.clear();
  for (const auto* v : {&st0, &st1})
    for (const auto& s : *v) {
      out.insert(out.end(), {s.p0, s.p1, s.p2, s.c01.re, s.c01.im, s.c02.re, s.c02.im, s.c12.re, s.c12.im});
    }
  for (const auto* v : {&Ex0, &Ex1, &Jr0, &Jr1}) out.insert(out.end(), v->begin(), v->end());
}

void BoundMedium::setState(const std::vector<double>& in) {
  const size_t nn = st0.size();
  if (in.size() != 2 * nn * 9 + 4 * nn) throw std::runtime_error("bound_electrons: state message of the wrong size");
  size_t o = 0;
  for (auto* v : {&st0, &st1})
    for (auto& s : *v) {
      s.p0 = in[o]; s.p1 = in[o + 1]; s.p2 = in[o + 2];
      s.c01 = cplx(in[o + 3], in[o + 4]); s.c02 = cplx(in[o + 5], in[o + 6]); s.c12 = cplx(in[o + 7], in[o + 8]);
      o += 9;
    }
  for (auto* v : {&Ex0, &Ex1, &Jr0, &Jr1}) {
    std::copy(in.begin() + long(o), in.begin() + long(o + nn), v->begin());
    o += nn;
  }
}

}  // namespace gem
