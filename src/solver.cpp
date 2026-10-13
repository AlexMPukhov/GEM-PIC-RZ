// ============================================================================
//  solver.cpp -- tridiagonal field update (verified against a dense solve of
//  all discrete equations with the verification script of the prototype)
// ============================================================================
#include "solver.hpp"

#include <map>

namespace gem {

void Tridiag::factor() {
  cp.resize(n); den.resize(n); dp.resize(n);
  double d = b[0];
  den[0] = 1.0 / d; cp[0] = cc[0] / d;
  for (int i = 1; i < n; ++i) {
    d = b[i] - a[i] * cp[i - 1];
    den[i] = 1.0 / d; cp[i] = cc[i] / d;
  }
}

void Tridiag::solve(const cd* rhs, cd* x) {
  dp[0] = rhs[i0] * den[0];
  for (int i = 1; i < n; ++i) dp[i] = (rhs[i0 + i] - a[i] * dp[i - 1]) * den[i];
  x[i0 + n - 1] = dp[n - 1];
  for (int i = n - 2; i >= 0; --i) x[i0 + i] = dp[i] - cp[i] * x[i0 + i + 1];
}

void FieldSolver::setup(const RadialGrid& grid, int M_) {
  g = &grid; M = M_; N = grid.N;
  D = -1.0; nCls = 0; hs.clear();
  for (auto* v : {&Er1, &Bp1, &Pr, &Qr, &tf1, &tf2, &rf}) v->assign(N + 1, 0.0);
  for (auto* v : {&Ep1, &Br1, &Pp, &Qp, &tc1, &tc2, &rc}) v->assign(N, 0.0);
  const auto& rho = grid.rho; const auto& c = grid.c;
  Tp.assign(M + 1, Tridiag());
  for (int m = 0; m <= M; ++m) {
    // Poisson operator -Lap^c_m on the centres (ghost psi = 0 beyond c_{N-1})
    Tridiag& Pm = Tp[m];
    Pm.i0 = 0; Pm.n = N;
    Pm.a.assign(N, 0.0); Pm.b.assign(N, 0.0); Pm.cc.assign(N, 0.0);
    for (int i = 0; i < N; ++i) {
      double a = rho[i + 1] / (c[i + 1] - c[i]) / grid.Ac[i];
      double Ld = -a, Lsub = 0, Lsup = (i + 1 < N) ? a : 0.0;
      if (i > 0) {
        double bb = rho[i] / (c[i] - c[i - 1]) / grid.Ac[i];
        Ld -= bb; Lsub = bb;
      }
      Ld -= double(m * m) / (c[i] * c[i]);
      Pm.b[i] = -Ld; Pm.a[i] = -Lsub; Pm.cc[i] = -Lsup;
    }
    Pm.factor();
  }
}

bool FieldSolver::prepare(double D_, const ZetaGrid& zg) {
  bool newWindow = (hs != zg.hs);
  if (newWindow) {
    hs = zg.hs;
    // classes of equal half-cell length (exact comparison; one class on a uniform grid)
    hCls.clear();
    cls.assign(hs.size(), 0);
    std::map<double, int> index;
    for (size_t n = 0; n < hs.size(); ++n) {
      auto it = index.find(hs[n]);
      if (it == index.end()) {
        it = index.emplace(hs[n], int(hCls.size())).first;
        hCls.push_back(hs[n]);
      }
      cls[n] = it->second;
    }
    nCls = int(hCls.size());
  }
  if (!newWindow && D_ == D) return false;
  D = D_;
  ++ver;
  const auto& rho = g->rho; const auto& c = g->c;
  Tf.assign(size_t(nCls) * (M + 1), Tridiag());
  Tc.assign(size_t(nCls) * (M + 1), Tridiag());
  for (int q = 0; q < nCls; ++q) {
    const double h = hCls[q];
    const double w = h * (2.0 * D + h);
    for (int m = 0; m <= M; ++m) {
      // --- faces: rows i0..N-1 (B_z(rho_N) = 0; for m>=1 also B_z(0) = 0)
      Tridiag& T = Tf[size_t(q) * (M + 1) + m];
      T.i0 = (m == 0) ? 0 : 1;
      T.n = N - T.i0;
      T.a.assign(T.n, 0.0); T.b.assign(T.n, 0.0); T.cc.assign(T.n, 0.0);
      for (int r = 0; r < T.n; ++r) {
        int i = T.i0 + r;
        double Ld, Lsub = 0, Lsup = 0;
        if (i == 0) {  // axis row, m = 0 only
          double a = c[0] / rho[1] / g->Af[0];
          Ld = -a; Lsup = a;
        } else {
          double a = c[i] / (rho[i + 1] - rho[i]) / g->Af[i];
          double bb = c[i - 1] / (rho[i] - rho[i - 1]) / g->Af[i];
          Ld = -a - bb - double(m * m) / (rho[i] * rho[i]);
          Lsup = (i + 1 <= N - 1) ? a : 0.0;
          Lsub = (i - 1 >= T.i0) ? bb : 0.0;
        }
        T.b[r] = 16.0 - w * Ld; T.a[r] = -w * Lsub; T.cc[r] = -w * Lsup;
      }
      T.factor();
      // --- centres: rows 0..N-1 (ghost E_z = 0 beyond c_{N-1})
      Tridiag& C = Tc[size_t(q) * (M + 1) + m];
      C.i0 = 0; C.n = N;
      C.a.assign(N, 0.0); C.b.assign(N, 0.0); C.cc.assign(N, 0.0);
      for (int i = 0; i < N; ++i) {
        double a = rho[i + 1] / (c[i + 1] - c[i]) / g->Ac[i];
        double Ld = -a, Lsub = 0, Lsup = (i + 1 < N) ? a : 0.0;
        if (i > 0) {
          double bb = rho[i] / (c[i] - c[i - 1]) / g->Ac[i];
          Ld -= bb; Lsub = bb;
        }
        Ld -= double(m * m) / (c[i] * c[i]);
        C.b[i] = 16.0 - w * Ld; C.a[i] = -w * Lsub; C.cc[i] = -w * Lsup;
      }
      C.factor();
    }
  }
  return true;
}

void FieldSolver::poisson(int m, const cd* rho, cd* psi) { Tp[m].solve(rho, psi); }

void FieldSolver::gradient(int m, const cd* psi, cd* Er, cd* Ep) const {
  const cd im(0.0, double(m));
  opG(psi, Er);
  for (int i = 0; i <= N; ++i) Er[i] = -Er[i];
  for (int i = 0; i < N; ++i) Ep[i] = -(im / g->c[i]) * psi[i];
}

void FieldSolver::opG(const cd* e, cd* out) const {  // ghost e_N = 0
  const auto& c = g->c;
  out[0] = 0.0;
  for (int i = 1; i <= N; ++i) {
    cd ei = (i < N) ? e[i] : cd(0.0);
    out[i] = (ei - e[i - 1]) / (c[i] - c[i - 1]);
  }
}
void FieldSolver::opg(const cd* b, cd* out) const {
  const auto& rho = g->rho;
  for (int i = 0; i < N; ++i) out[i] = (b[i + 1] - b[i]) / (rho[i + 1] - rho[i]);
}
void FieldSolver::opMf(const cd* b, cd* out, cd im) const {
  out[0] = 0.0;
  for (int i = 1; i <= N; ++i) out[i] = im / g->rho[i] * b[i];
}
void FieldSolver::opCE(const cd* Er, const cd* Ep, cd* out, cd im) const {
  const auto& c = g->c; const auto& rho = g->rho; const auto& Af = g->Af;
  out[0] = c[0] * Ep[0] / Af[0];
  for (int i = 1; i < N; ++i) out[i] = (c[i] * Ep[i] - c[i - 1] * Ep[i - 1]) / Af[i] - im / rho[i] * Er[i];
  out[N] = 0.0;
}
void FieldSolver::opCB(const cd* Bp, const cd* Br, cd* out, cd im) const {
  const auto& c = g->c; const auto& rho = g->rho; const auto& Ac = g->Ac;
  for (int i = 0; i < N; ++i) {
    cd inner = (i == 0) ? cd(0.0) : rho[i] * Bp[i];
    out[i] = (rho[i + 1] * Bp[i + 1] - inner) / Ac[i] - im / c[i] * Br[i];
  }
}

void FieldSolver::solveSlice(int m, int n, const cd* Rr0, const cd* Rp0, const cd* e0, const cd* b0, const cd* Rr1,
                             const cd* Lr1, const cd* Rp1, const cd* Lp1, const cd* e1, const cd* b1, const cd* jr,
                             const cd* jp, const cd* jz, const cd* Jr, const cd* Jp, cd* Rr, cd* Lr, cd* Rp,
                             cd* Lp, cd* e, cd* b) {
  const cd im(0.0, double(m));
  const auto& c = g->c;
  const double h = hs[n];
  const size_t T = size_t(cls[n]) * (M + 1) + m;
  for (int i = 0; i <= N; ++i) { Er1[i] = 0.5 * (Rr1[i] + Lr1[i]); Bp1[i] = 0.5 * (Rr1[i] - Lr1[i]); }
  for (int i = 0; i < N; ++i) { Ep1[i] = 0.5 * (Rp1[i] + Lp1[i]); Br1[i] = 0.5 * (Lp1[i] - Rp1[i]); }

  // known parts P (level s^0, trapezoid in s) and Q (from zeta_1)
  opG(e0, tf1.data()); opMf(b0, tf2.data(), im);
  for (int i = 0; i <= N; ++i) Pr[i] = Rr0[i] + 0.5 * D * (tf1[i] + tf2[i]) - 0.25 * D * Jr[i];
  opG(e1, tf1.data()); opMf(b1, tf2.data(), im);
  for (int i = 0; i <= N; ++i) Qr[i] = Lr1[i] - 0.25 * h * (tf1[i] - tf2[i]) - 0.5 * h * jr[i];
  opg(b0, tc1.data());
  for (int i = 0; i < N; ++i) Pp[i] = Rp0[i] + 0.5 * D * (im / c[i] * e0[i] - tc1[i]) - 0.25 * D * Jp[i];
  opg(b1, tc1.data());
  for (int i = 0; i < N; ++i) Qp[i] = Lp1[i] - 0.25 * h * (im / c[i] * e1[i] + tc1[i]) - 0.5 * h * jp[i];
  Pr[0] = 0.0; Qr[0] = 0.0;  // face 0 (axis): radial components not used / zero

  // B_z:  [16 - h(2D+h) Lap^f] b = 16 b1 - 8h [CE(zeta_1) + CE(bar)]
  for (int i = 0; i <= N; ++i) tf1[i] = 0.5 * (Pr[i] + Qr[i]);
  for (int i = 0; i < N; ++i) tc1[i] = 0.5 * (Pp[i] + Qp[i]);
  opCE(Er1.data(), Ep1.data(), rf.data(), im);
  opCE(tf1.data(), tc1.data(), tf2.data(), im);
  for (int i = 0; i <= N; ++i) rf[i] = 16.0 * b1[i] - 8.0 * h * (rf[i] + tf2[i]);
  std::fill(b, b + N + 1, cd(0.0));
  Tf[T].solve(rf.data(), b);

  // E_z:  [16 - h(2D+h) Lap^c] e = 16 e1 + 8h [CB(zeta_1) + CB(bar)] - 16 h j_z
  for (int i = 0; i <= N; ++i) tf1[i] = 0.5 * (Pr[i] - Qr[i]);
  for (int i = 0; i < N; ++i) tc1[i] = 0.5 * (Qp[i] - Pp[i]);
  opCB(Bp1.data(), Br1.data(), rc.data(), im);
  opCB(tf1.data(), tc1.data(), tc2.data(), im);
  for (int i = 0; i < N; ++i) rc[i] = 16.0 * e1[i] + 8.0 * h * (rc[i] + tc2[i]) - 16.0 * h * jz[i];
  Tc[T].solve(rc.data(), e);

  // explicit transverse update
  opG(e, tf1.data()); opMf(b, tf2.data(), im);
  for (int i = 0; i <= N; ++i) {
    Rr[i] = Pr[i] + 0.5 * D * (tf1[i] + tf2[i]);
    Lr[i] = Qr[i] - 0.25 * h * (tf1[i] - tf2[i]);
  }
  Rr[0] = 0.0; Lr[0] = 0.0;
  opg(b, tc1.data());
  for (int i = 0; i < N; ++i) {
    Rp[i] = Pp[i] + 0.5 * D * (im / c[i] * e[i] - tc1[i]);
    Lp[i] = Qp[i] - 0.25 * h * (im / c[i] * e[i] + tc1[i]);
  }
}

}  // namespace gem
