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
  for (auto* v : {&e0t, &b0t, &et, &bt, &sx}) v->assign(N + 1, 0.0);
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
    const double w = h * (4.0 * theta * D + h);  // h (2 Delta + h) for the trapezoid
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
  if (stiff()) {
    // the radial Laplacians (rows of Tf, Tc): Lap = (b - 16) / (-w) of the matrices above with w = 1
    Lc.assign(M + 1, Tridiag());
    Lf.assign(M + 1, Tridiag());
    for (int m = 0; m <= M; ++m) {
      Tridiag &F = Lf[m], &C = Lc[m];
      F.i0 = (m == 0) ? 0 : 1; F.n = N - F.i0;
      F.a.assign(F.n, 0.0); F.b.assign(F.n, 0.0); F.cc.assign(F.n, 0.0);
      for (int r = 0; r < F.n; ++r) {
        const int i = F.i0 + r;
        if (i == 0) {
          const double a = c[0] / rho[1] / g->Af[0];
          F.b[r] = -a; F.cc[r] = a;
        } else {
          const double a = c[i] / (rho[i + 1] - rho[i]) / g->Af[i];
          const double bb = c[i - 1] / (rho[i] - rho[i - 1]) / g->Af[i];
          F.b[r] = -a - bb - double(m * m) / (rho[i] * rho[i]);
          F.cc[r] = (i + 1 <= N - 1) ? a : 0.0;
          F.a[r] = (i - 1 >= F.i0) ? bb : 0.0;
        }
      }
      C.i0 = 0; C.n = N;
      C.a.assign(N, 0.0); C.b.assign(N, 0.0); C.cc.assign(N, 0.0);
      for (int i = 0; i < N; ++i) {
        const double a = rho[i + 1] / (c[i + 1] - c[i]) / g->Ac[i];
        C.b[i] = -a - double(m * m) / (c[i] * c[i]);
        C.cc[i] = (i + 1 < N) ? a : 0.0;
        if (i > 0) {
          const double bb = rho[i] / (c[i] - c[i - 1]) / g->Ac[i];
          C.b[i] -= bb; C.a[i] = bb;
        }
      }
    }
    // per class and mode: S = 1 - a Lap and the implicit matrix 16 - al1 Lap + al2 Lap^2,
    // al1 = 16 a + h^2 + 4 h D theta, al2 = a (h^2 + 4 h D)
    St.assign(size_t(nCls) * (M + 1), Stiff());
    for (int q = 0; q < nCls; ++q) {
      const double h = hCls[q], a = thetaStiff * h * D;
      const double al1 = 16.0 * a + h * h + 4.0 * h * D * theta, al2 = a * (h * h + 4.0 * h * D);
      for (int m = 0; m <= M; ++m) {
        Stiff& S = St[size_t(q) * (M + 1) + m];
        S.a = a;
        for (int k = 0; k < 2; ++k) {
          const Tridiag& L = k ? Lf[m] : Lc[m];
          Tridiag& T = k ? S.Sf : S.Sc;
          Banded5& B = k ? S.Bf : S.Bc;
          const int n = L.n;
          T.i0 = L.i0; T.n = n;
          T.a.resize(n); T.b.resize(n); T.cc.resize(n);
          for (int r = 0; r < n; ++r) { T.a[r] = -a * L.a[r]; T.b[r] = 1.0 - a * L.b[r]; T.cc[r] = -a * L.cc[r]; }
          T.factor();
          B.n = n;
          B.A.assign(size_t(n) * 5, cd(0.0));
          auto La = [&](int r) { return (r >= 0 && r < n) ? L.a[r] : 0.0; };
          auto Lb = [&](int r) { return (r >= 0 && r < n) ? L.b[r] : 0.0; };
          auto Lu = [&](int r) { return (r >= 0 && r < n) ? L.cc[r] : 0.0; };
          for (int r = 0; r < n; ++r) {  // rows of Lap^2
            const double q2[5] = {La(r) * La(r - 1), La(r) * Lb(r - 1) + Lb(r) * La(r),
                                  La(r) * Lu(r - 1) + Lb(r) * Lb(r) + Lu(r) * La(r + 1), Lb(r) * Lu(r) + Lu(r) * Lb(r + 1),
                                  Lu(r) * Lu(r + 1)};
            const double l1[5] = {0.0, La(r), Lb(r), Lu(r), 0.0};
            for (int k5 = 0; k5 < 5; ++k5) {
              if (r + k5 - 2 < 0 || r + k5 - 2 >= n) continue;
              B.at(r, k5 - 2) = (k5 == 2 ? 16.0 : 0.0) - al1 * l1[k5] + al2 * q2[k5];
            }
          }
          B.factor();
        }
      }
    }
  }
  return true;
}

void FieldSolver::lapRows(const Tridiag& L, const cd* x, cd* y) const {
  for (int r = 0; r < L.n; ++r) {
    const int i = L.i0 + r;
    cd s = L.b[r] * x[i];
    if (r > 0) s += L.a[r] * x[i - 1];
    if (r + 1 < L.n) s += L.cc[r] * x[i + 1];
    y[i] = s;
  }
}

void FieldSolver::thetaApply(Stiff& S, bool faces, int m, bool one, const cd* x, cd* y) {
  const Tridiag& L = faces ? Lf[m] : Lc[m];
  Tridiag& T = faces ? S.Sf : S.Sc;
  const int len = faces ? N + 1 : N;
  std::fill(sx.begin(), sx.end(), cd(0.0));
  if (one) {
    for (int r = 0; r < L.n; ++r) sx[L.i0 + r] = (1.0 - theta) * x[L.i0 + r];
  } else {
    lapRows(L, x, sx.data());
    for (int r = 0; r < L.n; ++r) sx[L.i0 + r] = theta * x[L.i0 + r] - S.a * sx[L.i0 + r];
  }
  std::fill(y, y + len, cd(0.0));
  T.solve(sx.data(), y);
}

void FieldSolver::stiffSolve(Stiff& S, bool faces, int m, const cd* rhs, cd* x) {
  const Tridiag& L = faces ? Lf[m] : Lc[m];
  const Banded5& B = faces ? S.Bf : S.Bc;
  lapRows(L, rhs, sx.data());  // right-hand side times S = 1 - a Lap
  for (int r = 0; r < L.n; ++r) x[L.i0 + r] = rhs[L.i0 + r] - S.a * sx[L.i0 + r];
  B.solve(x + L.i0);
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
  Stiff* S = stiff() ? &St[T] : nullptr;
  for (int i = 0; i <= N; ++i) { Er1[i] = 0.5 * (Rr1[i] + Lr1[i]); Bp1[i] = 0.5 * (Rr1[i] - Lr1[i]); }
  for (int i = 0; i < N; ++i) { Ep1[i] = 0.5 * (Rp1[i] + Lp1[i]); Br1[i] = 0.5 * (Lp1[i] - Rp1[i]); }

  // known parts P (level s^0, trapezoid in s) and Q (from zeta_1); thetaStiff: the fields of s^0 as (1 - Theta) e0,
  // (1 - Theta) b0
  double w0 = 1.0 - theta;
  if (S) {
    thetaApply(*S, false, m, true, e0, e0t.data());
    thetaApply(*S, true, m, true, b0, b0t.data());
    e0 = e0t.data(); b0 = b0t.data();
    w0 = 1.0;
  }
  opG(e0, tf1.data()); opMf(b0, tf2.data(), im);
  for (int i = 0; i <= N; ++i) Pr[i] = Rr0[i] + w0 * D * (tf1[i] + tf2[i]) - 0.25 * D * Jr[i];
  opG(e1, tf1.data()); opMf(b1, tf2.data(), im);
  for (int i = 0; i <= N; ++i) Qr[i] = Lr1[i] - 0.25 * h * (tf1[i] - tf2[i]) - 0.5 * h * jr[i];
  opg(b0, tc1.data());
  for (int i = 0; i < N; ++i) Pp[i] = Rp0[i] + w0 * D * (im / c[i] * e0[i] - tc1[i]) - 0.25 * D * Jp[i];
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
  if (S) stiffSolve(*S, true, m, rf.data(), b);
  else Tf[T].solve(rf.data(), b);

  // E_z:  [16 - h(2D+h) Lap^c] e = 16 e1 + 8h [CB(zeta_1) + CB(bar)] - 16 h j_z
  for (int i = 0; i <= N; ++i) tf1[i] = 0.5 * (Pr[i] - Qr[i]);
  for (int i = 0; i < N; ++i) tc1[i] = 0.5 * (Qp[i] - Pp[i]);
  opCB(Bp1.data(), Br1.data(), rc.data(), im);
  opCB(tf1.data(), tc1.data(), tc2.data(), im);
  for (int i = 0; i < N; ++i) rc[i] = 16.0 * e1[i] + 8.0 * h * (rc[i] + tc2[i]) - 16.0 * h * jz[i];
  if (S) stiffSolve(*S, false, m, rc.data(), e);
  else Tc[T].solve(rc.data(), e);

  // explicit transverse update (thetaStiff: R with Theta e, Theta b)
  const cd *eR = e, *bR = b;
  double w1 = theta;
  if (S) {
    thetaApply(*S, false, m, false, e, et.data());
    thetaApply(*S, true, m, false, b, bt.data());
    eR = et.data(); bR = bt.data();
    w1 = 1.0;
  }
  opG(eR, tf1.data()); opMf(bR, tf2.data(), im);
  for (int i = 0; i <= N; ++i) Rr[i] = Pr[i] + w1 * D * (tf1[i] + tf2[i]);
  opG(e, tf1.data()); opMf(b, tf2.data(), im);
  for (int i = 0; i <= N; ++i) Lr[i] = Qr[i] - 0.25 * h * (tf1[i] - tf2[i]);
  Rr[0] = 0.0; Lr[0] = 0.0;
  opg(bR, tc1.data());
  for (int i = 0; i < N; ++i) Rp[i] = Pp[i] + w1 * D * (im / c[i] * eR[i] - tc1[i]);
  opg(b, tc1.data());
  for (int i = 0; i < N; ++i) Lp[i] = Qp[i] - 0.25 * h * (im / c[i] * e[i] + tc1[i]);
}

}  // namespace gem
