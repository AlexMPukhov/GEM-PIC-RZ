// ============================================================================
//  dielectric.cpp -- the slice equations with a dielectric eps(r) (block
//  'dielectric'), see FieldSolver::setDielectric in solver.hpp.
//
//  Per node of slice n (class h, Delta = D, tD = theta D), with A, B the parts
//  of R linear in E_z, B_z (faces: A = G e, B = (im/rho) b; centres:
//  A = (im/c) e, B = -g b) and the known parts Pk, Qk of R and L:
//     R = Pk + tD (A + B) - (2 tD chi / h) E        (J^1_n = Kn + (2 chi/h) E_n)
//     L = Qk - (h/4)(A - B) - (chi/2) E              (half-cell current -chi (E_{n+1} - E_n)/h)
//     E = (R + L)/2  =>  E = mu [(Pk + Qk)/2 + alpha A + beta B],
//     mu = 1 / (1 + tD chi/h + chi/4), alpha = tD/2 - h/8, beta = tD/2 + h/8.
//  The half-cell equations
//     16 b + 8h CE(E_n) = 16 b1 - 8h CE(E_{n+1}),
//     16 eps e - 8h CB(B_n) = 16 eps e1 + 8h CB(B_{n+1}) - 16 h j_z
//  are linear in (e, b). In vacuum (chi = 0) the identities CE(G e, (im/c) e) = 0,
//  CB((im/rho) b, g b) = 0 decouple them into the two tridiagonal systems of
//  solveSlice; with chi(r) the e and b of the modes m >= 1 couple where chi
//  changes. With (b_0, e_0, b_1, e_1, ..., b_N) interleaved the system is
//  pentadiagonal; its matrix is built by probing the same operator that gives
//  the right-hand side (5 columns at a time), so both agree by construction.
// ============================================================================
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

#include "solver.hpp"

namespace gem {

bool FieldSolver::setDielectric(const std::vector<double>& chiC1, const std::vector<double>& chiF1,
                                const std::vector<double>& chiC0, const std::vector<double>& chiF0) {
  if (!diel) {
    for (auto* v : {&dA, &dB, &dAp, &dBp, &dEr, &dEp, &dRr, &dLr, &dRp, &dLp, &dPkr, &dQkr, &dPkp, &dQkp, &dE, &dBb,
                    &dU, &dV})
      v->assign(N + 1, 0.0);
    dX.assign(2 * N + 1, 0.0);
    for (auto* v : {&dJ0r, &dJ1r, &dJ0rc, &dJ1rc}) v->assign(size_t(M + 1) * (N + 1), 0.0);
    for (auto* v : {&dJ0p, &dJ1p, &dJ0pc, &dJ1pc}) v->assign(size_t(M + 1) * N, 0.0);
  }
  const bool changed = !diel || chiC1 != dChiC1 || chiF1 != dChiF1;
  diel = true;
  dChiC1 = chiC1; dChiF1 = chiF1; dChiC0 = chiC0; dChiF0 = chiF0;
  if (changed) dielChanged = true;
  return changed;
}

void FieldSolver::dielBegin() {
  for (auto* v : {&dJ0r, &dJ1r, &dJ0p, &dJ1p}) std::fill(v->begin(), v->end(), cd(0.0));
}

void FieldSolver::dielNext() {
  dJ0r.swap(dJ0rc); dJ1r.swap(dJ1rc); dJ0p.swap(dJ0pc); dJ1p.swap(dJ1pc);
}

void FieldSolver::dielGet(std::vector<cd>& out) const {
  out.clear();
  for (const auto* v : {&dJ0r, &dJ0p, &dJ1r, &dJ1p}) out.insert(out.end(), v->begin(), v->end());
}

void FieldSolver::dielSet(const std::vector<cd>& in) {
  if (in.size() != dJ0r.size() + dJ0p.size() + dJ1r.size() + dJ1p.size())
    throw std::runtime_error("dielectric: polarization message of the wrong size");
  size_t o = 0;
  for (auto* v : {&dJ0r, &dJ0p, &dJ1r, &dJ1p}) {
    std::copy(in.begin() + long(o), in.begin() + long(o + v->size()), v->begin());
    o += v->size();
  }
}

void FieldSolver::dielLevel0(int m, int n, const cd* Rr0, const cd* Lr0, const cd* Rp0, const cd* Lp0, const cd* Rr0n,
                             const cd* Lr0n, const cd* Rp0n, const cd* Lp0n) {
  const double h = hs[n];
  const size_t of = size_t(m) * (N + 1), oc = size_t(m) * N;
  // J_n = -J_{n+1} - (2 chi/h)(E_{n+1} - E_n), E = (R + L)/2
  for (int i = 0; i <= N; ++i) {
    const double x = dChiF0[i];
    dJ0rc[of + i] = (x == 0.0 || i == 0) ? cd(0.0)
                                         : -dJ0r[of + i] - (x / h) * ((Rr0n[i] + Lr0n[i]) - (Rr0[i] + Lr0[i]));
  }
  for (int i = 0; i < N; ++i) {
    const double x = dChiC0[i];
    dJ0pc[oc + i] = (x == 0.0) ? cd(0.0) : -dJ0p[oc + i] - (x / h) * ((Rp0n[i] + Lp0n[i]) - (Rp0[i] + Lp0[i]));
  }
}

void FieldSolver::dielCorrection(int q, bool face, int i, double& fR, double& fL) const {
  const double h = hCls[q], tD = theta * D, x = face ? dChiF1[i] : dChiC1[i];
  const double mu = 1.0 / (1.0 + tD * x / h + 0.25 * x);
  fR = 1.0 - tD * x * mu / h;
  fL = -0.25 * x * mu;
}

void FieldSolver::dielTransverse(int m, double h, const cd* e, const cd* b, const cd* Pkr, const cd* Qkr,
                                 const cd* Pkp, const cd* Qkp, cd* Rr, cd* Lr, cd* Rp, cd* Lp, cd* Er, cd* Ep) {
  const cd im(0.0, double(m));
  const auto& c = g->c;
  const double tD = theta * D, al = 0.5 * tD - 0.125 * h, be = 0.5 * tD + 0.125 * h;
  opG(e, dA.data());
  opMf(b, dB.data(), im);
  for (int i = 0; i <= N; ++i) {
    const double x = dChiF1[i], mu = 1.0 / (1.0 + tD * x / h + 0.25 * x);
    const cd P = Pkr ? Pkr[i] : cd(0.0), Q = Qkr ? Qkr[i] : cd(0.0);
    const cd E = mu * (0.5 * (P + Q) + al * dA[i] + be * dB[i]);
    Er[i] = E;
    Rr[i] = P + tD * (dA[i] + dB[i]) - (2.0 * tD * x / h) * E;
    Lr[i] = Q - 0.25 * h * (dA[i] - dB[i]) - 0.5 * x * E;
  }
  Rr[0] = 0.0; Lr[0] = 0.0; Er[0] = 0.0;
  opg(b, dBp.data());
  for (int i = 0; i < N; ++i) {
    const double x = dChiC1[i], mu = 1.0 / (1.0 + tD * x / h + 0.25 * x);
    const cd A = im / c[i] * e[i], B = -dBp[i];
    const cd P = Pkp ? Pkp[i] : cd(0.0), Q = Qkp ? Qkp[i] : cd(0.0);
    const cd E = mu * (0.5 * (P + Q) + al * A + be * B);
    Ep[i] = E;
    Rp[i] = P + tD * (A + B) - (2.0 * tD * x / h) * E;
    Lp[i] = Q - 0.25 * h * (A - B) - 0.5 * x * E;
  }
}

void FieldSolver::dielOperator(int m, double h, const cd* e, const cd* b, cd* outE, cd* outB) {
  const cd im(0.0, double(m));
  dielTransverse(m, h, e, b, nullptr, nullptr, nullptr, nullptr, dRr.data(), dLr.data(), dRp.data(), dLp.data(),
                 dEr.data(), dEp.data());
  opCE(dEr.data(), dEp.data(), dE.data(), im);
  for (int i = 0; i <= N; ++i) outB[i] = 16.0 * b[i] + 8.0 * h * dE[i];
  for (int i = 0; i <= N; ++i) dBb[i] = 0.5 * (dRr[i] - dLr[i]);  // B_phi
  for (int i = 0; i < N; ++i) dV[i] = 0.5 * (dLp[i] - dRp[i]);    // B_r
  opCB(dBb.data(), dV.data(), dE.data(), im);
  for (int i = 0; i < N; ++i) outE[i] = 16.0 * (1.0 + dChiC1[i]) * e[i] - 8.0 * h * dE[i];
}

void FieldSolver::buildDiel() {
  const int n2 = 2 * N + 1;
  Bd.assign(size_t(nCls) * (M + 1), Banded5());
  std::vector<cd> e(N + 1), b(N + 1), oE(N + 1), oB(N + 1), x(n2);
  int mm = 0;
  auto inRow = [&](int r) { return (r % 2 == 0) ? rowInB(mm, r / 2) : rowInE((r - 1) / 2); };
  auto out = [&](int r) { return (r % 2 == 0) ? oB[r / 2] : oE[(r - 1) / 2]; };
  auto slot = [&](int j) -> cd& { return (j % 2 == 0) ? b[j / 2] : e[(j - 1) / 2]; };
  for (int q = 0; q < nCls; ++q) {
    const double h = hCls[q];
    for (int m = 0; m <= M; ++m) {
      mm = m;
      Banded5& B = Bd[size_t(q) * (M + 1) + m];
      B.n = n2;
      B.A.assign(size_t(n2) * 5, cd(0.0));
      // probe with the columns j = col (mod 5): a row couples to the columns r - 2 ... r + 2 only
      for (int col = 0; col < 5; ++col) {
        std::fill(e.begin(), e.end(), cd(0.0));
        std::fill(b.begin(), b.end(), cd(0.0));
        for (int j = col; j < n2; j += 5)
          if (inRow(j)) slot(j) = 1.0;
        dielOperator(m, h, e.data(), b.data(), oE.data(), oB.data());
        for (int r = 0; r < n2; ++r) {
          if (!inRow(r)) continue;
          for (int k = -2; k <= 2; ++k) {
            const int j = r + k;
            if (j >= 0 && j < n2 && j % 5 == col && inRow(j)) B.at(r, k) = out(r);
          }
        }
      }
      for (int r = 0; r < n2; ++r)
        if (!inRow(r)) B.at(r, 0) = 16.0;  // decoupled rows (outside the wall, B_z(rho_N) = 0, B_z(0) for m >= 1)
      // check the band against the operator (the probing assumes the bandwidth)
      std::fill(e.begin(), e.end(), cd(0.0));
      std::fill(b.begin(), b.end(), cd(0.0));
      for (int j = 0; j < n2; ++j) {
        x[j] = inRow(j) ? cd(std::sin(1.3 * j + 0.2 * m) + 0.1, std::cos(0.7 * j)) : cd(0.0);
        slot(j) = x[j];
      }
      dielOperator(m, h, e.data(), b.data(), oE.data(), oB.data());
      double err = 0.0, nrm = 0.0;
      for (int r = 0; r < n2; ++r) {
        if (!inRow(r)) continue;
        cd s = 0.0;
        for (int k = -2; k <= 2; ++k)
          if (r + k >= 0 && r + k < n2) s += B.at(r, k) * x[r + k];
        err = std::max(err, std::abs(s - out(r)));
        nrm = std::max(nrm, std::abs(out(r)));
      }
      if (err > 1e-10 * nrm)
        throw std::runtime_error("dielectric: the slice matrix is not pentadiagonal (internal error), mode " +
                                 std::to_string(m));
      B.factor();
    }
  }
}

void FieldSolver::solveSliceDiel(int m, int n, const cd* Rr0, const cd* Lr0, const cd* Rp0, const cd* e0,
                                 const cd* b0, const cd* Rr1, const cd* Lr1, const cd* Rp1, const cd* Lp1,
                                 const cd* e1, const cd* b1, const cd* jr, const cd* jp, const cd* jz, const cd* Jr,
                                 const cd* Jp, cd* Rr, cd* Lr, cd* Rp, cd* Lp, cd* e, cd* b) {
  const cd im(0.0, double(m));
  const auto& c = g->c;
  const double h = hs[n];
  const size_t T = size_t(cls[n]) * (M + 1) + m;
  const size_t of = size_t(m) * (N + 1), oc = size_t(m) * N;
  // known parts as in vacuum (level s^0 with weight 1 - theta, slice n+1 of s^1, currents)
  const double w0 = 1.0 - theta;
  opG(e0, tf1.data()); opMf(b0, tf2.data(), im);
  for (int i = 0; i <= N; ++i) Pr[i] = Rr0[i] + w0 * D * (tf1[i] + tf2[i]) - 0.25 * D * Jr[i];
  opG(e1, tf1.data()); opMf(b1, tf2.data(), im);
  for (int i = 0; i <= N; ++i) Qr[i] = Lr1[i] - 0.25 * h * (tf1[i] - tf2[i]) - 0.5 * h * jr[i];
  opg(b0, tc1.data());
  for (int i = 0; i < N; ++i) Pp[i] = Rp0[i] + w0 * D * (im / c[i] * e0[i] - tc1[i]) - 0.25 * D * Jp[i];
  opg(b1, tc1.data());
  for (int i = 0; i < N; ++i) Qp[i] = Lp1[i] - 0.25 * h * (im / c[i] * e1[i] + tc1[i]) - 0.5 * h * jp[i];
  // polarization: J^1_n = Kn + (2 chi/h) E_n, Kn = -J^1_{n+1} - (2 chi/h) E^1_{n+1}; R gets
  // -D [(1 - theta) J^0_n + theta J^1_n], L the half-cell current with -h/2: + (chi/2)(E_{n+1} - E_n)
  for (int i = 0; i <= N; ++i) {
    const double x = dChiF1[i];
    const cd E1 = 0.5 * (Rr1[i] + Lr1[i]);
    const cd Kn = -dJ1r[of + i] - (2.0 * x / h) * E1;
    dPkr[i] = Pr[i] - D * (w0 * dJ0rc[of + i] + theta * Kn);
    dQkr[i] = Qr[i] + 0.5 * x * E1;
    Er1[i] = E1;
    Bp1[i] = 0.5 * (Rr1[i] - Lr1[i]);
    tf2[i] = Kn;
  }
  dPkr[0] = 0.0; dQkr[0] = 0.0;
  for (int i = 0; i < N; ++i) {
    const double x = dChiC1[i];
    const cd E1 = 0.5 * (Rp1[i] + Lp1[i]);
    const cd Kn = -dJ1p[oc + i] - (2.0 * x / h) * E1;
    dPkp[i] = Pp[i] - D * (w0 * dJ0pc[oc + i] + theta * Kn);
    dQkp[i] = Qp[i] + 0.5 * x * E1;
    Ep1[i] = E1;
    Br1[i] = 0.5 * (Lp1[i] - Rp1[i]);
    tc2[i] = Kn;
  }
  // sloped wall: E_z at the wall (known)
  cd ew = 0.0;
  if (wNeu && wSlope != 0.0) ew = -wSlope * 0.5 * (Rr0[wK] + Lr0[wK]);
  std::fill(dU.begin(), dU.end(), cd(0.0));
  std::fill(dV.begin(), dV.end(), cd(0.0));
  if (wNeu) dU[wK] = ew;
  // transverse fields of the known parts (e = the wall value only, b = 0): R, L in rf, Qr (faces), Pp, Qp
  // (centres), E in Pr (faces), rc (centres)
  dielTransverse(m, h, dU.data(), dV.data(), dPkr.data(), dQkr.data(), dPkp.data(), dQkp.data(), rf.data(),
                 Qr.data(), Pp.data(), Qp.data(), Pr.data(), rc.data());
  // right-hand sides, interleaved (b_i at 2i, e_i at 2i+1); decoupled rows 0
  std::vector<cd>& x = dX;
  opCE(Er1.data(), Ep1.data(), tf1.data(), im);
  opCE(Pr.data(), rc.data(), dE.data(), im);
  for (int i = 0; i <= N; ++i) x[2 * i] = rowInB(m, i) ? 16.0 * b1[i] - 8.0 * h * (tf1[i] + dE[i]) : cd(0.0);
  for (int i = 0; i <= N; ++i) dA[i] = 0.5 * (rf[i] - Qr[i]);   // known B_phi
  for (int i = 0; i < N; ++i) dAp[i] = 0.5 * (Qp[i] - Pp[i]);   // known B_r
  opCB(Bp1.data(), Br1.data(), tc1.data(), im);
  opCB(dA.data(), dAp.data(), dE.data(), im);
  for (int i = 0; i < N; ++i)
    x[2 * i + 1] = rowInE(i) ? 16.0 * (1.0 + dChiC1[i]) * e1[i] + 8.0 * h * (tc1[i] + dE[i]) - 16.0 * h * jz[i]
                             : cd(0.0);
  Bd[T].solve(x.data());
  for (int i = 0; i <= N; ++i) b[i] = x[2 * i];
  for (int i = 0; i < N; ++i) e[i] = x[2 * i + 1];
  if (wNeu) e[wK] = ew;
  // transverse fields of slice n and its polarization current
  dielTransverse(m, h, e, b, dPkr.data(), dQkr.data(), dPkp.data(), dQkp.data(), Rr, Lr, Rp, Lp, dEr.data(),
                 dEp.data());
  for (int i = 0; i <= N; ++i)
    dJ1rc[of + i] = (i == 0 || i > wK) ? cd(0.0) : tf2[i] + (2.0 * dChiF1[i] / h) * dEr[i];
  for (int i = 0; i < N; ++i) dJ1pc[oc + i] = (i >= wK) ? cd(0.0) : tc2[i] + (2.0 * dChiC1[i] / h) * dEp[i];
  maskWall(Rr, Lr, Rp, Lp, e, b);
}

}  // namespace gem
