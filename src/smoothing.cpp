// ============================================================================
//  smoothing.cpp -- radial smoothing of the plasma sources (see smoothing.hpp)
// ============================================================================
#include "smoothing.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "solver.hpp"

namespace gem {

void Banded5::factor() {
  // Doolittle LU in band storage: L has unit diagonal and two sub-diagonals, U the diagonal
  // and two super-diagonals (the matrix is diagonally dominant: no pivoting needed)
  for (int i = 0; i < n; ++i) {
    for (int k = -2; k <= -1; ++k) {  // L(i, i+k)
      const int j = i + k;
      if (j < 0) continue;
      cd s = at(i, k);
      for (int p = -2; p < k; ++p) {  // sum_q L(i,q) U(q,j), q = i+p < j
        const int q = i + p;
        if (q < 0 || j - q > 2) continue;
        s -= at(i, p) * at(q, j - q);
      }
      at(i, k) = s / at(j, 0);  // (U(j,j) not yet inverted: j < i)
    }
    for (int k = 0; k <= 2; ++k) {  // U(i, i+k)
      const int j = i + k;
      if (j >= n) continue;
      cd s = at(i, k);
      for (int p = -2; p <= -1; ++p) {
        const int q = i + p;
        if (q < 0 || j - q > 2) continue;
        s -= at(i, p) * at(q, j - q);
      }
      at(i, k) = s;
    }
  }
  for (int i = 0; i < n; ++i) at(i, 0) = 1.0 / at(i, 0);  // solve() multiplies by the reciprocal
}

void Banded5::solve(cd* x) const {
  for (int i = 0; i < n; ++i)
    for (int p = -2; p <= -1; ++p)
      if (i + p >= 0) x[i] -= at(i, p) * x[i + p];
  for (int i = n - 1; i >= 0; --i) {
    for (int k = 1; k <= 2; ++k)
      if (i + k < n) x[i] -= at(i, k) * x[i + k];
    x[i] *= at(i, 0);
  }
}

// Lap_vec x = grad div x - rot curl x with the operators of the field solver:
//   div_i  = (rho_{i+1} x_r,i+1 - rho_i x_r,i) / A^c_i + (i m / c_i) x_phi,i      (centres; 0 beyond the wall)
//   curl_i = (c_i x_phi,i - c_{i-1} x_phi,i-1) / A^f_i - (i m / rho_i) x_r,i      (faces; opCE)
//   face i:   (div_i - div_{i-1}) / (c_i - c_{i-1}) - (i m / rho_i) curl_i
//   centre i: (i m / c_i) div_i + (curl_{i+1} - curl_i) / (rho_{i+1} - rho_i)
// Storage: x[2i] = x_r at face i (i = 0..N), x[2i+1] = x_phi at centre i. Face 0 is not used;
// curl_0 = 0 for m >= 1 (regularity), which makes div Lap_vec = Lap_m div hold at the axis cell too.
void SourceFilter::applyVecLap(int m, const cd* x, cd* y) const {
  const auto& rho = g->rho;
  const auto& c = g->c;
  const cd im(0.0, double(m));
  std::vector<cd> dv(size_t(N) + 1, 0.0), cu(size_t(N) + 1, 0.0);
  for (int i = 0; i < N; ++i)
    dv[i] = (rho[i + 1] * x[2 * (i + 1)] - rho[i] * x[2 * i]) / g->Ac[i] + (im / c[i]) * x[2 * i + 1];
  dv[N] = 0.0;
  cu[0] = (m == 0) ? c[0] * x[1] / g->Af[0] : cd(0.0);  // a scalar of mode m >= 1 vanishes on the axis
  for (int i = 1; i < N; ++i)
    cu[i] = (c[i] * x[2 * i + 1] - c[i - 1] * x[2 * i - 1]) / g->Af[i] - (im / rho[i]) * x[2 * i];
  cu[N] = 0.0;
  y[0] = 0.0;
  for (int i = 1; i <= N; ++i) y[2 * i] = (dv[i] - dv[i - 1]) / (c[i] - c[i - 1]) - (im / rho[i]) * cu[i];
  for (int i = 0; i < N; ++i)
    y[2 * i + 1] = (im / c[i]) * dv[i] + (cu[i + 1] - cu[i]) / (rho[i + 1] - rho[i]);
}

void SourceFilter::setup(const RadialGrid& grid, const FieldSolver& solver, int M_, double a) {
  g = &grid; N = grid.N; M = M_;
  a2 = (a > 0.0) ? a * a : 0.0;
  V.clear(); Sa.clear(); Sb.clear(); Sc.clear(); Scp.clear(); Sden.clear();
  if (!on()) return;
  const int n = 2 * N + 1;
  work.assign(size_t(n), 0.0); xv.assign(size_t(n), 0.0);
  std::vector<cd> e(size_t(n), 0.0), col(size_t(n), 0.0);
  V.assign(size_t(M) + 1, Banded5());
  for (int m = 0; m <= M; ++m) {
    Banded5& B = V[m];
    B.n = n;
    B.A.assign(size_t(n) * 5, 0.0);
    // columns of 1 - a^2 Lap_vec by probing with unit vectors (the operator couples
    // only neighbours: rows k-2 .. k+2 of column k)
    for (int k = 0; k < n; ++k) {
      e[k] = 1.0;
      applyVecLap(m, e.data(), col.data());
      e[k] = 0.0;
      for (int i = std::max(0, k - 2); i <= std::min(n - 1, k + 2); ++i) {
        cd v = -a2 * col[i];
        if (i == k) v += 1.0;
        if (k == 0 || i == 0) v = (i == k) ? cd(1.0) : cd(0.0);  // face 0: pass-through
        B.at(i, k - i) = v;
      }
      for (int i = 0; i < n; ++i)
        if ((i < k - 2 || i > k + 2) && std::abs(col[i]) > 0.0) throw std::logic_error("smoothing: bandwidth > 2");
    }
    B.factor();
  }
  // scalar filter: 1 - a^2 Lap_m = 1 + a^2 (-Lap_m), the Poisson matrix of the field solver
  Sa.resize(M + 1); Sb.resize(M + 1); Sc.resize(M + 1); Scp.resize(M + 1); Sden.resize(M + 1);
  for (int m = 0; m <= M; ++m) {
    const Tridiag& P = solver.poissonMatrix(m);
    Sa[m].resize(N); Sb[m].resize(N); Sc[m].resize(N); Scp[m].resize(N); Sden[m].resize(N);
    for (int i = 0; i < N; ++i) { Sa[m][i] = a2 * P.a[i]; Sb[m][i] = 1.0 + a2 * P.b[i]; Sc[m][i] = a2 * P.cc[i]; }
    double d = Sb[m][0];
    Sden[m][0] = 1.0 / d; Scp[m][0] = Sc[m][0] / d;
    for (int i = 1; i < N; ++i) {
      d = Sb[m][i] - Sa[m][i] * Scp[m][i - 1];
      Sden[m][i] = 1.0 / d; Scp[m][i] = Sc[m][i] / d;
    }
  }
}

void SourceFilter::vector(int m, cd* jr, cd* jp) {
  if (!on()) return;
  for (int i = 0; i <= N; ++i) xv[2 * i] = jr[i];
  for (int i = 0; i < N; ++i) xv[2 * i + 1] = jp[i];
  V[m].solve(xv.data());
  for (int i = 1; i <= N; ++i) jr[i] = xv[2 * i];
  for (int i = 0; i < N; ++i) jp[i] = xv[2 * i + 1];
}

void SourceFilter::scalar(int m, cd* f) {
  if (!on()) return;
  const auto &a = Sa[m], &cp = Scp[m], &den = Sden[m];
  work[0] = f[0] * den[0];
  for (int i = 1; i < N; ++i) work[i] = (f[i] - a[i] * work[i - 1]) * den[i];
  f[N - 1] = work[N - 1];
  for (int i = N - 2; i >= 0; --i) f[i] = work[i] - cp[i] * f[i + 1];
}

void SourceFilter::scalar0(double* f) {
  if (!on()) return;
  for (int i = 0; i < N; ++i) xv[i] = f[i];
  scalar(0, xv.data());
  for (int i = 0; i < N; ++i) f[i] = xv[i].real();
}

}  // namespace gem
