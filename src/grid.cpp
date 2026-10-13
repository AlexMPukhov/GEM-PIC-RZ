// ============================================================================
//  grid.cpp -- radial grid construction and interpolation weights
// ============================================================================
#include "grid.hpp"

#include <algorithm>
#include <cmath>

namespace gem {

void RadialGrid::build(const GridConfig& g) {
  N = g.Nr;
  rho.resize(N + 1);
  c.resize(N + 1);
  for (int i = 0; i <= N; ++i) rho[i] = g.mapping(i);
  for (int i = 0; i < N; ++i) c[i] = g.mapping(i + 0.5);
  c[N] = 2.0 * rho[N] - c[N - 1];

  Af.assign(N + 1, 0.0);
  Af[0] = 0.5 * c[0] * c[0];  // axis: disk r < c_0
  for (int i = 1; i <= N; ++i) Af[i] = rho[i] * (c[i] - c[i - 1]);
  Ac.resize(N);
  for (int i = 0; i < N; ++i) Ac[i] = c[i] * (rho[i + 1] - rho[i]);

  // ---------------------------------------------------------------------------
  // Deposition shape. In the interval [a, b] between two neighbouring nodes a
  // particle at r gives the weight
  //     w_b = x + kappa x (1 - x),   w_a = 1 - w_b,   x = (r^2 - a^2)/(b^2 - a^2),
  // i.e. area weighting (kappa = 0) plus a quadratic correction (|kappa| <= 1
  // keeps the weights monotonic and in [0, 1]). A uniform density puts the
  // share alpha = 1/2 - kappa/6 of the interval area on the left node. kappa is
  // chosen node by node, from the axis outwards, so that the deposited volume
  // of every node equals the measure of the field solver (A^c or A^f). Then
  //  * a uniform density is deposited exactly, and
  //  * the solver sees exactly the deposited charge (discrete Gauss law).
  // Plain area weighting violates the second point on a stretched grid by
  // O(R''/r) (e.g. 1.4 % in the enclosed charge of a bunch in the transition).
  // ---------------------------------------------------------------------------
  auto build = [](const std::vector<double>& node, int nNodes, double V0inner, const std::vector<double>& meas,
                  int firstMatched, double alpha0, std::vector<double>& kap, std::vector<double>& V) {
    const int nI = nNodes - 1;  // intervals
    kap.assign(nI, 0.0);
    V.assign(nNodes, 0.0);
    std::vector<double> A(nI);
    for (int j = 0; j < nI; ++j) A[j] = 0.5 * (node[j + 1] * node[j + 1] - node[j] * node[j]);
    double prevRight = V0inner;  // volume that node j receives from the interval on its left
    for (int j = 0; j < nI; ++j) {
      double alpha = (j < firstMatched) ? alpha0 : (meas[j] - prevRight) / A[j];
      double k = 6.0 * (0.5 - alpha);
      k = std::max(-1.0, std::min(1.0, k));
      alpha = 0.5 - k / 6.0;
      kap[j] = k;
      V[j] = prevRight + alpha * A[j];
      prevRight = (1.0 - alpha) * A[j];
    }
    V[nNodes - 1] = prevRight;
  };
  // centres: node 0 also collects everything inside r < c_0; the last centre
  // also collects c_{N-1} < r < rho_N
  build(c, N, 0.5 * c[0] * c[0], Ac, 0, 0.5, kC, Vc);
  Vc[N - 1] += 0.5 * (rho[N] * rho[N] - c[N - 1] * c[N - 1]);
  // faces: face 0 (axis) is not used by the solver; alpha = 1/2 in the first interval
  std::vector<double> AfFaces(Af.begin(), Af.end());
  build(rho, N + 1, 0.0, AfFaces, 1, 0.5, kF, Vf);

  maxVolumeMismatch = 0.0;
  for (int i = 0; i < N - 1; ++i) maxVolumeMismatch = std::max(maxVolumeMismatch, std::fabs(Vc[i] / Ac[i] - 1.0));
  for (int i = 1; i < N; ++i) maxVolumeMismatch = std::max(maxVolumeMismatch, std::fabs(Vf[i] / Af[i] - 1.0));
}

RadialGrid::W RadialGrid::locC(double r) const {
  if (r < c[0]) return {0, 1.0, 0.0, 1};
  if (r >= c[N - 1]) return {N - 1, 1.0, 0.0, 2};
  int a = int(std::upper_bound(c.begin(), c.begin() + N, r) - c.begin()) - 1;
  const double x = (r * r - c[a] * c[a]) / (c[a + 1] * c[a + 1] - c[a] * c[a]);
  const double wb = x + kC[a] * x * (1.0 - x);
  return {a, 1.0 - wb, wb, 0};
}

RadialGrid::W RadialGrid::locF(double r) const {
  if (r >= rho[N]) return {N, 1.0, 0.0, 2};
  int a = int(std::upper_bound(rho.begin(), rho.end(), r) - rho.begin()) - 1;
  const double x = (r * r - rho[a] * rho[a]) / (rho[a + 1] * rho[a + 1] - rho[a] * rho[a]);
  const double wb = x + kF[a] * x * (1.0 - x);
  return {a, 1.0 - wb, wb, 0};
}

RadialGrid::W RadialGrid::locC_lin(double r) const {
  if (r < c[0]) return {0, 1.0, 0.0, 1};
  if (r >= c[N - 1]) return {N - 1, 1.0, 0.0, 2};
  int a = int(std::upper_bound(c.begin(), c.begin() + N, r) - c.begin()) - 1;
  double wb = (r - c[a]) / (c[a + 1] - c[a]);
  return {a, 1.0 - wb, wb, 0};
}

RadialGrid::W RadialGrid::locF_lin(double r) const {
  if (r >= rho[N]) return {N, 1.0, 0.0, 2};
  int a = int(std::upper_bound(rho.begin(), rho.end(), r) - rho.begin()) - 1;
  double wb = (r - rho[a]) / (rho[a + 1] - rho[a]);
  return {a, 1.0 - wb, wb, 0};
}

}  // namespace gem
