// ============================================================================
//  eigen.cpp -- symmetric tridiagonal eigenproblem (see eigen.hpp)
// ============================================================================
#include "eigen.hpp"

#include <cmath>

namespace gem {

bool symTridiagEigen(std::vector<double>& d, std::vector<double> e, std::vector<double>& Z) {
  const int n = int(d.size());
  Z.assign(size_t(n) * n, 0.0);
  for (int i = 0; i < n; ++i) Z[size_t(i) * n + i] = 1.0;
  e.resize(n, 0.0);
  e[n - 1] = 0.0;
  for (int l = 0; l < n; ++l) {
    int iter = 0, m;
    do {
      for (m = l; m < n - 1; ++m) {
        const double dd = std::fabs(d[m]) + std::fabs(d[m + 1]);
        if (std::fabs(e[m]) <= 1e-15 * dd) break;
      }
      if (m != l) {
        if (++iter > 60) return false;
        double g = (d[l + 1] - d[l]) / (2.0 * e[l]);
        double r = std::hypot(g, 1.0);
        g = d[m] - d[l] + e[l] / (g + (g >= 0 ? r : -r));
        double s = 1.0, c = 1.0, p = 0.0;
        int i;
        for (i = m - 1; i >= l; --i) {
          double f = s * e[i], b = c * e[i];
          r = std::hypot(f, g);
          e[i + 1] = r;
          if (r == 0.0) {
            d[i + 1] -= p;
            e[m] = 0.0;
            break;
          }
          s = f / r;
          c = g / r;
          g = d[i + 1] - p;
          r = (d[i] - g) * s + 2.0 * c * b;
          p = s * r;
          d[i + 1] = g + p;
          g = c * r - b;
          for (int k = 0; k < n; ++k) {  // rotate the eigenvector columns i, i+1
            f = Z[size_t(k) * n + i + 1];
            Z[size_t(k) * n + i + 1] = s * Z[size_t(k) * n + i] + c * f;
            Z[size_t(k) * n + i] = c * Z[size_t(k) * n + i] - s * f;
          }
        }
        if (r == 0.0 && i >= l) continue;
        d[l] -= p;
        e[l] = g;
        e[m] = 0.0;
      }
    } while (m != l);
  }
  return true;
}

}  // namespace gem
