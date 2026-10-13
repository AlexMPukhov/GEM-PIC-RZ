// ============================================================================
//  mapping.cpp -- free-hand grid mapping (see mapping.hpp)
// ============================================================================
#include "mapping.hpp"

#include <algorithm>
#include <cmath>

namespace gem {

namespace {
// ln(1 + y) / y and (exp(y) - 1) / y, accurate for small y
double log1pOver(double y) { return std::fabs(y) < 1e-8 ? 1.0 - 0.5 * y : std::log1p(y) / y; }
double expm1Over(double y) { return std::fabs(y) < 1e-8 ? 1.0 + 0.5 * y : std::expm1(y) / y; }
}  // namespace

std::string StepMapping::set(const std::vector<double>& x, const std::vector<double>& h) {
  if (x.empty() || x.size() != h.size()) return "needs as many steps as reference points (at least one)";
  for (size_t k = 0; k < x.size(); ++k) {
    if (!(x[k] >= 0)) return "reference points must not be negative";
    if (k > 0 && !(x[k] > x[k - 1])) return "reference points must be strictly increasing";
    if (!(h[k] > 0)) return "grid steps must be positive";
  }
  // prepend x = 0 with the first step (constant up to the first point)
  xk.clear();
  hk.clear();
  if (x[0] > 0) { xk.push_back(0.0); hk.push_back(h[0]); }
  xk.insert(xk.end(), x.begin(), x.end());
  hk.insert(hk.end(), h.begin(), h.end());
  const size_t K = xk.size();
  bk.assign(K, 0.0);
  xik.assign(K, 0.0);
  maxSlope = 0.0;
  for (size_t k = 0; k + 1 < K; ++k) {
    const double L = xk[k + 1] - xk[k];
    bk[k] = (hk[k + 1] - hk[k]) / L;
    maxSlope = std::max(maxSlope, std::fabs(bk[k]));
    // xi increment over the segment: int dx / (h_k + b (x - x_k)) = ln(h_{k+1}/h_k) / b
    xik[k + 1] = xik[k] + (L / hk[k]) * log1pOver(bk[k] * L / hk[k]);
  }
  bk[K - 1] = 0.0;  // constant beyond the last point
  return "";
}

double StepMapping::step(double x) const {
  if (x <= xk.front()) return hk.front();
  if (x >= xk.back()) return hk.back();
  const size_t k = size_t(std::upper_bound(xk.begin(), xk.end(), x) - xk.begin()) - 1;
  return hk[k] + bk[k] * (x - xk[k]);
}

double StepMapping::index(double x) const {
  if (x <= 0) return x / hk.front();
  const size_t k = std::min(size_t(std::upper_bound(xk.begin(), xk.end(), x) - xk.begin()) - 1, xk.size() - 1);
  const double d = x - xk[k];
  return xik[k] + (d / hk[k]) * log1pOver(bk[k] * d / hk[k]);
}

double StepMapping::coord(double xi) const {
  if (xi <= 0) return xi * hk.front();
  const size_t k = std::min(size_t(std::upper_bound(xik.begin(), xik.end(), xi) - xik.begin()) - 1, xik.size() - 1);
  const double t = xi - xik[k];
  return xk[k] + hk[k] * t * expm1Over(bk[k] * t);
}

double StepMapping::minStep() const { return *std::min_element(hk.begin(), hk.end()); }
double StepMapping::maxStep() const { return *std::max_element(hk.begin(), hk.end()); }

}  // namespace gem
