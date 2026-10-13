// ============================================================================
//  mapping.hpp -- free-hand grid: coordinate mapping x = X(xi) defined by
//  reference points (x_k, h_k) of the local grid step
//
//  The step h(x) is linear in x between the reference points and constant
//  outside them. The grid index xi is the "number of cells" up to x,
//      dx/dxi = h(x)   <=>   xi(x) = int_0^x dx' / h(x'),
//  so that node i lies at X(i) and the cell size near x is h(x). On a segment
//  with h = h_k + b (x - x_k) both directions are analytic:
//      xi(x) = xi_k + ln(h(x)/h_k) / b          (b != 0)
//      X(xi) = x_k + h_k (exp(b (xi - xi_k)) - 1) / b
//  The mapping is smooth within each segment; h is continuous, its slope may
//  jump at the reference points. |b| is the relative change of the step from
//  one cell to the next.
//
//  Used for the radial grid; written for any axis starting at x = 0 (the
//  longitudinal window can use it in the same way).
// ============================================================================
#pragma once
#include <string>
#include <vector>

namespace gem {

class StepMapping {
 public:
  // reference positions x (>= 0, strictly increasing) and steps h (> 0);
  // returns an empty string or an error message
  std::string set(const std::vector<double>& x, const std::vector<double>& h);

  double coord(double xi) const;  // X(xi), X(0) = 0
  double index(double x) const;   // xi(x) = inverse of coord
  double step(double x) const;    // h(x)

  double maxRelativeChange() const { return maxSlope; }  // max |dh/dx| = max cell-to-cell change
  double minStep() const;
  double maxStep() const;
  const std::vector<double>& points() const { return xk; }
  const std::vector<double>& steps() const { return hk; }

 private:
  // segments: [0, x0] constant h0, [x_k, x_{k+1}] linear, [x_last, inf) constant
  std::vector<double> xk, hk;   // reference points, starting with x = 0
  std::vector<double> xik;      // xi at the reference points
  std::vector<double> bk;       // slope of h on segment k
  double maxSlope = 0.0;
};

}  // namespace gem
