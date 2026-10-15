// ============================================================================
//  smoothing.hpp -- radial smoothing of the plasma sources (plasma.smooth_length)
//
//  Regularizes the singular density and E_z spike where the plasma electrons
//  close the bubble (a caustic of the cold-plasma flow), as plasma.smooth_length
//  in QUARZ: every source of the streaming plasma is replaced by
//
//      j_perp <- (1 - a^2 Lap_vec)^-1 j_perp,     rho, j_z <- (1 - a^2 Lap_m)^-1 (rho, j_z),
//
//  per azimuthal mode m, which gives every ring a finite width a. Lap_vec is the
//  discrete vector Laplacian grad div - rot curl on the staggered grid (j_r on the
//  faces, j_phi on the centres) and Lap_m = div grad the scalar Laplacian of the
//  centres (the operator of the Poisson problem, zero beyond the wall). With the
//  operators of the field solver, div Lap_vec = Lap_m div holds exactly, so the
//  filter commutes with the divergence: smoothed currents conserve the smoothed
//  charge. The vector system is pentadiagonal (j_r and j_phi interleaved), the
//  scalar one tridiagonal; both are factorized once per mode.
// ============================================================================
#pragma once
#include <cmath>
#include <vector>

#include "fields.hpp"
#include "grid.hpp"

namespace gem {

class FieldSolver;

// banded matrix with two sub- and two super-diagonals, complex, LU without pivoting
struct Banded5 {
  int n = 0;
  std::vector<cd> A;  // row-major, 5 entries per row: offsets -2 .. +2 (LU factors after factor(), 1/U on the diagonal)
  cd& at(int i, int k) { return A[size_t(i) * 5 + size_t(k + 2)]; }
  cd at(int i, int k) const { return A[size_t(i) * 5 + size_t(k + 2)]; }
  void factor();
  void solve(cd* x) const;  // in place
};

class SourceFilter {
 public:
  // a: smoothing length (code units); a <= 0 switches the filter off
  void setup(const RadialGrid& grid, const FieldSolver& solver, int M, double a);
  bool on() const { return a2 > 0.0; }
  double length() const { return a2 > 0.0 ? std::sqrt(a2) : 0.0; }
  // j_r on the faces (N+1 values; face 0 untouched), j_phi on the centres (N values)
  void vector(int m, cd* jr, cd* jp);
  // centre quantity of mode m (j_z, or the m = 0 densities with m = 0)
  void scalar(int m, cd* f);
  void scalar0(double* f);  // real m = 0 centre quantity

  // factors for the parallel engine
  const std::vector<Banded5>& vectorFactors() const { return V; }
  // k = 0: sub-diagonal a, 1: cp, 2: 1/den (Thomas factors of the scalar filter), per mode
  const std::vector<std::vector<double>>& scalarFactors(int k) const { return k == 0 ? Sa : (k == 1 ? Scp : Sden); }
  int N = 0, M = 0;

 private:
  double a2 = 0.0;
  const RadialGrid* g = nullptr;
  std::vector<Banded5> V;            // per mode: 1 - a^2 Lap_vec
  std::vector<std::vector<double>> Sa, Sb, Sc;  // per mode: 1 - a^2 Lap_m (sub, diag, super), factorized
  std::vector<std::vector<double>> Scp, Sden;
  std::vector<cd> work, xv;
  void applyVecLap(int m, const cd* x, cd* y) const;  // y = Lap_vec x (interleaved storage)
};

}  // namespace gem
