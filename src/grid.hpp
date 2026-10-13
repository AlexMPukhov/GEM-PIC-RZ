// ============================================================================
//  grid.hpp -- non-uniform staggered radial grid (Section 5.1 of the paper)
//
//    faces   rho_i = R(i),       i = 0..N      (rho_0 = 0 is the axis)
//    centres c_i   = R(i + 1/2), i = 0..N-1    (ghost c_N = 2 rho_N - c_{N-1})
//
//    faces   : B_z, E_r, B_phi, R_r, L_r, j_r
//    centres : E_z, E_phi, B_r, R_phi, L_phi, j_phi, j_z
// ============================================================================
#pragma once
#include <vector>

#include "config.hpp"

namespace gem {

struct RadialGrid {
  int N = 0;
  std::vector<double> rho;     // faces rho_0..rho_N
  std::vector<double> c;       // centres c_0..c_{N-1} and ghost c_N
  std::vector<double> Af, Ac;  // field-solver measures, Eq. (measures)
  std::vector<double> Vf, Vc;  // deposition volumes int S_i r dr
  std::vector<double> kC, kF;  // shape parameters of the deposition intervals (see grid.cpp)
  double maxVolumeMismatch = 0; // max |V/A - 1| over the nodes used by the solver

  void build(const GridConfig& g);

  // interpolation weights: point lies between index a and a+1
  struct W {
    int a;
    double wa, wb;
    int mode;  // 0: between a and a+1;  1: inside r < c_0 (centres);  2: outside
  };
  // deposition: area weighting with a quadratic correction, chosen so that
  // uniform density deposits exactly the solver measures A^c, A^f
  W locC(double r) const;
  W locF(double r) const;
  // field gather: linear in r (exact for fields ~ r near the axis)
  W locC_lin(double r) const;
  W locF_lin(double r) const;
};

}  // namespace gem
