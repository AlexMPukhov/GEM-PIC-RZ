// ============================================================================
//  bound.hpp -- the bound-electron medium of a run (block 'bound_electrons'):
//  Bloch nodes on the cell centres x n_phi angles, their field E_x from the
//  modes, their current projected back to the modes (bloch.hpp: physics of
//  one node).
//
//  The medium current enters the field solver as an ordinary current:
//    slice current J_n of the trapezoidal rule over the half-cells, (J_n + J_{n+1})/2 = n_a (P_n - P_{n+1})/h
//    (as for the dielectric; the local d_t P of the discrete state would be the difference of two large terms):
//      R update -D [(1 - theta) J^0_n + theta J^1_n]   (in Jr, Jp: +4 [...]),
//    half-cell current (n_a P_n - n_a P_{n+1}) / h: L update -h/2 j (in jr, jp).
//  It has no z component (scalar dipole along x). Level s^1 is advanced with
//  the fields of the current pass and committed with the final fields of the
//  slice (R, L of slice n are then corrected for the change); level s^0 is
//  recomputed from the stored slices (no state is kept between the steps:
//  the atoms of the plane z = s^1 are fresh ahead of the window).
// ============================================================================
#pragma once
#include <vector>

#include "bloch.hpp"
#include "config.hpp"
#include "fields.hpp"
#include "grid.hpp"
#include "ionization.hpp"

namespace gem {

class BoundMedium {
 public:
  bool on = false;
  BlochLines L;
  int nphi = 0, N = 0, M = 0, nc = 0;  // nc: centres with atoms (c_i < r_out, inside the wall)
  // T: longest time in the window (ADK thresholds)
  void setup(const BoundElectronsConfig& c, const RadialGrid& g, int M, const Units& u, double T);
  // start of a march: atom densities of the levels s^0, s^1 per centre; states ground (ahead of the window); K: wall
  void begin(const std::vector<double>& na0, const std::vector<double>& na1, int K);
  bool ionizing() const { return ion; }
  // normalized populations (ground, excited 1, 2) of level s^1 at the node nearest to (x, y), at the slice before the
  // half-cell being pushed (the particles of the linked ion species ionize with sum_k p_k w_k)
  void populations(double x, double y, double& p0, double& p1, double& p2) const;
  // ADK rates of the levels [1/omega0] in the field E [m c omega0 / e]
  double rateLevel(int k, double E) const { return adkL[k].rate(0, E); }
  // level s^0: advance from slice n+1 to n with the stored fields; J^0 of slice n (modes) in J0r, J0p
  void level0(int n, double h, const Fields& fold);
  // level s^1: slice n from n+1 with the fields of fnew (slice n: current pass); J^1 of slice n and the half-cell
  // current in J1r, J1p, Hr, Hp. commit: the state of slice n becomes the running state (the slice is final)
  void level1(int n, double h, const Fields& fnew, bool commit);
  // modes: faces (M+1)(N+1), centres (M+1) N
  std::vector<cd> J0r, J0p, J1r, J1p, Hr, Hp;
  // running states at slice n+1 for the MPI messages (both levels; nodes x 9 doubles each + P of level 1)
  void getState(std::vector<double>& out) const;
  void setState(const std::vector<double>& in);
  // diagnostics of the last committed slice of level s^1 (m = 0 average over the angles): excited population
  double excitedMax() const { return exMax; }
  // excited population on the axis at the rear slice of the last march (centre 0, mean over the angles)
  double excitedAxis() const { return exAxis; }
  // neutral fraction (trace of the density matrix) on the axis and its area mean over the medium, rear slice
  double neutralAxis() const { return nuAxis; }
  double neutralMean() const { return nuMean; }

 private:
  const RadialGrid* g = nullptr;
  BoundElectronsConfig cfg;
  std::vector<double> cs, sn;          // cos, sin of the angles
  std::vector<double> na0, na1;        // atom density of the levels per centre
  std::vector<BlochNode> st0, st1;     // running states (slice n+1)
  std::vector<BlochNode> tmp1;         // tentative states of level 1 (slice n)
  std::vector<double> Ex0, Ex1, Exa;   // E_x at the nodes: slice n+1 of both levels, slice n
  std::vector<double> Jx, Hx;          // node currents
  std::vector<double> Jr0, Jr1;        // running slice currents of both levels (slice n+1)
  double exMax = 0.0, exAxis = 0.0, nuAxis = 1.0, nuMean = 1.0;
  bool ion = false;
  AdkRates adkL[3];  // ground, excited 1, 2
  void fieldAt(const Fields& F, int n, std::vector<double>& Ex) const;
  void rates(double E0, double E1, double& g0, double& g1, double& g2) const;  // mean of the half-cell ends
  void project(const std::vector<double>& x, cd* fr, cd* fp) const;  // node values of J_x -> modes (r faces, phi centres)
};

}  // namespace gem
