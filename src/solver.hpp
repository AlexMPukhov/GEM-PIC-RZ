// ============================================================================
//  solver.hpp -- field solver for one zeta-slice and one azimuthal mode
//                (finite differences in r, Section 5.3 of the paper)
// ============================================================================
#pragma once
#include <vector>

#include "fields.hpp"
#include "grid.hpp"
#include "zgrid.hpp"

namespace gem {

struct EngineState;

// Tridiagonal system with real coefficients and complex right-hand side,
// rows i0..i0+n-1 of the full index range (other entries stay untouched).
struct Tridiag {
  int i0 = 0, n = 0;
  std::vector<double> a, b, cc;  // sub, diag, super
  std::vector<double> cp, den;   // Thomas factors
  std::vector<cd> dp;            // scratch
  void factor();
  void solve(const cd* rhs, cd* x);
};

class FieldSolver {
 public:
  // radial operators and the Poisson matrices
  void setup(const RadialGrid& grid, int M);
  // factorize the matrices [16 - h_n(2D+h_n) Lap_m] of all modes for the step
  // D and the half-cells h_n of the window (one set per distinct h_n; nothing
  // to do if D and the window are the ones of the last call). Returns true if
  // the matrices were (re)computed.
  bool prepare(double D, const ZetaGrid& zg);
  int classes() const { return nCls; }
  int version() const { return ver; }  // increases whenever the matrices change
  int classOf(int n) const { return cls[n]; }

  // Unknowns at (zeta_n, s^1).
  //  level s^0, slice n     : Rr0, Rp0, e0 (E_z), b0 (B_z)
  //  level s^1, slice n+1   : Rr1, Lr1, Rp1, Lp1, e1, b1
  //  currents: jr, jp, jz = j^1 at half-cell n+1/2;  Jr, Jp = four-point sums
  //  output (level s^1, slice n): Rr, Lr, Rp, Lp, e, b
  // Discrete Poisson problem  -Lap^c_m psi = rho  on the centres (ghost psi = 0
  // beyond c_{N-1}), with the same operator as the field update, so that the
  // fields E = -grad psi, B = z x E are an exact discrete equilibrium.
  void poisson(int m, const cd* rho, cd* psi);
  // the tridiagonal matrix -Lap^c_m of the Poisson problem (sub a, diagonal b, super cc)
  const Tridiag& poissonMatrix(int m) const { return Tp[m]; }
  // E_r (faces) = -G psi,  E_phi (centres) = -(i m / c) psi
  void gradient(int m, const cd* psi, cd* Er, cd* Ep) const;

  // slice n (step h_n to slice n+1)
  void solveSlice(int m, int n, const cd* Rr0, const cd* Rp0, const cd* e0, const cd* b0, const cd* Rr1, const cd* Lr1,
                  const cd* Rp1, const cd* Lp1, const cd* e1, const cd* b1, const cd* jr, const cd* jp,
                  const cd* jz, const cd* Jr, const cd* Jp, cd* Rr, cd* Lr, cd* Rp, cd* Lp, cd* e, cd* b);

 private:
  friend class Engine;        // uploads the factorized matrices
  friend struct EngineState;
  const RadialGrid* g = nullptr;
  double D = -1.0;
  int M = 0, N = 0;
  std::vector<double> hs;       // half-cell lengths of the window
  std::vector<double> hCls;     // distinct half-cell lengths
  std::vector<int> cls;         // class of half-cell n
  int nCls = 0;
  int ver = 0;
  std::vector<Tridiag> Tf, Tc;  // per class and mode (index cls * (M+1) + m): faces (B_z), centres (E_z)
  std::vector<Tridiag> Tp;      // per mode: -Lap^c_m (Poisson problem)
  std::vector<cd> Er1, Ep1, Bp1, Br1, Pr, Pp, Qr, Qp, tf1, tf2, tc1, tc2, rf, rc;

  // discrete operators (Section 5.2)
  void opG(const cd* e, cd* out) const;                            // centres -> faces
  void opg(const cd* b, cd* out) const;                            // faces -> centres
  void opMf(const cd* b, cd* out, cd im) const;                    // (i m / rho) on faces
  void opCE(const cd* Er, const cd* Ep, cd* out, cd im) const;     // curl_z E on faces
  void opCB(const cd* Bp, const cd* Br, cd* out, cd im) const;     // curl_z B on centres
};

}  // namespace gem
