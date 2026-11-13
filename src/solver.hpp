// ============================================================================
//  solver.hpp -- field solver for one zeta-slice and one azimuthal mode
//                (finite differences in r, Section 5.3 of the paper)
// ============================================================================
#pragma once
#include <vector>

#include "fields.hpp"
#include "grid.hpp"
#include "smoothing.hpp"
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
  // implicitness of the R update in s: R^1 = R^0 + D[(1-theta) F^0 + theta F^1] (0.5: trapezoid)
  double theta = 0.5;
  // numerics.theta_stiff = c: the fields E_z, B_z enter the R update with an implicitness per radial eigenmode,
  //     theta_k = (theta + mu_k) / (1 + mu_k),   mu_k = -a lambda_k,   a = c h Delta,
  // lambda_k < 0 the eigenvalues of the radial Laplacian of E_z (centres) or B_z (faces) of the mode: theta for the
  // modes the step resolves (c h Delta |lambda| << 1, all radiation), -> 1 (L-stable) for the stiff near-field modes
  // of thin beam structures. As operators, Theta = S^-1 (theta - a Lap) with S = 1 - a Lap; the implicit matrices
  // become S (16 - h^2 Lap) - 4 h Delta Lap (theta - a Lap) (pentadiagonal). 0: off (scalar theta).
  double thetaStiff = 0.0;
  bool stiff() const { return thetaStiff > 0.0; }
  // perfectly conducting wall at r = rw (block 'wall'; rw < 0: none). Inside lie the faces 0..K and the centres
  // 0..K-1, rho_K <= rw < rho_{K+1}: E_z = 0 at rw (the ghost centre c_K := rw), zero flux of B_z through rw
  // (E_phi = 0: the control volume of face K ends at rw), all fields outside are zero. Without a wall K = N,
  // E_z = 0 at the ghost c_N and B_z(rho_N) = 0, as before. Returns true if the matrices change.
  // slope = d rw / ds of the step: on a sloped wall the tangential field E_z + slope E_r vanishes, so E_z at the
  // wall is -slope E_r (E_r of face K at level s^0, same slice), kept in the ghost centre e[K]. With E_z = 0 the
  // moving wall would swallow the radial Poynting flux S_r = slope S_z that pushes the field inwards.
  bool setWall(double rw, double slope = 0.0);
  bool hasWall() const { return wNeu; }
  int wallK() const { return wK; }    // last face inside (N without a wall)
  double wallR() const { return wR; } // ghost centre of E_z (the wall radius; c_N without a wall)
  double wallAf() const { return wAf; }
  // zero the fields of one slice outside the wall (faces > K, centres >= K; E_z: centres > K, e[K] = E_z at the wall)
  void maskWall(cd* Rr, cd* Lr, cd* Rp, cd* Lp, cd* e, cd* b) const;
  // ---- dielectric (block 'dielectric'): chi = eps - 1 of the level s^1 (chi1) and s^0 (chi0), on the centres
  // (cell average: E_z, E_phi are tangential) and the faces (harmonic mean of the neighbouring centres: E_r is
  // normal). The medium enters as the polarization current J = d_t P = -d_zeta (chi E), integrated over every
  // half-cell with the trapezoidal rule like L and E_z,
  //   (J_n + J_{n+1}) / 2 = -chi (E_{n+1} - E_n) / h_n,
  // in the R update (s-trapezoid, weights 1 - theta, theta), in the L update (the half-cell current) and in E_z
  // (eps d_zeta E_z). J of the transverse components is a running value along zeta (J = 0 ahead of the window):
  // J^1 from the slices of the march, J^0 recomputed from the stored slices of s^0, so there is no extra state.
  // With a dielectric every mode has one pentadiagonal system for (B_z, E_z) interleaved (b_0, e_0, b_1, ...):
  // where chi varies with r, E_z and B_z of the modes m >= 1 couple (hybrid modes). chiC/chiF empty: vacuum.
  // Returns true if the matrices change.
  bool setDielectric(const std::vector<double>& chiC1, const std::vector<double>& chiF1,
                     const std::vector<double>& chiC0, const std::vector<double>& chiF0);
  bool dielectric() const { return diel; }
  double chiC(int i) const { return dChiC1[i]; }
  double chiF(int i) const { return dChiF1[i]; }
  // running J at the front of a march (zero ahead of the window, or the values of the rank ahead)
  void dielBegin();
  // J^0 at slice n of mode m from slices n, n+1 of level s^0 (the fields R, L of both)
  void dielLevel0(int m, int n, const cd* Rr0, const cd* Lr0, const cd* Rp0, const cd* Lp0, const cd* Rr0n,
                  const cd* Lr0n, const cd* Rp0n, const cd* Lp0n);
  // slice n is final: its J become the values at n+1 for the next slice
  void dielNext();
  // the running J at slice n+1 as one vector (J^0_r, J^0_phi, J^1_r, J^1_phi): the MPI message to the rank behind
  void dielGet(std::vector<cd>& out) const;
  void dielSet(const std::vector<cd>& in);
  // the R correction of a final half-cell (predicted -> final current, change dP of the known part of R at a node)
  // in a medium: R and L change by fR dP and fL dP (E = (R + L)/2 follows the local implicit relation)
  void dielCorrection(int cls, bool face, int i, double& fR, double& fL) const;
  // running J (for the MPI messages): r components (M+1)(N+1), phi components (M+1)N
  std::vector<cd> dJ0r, dJ0p, dJ1r, dJ1p;     // values at slice n+1 (ahead)
  std::vector<cd> dJ0rc, dJ0pc, dJ1rc, dJ1pc;  // values at slice n (being solved)
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
  void solveSlice(int m, int n, const cd* Rr0, const cd* Lr0, const cd* Rp0, const cd* e0, const cd* b0, const cd* Rr1, const cd* Lr1,
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
  // wall (setWall): last face inside, ghost centre, control volume of face K, wall present (B_z: zero flux)
  int wK = 0;
  double wR = 0.0, wAf = 0.0, wSlope = 0.0;
  std::vector<double> wA;  // per class and mode: w a of centre row K-1 (the wall value of E_z enters its right-hand side)
  bool wNeu = false, wallChanged = false;
  // coefficients of the radial Laplacians in row i (flux coefficient a to the next node or the wall, bb to the
  // previous node, m^2 / r^2; up / low: the neighbour is an unknown; in: the row lies inside the wall)
  struct Row {
    double a = 0.0, bb = 0.0, mm = 0.0;
    bool up = false, low = false, in = true;
  };
  Row rowC(int m, int i) const;  // centres (E_z, Poisson)
  Row rowF(int m, int i) const;  // faces (B_z)
  void buildPoisson();
  std::vector<Tridiag> Tf, Tc;  // per class and mode (index cls * (M+1) + m): faces (B_z), centres (E_z)
  std::vector<Tridiag> Tp;      // per mode: -Lap^c_m (Poisson problem)
  std::vector<cd> Er1, Ep1, Bp1, Br1, Pr, Pp, Qr, Qp, tf1, tf2, tc1, tc2, rf, rc;
  // thetaStiff: per mode the radial Laplacians (unfactorized, rows as Tc, Tf); per class and mode a = c h Delta,
  // the factorized implicit matrices (Bc centres, Bf faces) and S = 1 - a Lap (Sc, Sf)
  std::vector<Tridiag> Lc, Lf;
  struct Stiff {
    double a = 0.0;
    Banded5 Bc, Bf;
    Tridiag Sc, Sf;
  };
  std::vector<Stiff> St;
  std::vector<cd> e0t, b0t, et, bt, sx;
  void lapRows(const Tridiag& L, const cd* x, cd* y) const;  // y = Lap x on the rows of L (absolute indices)
  // y = Theta x (one = false) or (1 - Theta) x (one = true) on the rows of L; other entries of y are zero
  void thetaApply(Stiff& S, bool faces, int m, bool one, const cd* x, cd* y);
  void stiffSolve(Stiff& S, bool faces, int m, const cd* rhs, cd* x);  // the pentadiagonal implicit solve

  // dielectric: chi per node of s^1 and s^0, matrices per class and mode (pentadiagonal, interleaved b, e)
  bool diel = false, dielChanged = false;
  std::vector<double> dChiC1, dChiF1, dChiC0, dChiF0;
  std::vector<Banded5> Bd;
  std::vector<cd> dA, dB, dAp, dBp, dEr, dEp, dRr, dLr, dRp, dLp, dPkr, dQkr, dPkp, dQkp, dE, dBb, dU, dV, dX;
  void buildDiel();
  // transverse fields of slice n (class h) from e (centres, e[K] = the wall value), b (faces) and the known parts
  // Pk, Qk of R and L (null: zero): R, L, E of both components
  void dielTransverse(int m, double h, const cd* e, const cd* b, const cd* Pkr, const cd* Qkr, const cd* Pkp,
                      const cd* Qkp, cd* Rr, cd* Lr, cd* Rp, cd* Lp, cd* Er, cd* Ep);
  // left-hand side of the slice equations: 16 b + 8h CE(E), 16 eps e - 8h CB(B) for the transverse fields of
  // (e, b) without known parts
  void dielOperator(int m, double h, const cd* e, const cd* b, cd* outE, cd* outB);
  bool rowInE(int i) const { return i < wK; }
  bool rowInB(int m, int i) const { return !(i > wK || (i == wK && !wNeu) || (m > 0 && i == 0)); }
  void solveSliceDiel(int m, int n, const cd* Rr0, const cd* Lr0, const cd* Rp0, const cd* e0, const cd* b0,
                      const cd* Rr1, const cd* Lr1, const cd* Rp1, const cd* Lp1, const cd* e1, const cd* b1,
                      const cd* jr, const cd* jp, const cd* jz, const cd* Jr, const cd* Jp, cd* Rr, cd* Lr, cd* Rp,
                      cd* Lp, cd* e, cd* b);

  // discrete operators (Section 5.2)
  void opG(const cd* e, cd* out) const;                            // centres -> faces
  void opg(const cd* b, cd* out) const;                            // faces -> centres
  void opMf(const cd* b, cd* out, cd im) const;                    // (i m / rho) on faces
  void opCE(const cd* Er, const cd* Ep, cd* out, cd im) const;     // curl_z E on faces
  void opCB(const cd* Bp, const cd* Br, cd* out, cd im) const;     // curl_z B on centres
};

}  // namespace gem
