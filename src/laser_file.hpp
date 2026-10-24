// ============================================================================
//  laser_file.hpp -- laser pulses from and to LASY / openPMD envelope files
//
//  File format (LASY, https://github.com/LASY-org/lasy; openPMD extension
//  LaserEnvelope): mesh record "laserEnvelope" (a scalar, complex dataset) of
//  one iteration, the complex envelope of the transverse electric field at a
//  plane z = const, as a function of time t and of (r, theta) or (x, y):
//      E_x + i E_y  ->  E_{x,y}(x, y, t) = Re[ p_{x,y} env(x, y, t) e^{-i omega_L t} ]
//  attributes: angularFrequency (omega_L, rad/s), polarization (p_x, p_y, real or
//  complex), envelopeField ("electric_field" [V/m] or "normalized_vector_potential"),
//  geometry thetaMode (axes t, r; components m = 0, cos(m theta), sin(m theta):
//  env = C_0 + sum_m C_m cos(m theta) + S_m sin(m theta)) or cartesian (axes t, y, x).
//  The time axis is the lab time at the plane (LASY adds the position / c).
//
//  Mapping onto GEM-PIC-RZ (GT2): a level s is the plane z = s and its slice zeta
//  is the time t = s - zeta. A file written at the plane z = z_f therefore IS a
//  level of the code; its time t_c is put at zeta_c:
//      tau(zeta) = t_c + zeta_c - zeta     (lab time of slice zeta at the file's plane)
//  At another level s the field follows from the exact propagation of every
//  plane wave (omega, k_perp) in vacuum or a uniform plasma,
//      k_z = sqrt(omega^2 - omega_p^2 - k_perp^2),   phase e^{i (k_z - omega)(s - z_f)}
//  at fixed tau (azimuthal modes: Hankel transforms, time: FFT).
//  The forward field of the code is R = 2 E_perp, L = 0.
// ============================================================================
#pragma once
#include <algorithm>
#include <cmath>
#include <complex>
#include <string>
#include <vector>

#include "units.hpp"

namespace gem {

using cd = std::complex<double>;

// content of an envelope file, converted to code units
struct LaserFileData {
  std::string path, mesh, software, field;  // field: envelopeField of the file
  long iteration = -1;
  bool cartesian = false;
  double omegaL = 1.0;        // carrier frequency [omega0]
  cd px = 1.0, py = 0.0;      // polarization vector (normalized)
  int nt = 0, nr = 0, nx = 0, ny = 0, ncomp = 1;  // ncomp: thetaMode components 2 K_f + 1
  double t0 = 0, dt = 1;      // time axis [1/omega0] (lab time at the file's plane)
  double r0 = 0, dr = 1;      // thetaMode: radial axis [1/k0]
  double x0 = 0, dx = 1, y0 = 0, dy = 1;  // cartesian axes [1/k0]
  // envelope of the electric field [m c omega0 / e]; thetaMode [comp][t][r], cartesian [t][y][x]
  std::vector<cd> env;
  double peak = 0;            // largest |env| over the file (all theta)
};

// Reads the mesh `mesh` (default "laserEnvelope") of `iteration` (< 0: the last).
// Throws std::runtime_error with the file name in the message.
LaserFileData readLaserFile(const std::string& path, const std::string& mesh, long iteration, const Units& u);

// in-place FFT of a power-of-two length: a_q <- sum_l a_l exp(sign 2 pi i q l / N)
void fft(std::vector<cd>& a, int sign);

// Envelope modes a_k, k = -K..K (a[k + K]; env(theta) = sum_k a_k e^{-i k theta}) of a pulse with the
// polarization (px, py) and carrier omegaL, at the lab time tau -> R_r,m and R_phi,m, m = 0..M of the code
// (F = F_0 + 2 Re sum F_m e^{i m theta}, forward wave R = 2 E_perp)
void envelopeToR(const cd* a, int K, cd px, cd py, double omegaL, double tau, int M, cd* Rr, cd* Rp);

// The field of a file laser at any level s and radius r (see the top of this file).
class LaserFileField {
 public:
  struct Params {
    double scale = 1.0;        // factor on the envelope
    double zeta_center = 0.0;  // window position of the file time t_center
    double t_center = 0.0;     // [1/omega0]
    double file_z = 0.0;       // plane of the file along z (simulation coordinates)
    double background = 0.0;   // omega_p^2 [n_c] of the medium of the propagation
    int n_theta = 0;           // cartesian files: azimuthal samples (0: automatic)
    int padding = 2;           // time window of the FFT / length of the file's
  };
  LaserFileField(const LaserFileData& d, const Params& p, int M, double rNeed);

  // R_r,m (comp 0) or R_phi,m (comp 1) at the radii x window positions zetas, level s:
  // out[(i * zetas.size() + n) * (M + 1) + m] (overwritten). direct: interpolate the file
  // (only at its own plane s = file_z); else exact propagation from the file's plane.
  void field(double s, const std::vector<double>& radii, const std::vector<double>& zetas, int comp,
             std::vector<cd>& out, bool direct) const;
  bool atFilePlane(double s) const { return std::fabs(s - P.file_z) <= 1e-9 * std::max(1.0, std::fabs(P.file_z)); }

  mutable std::vector<std::string> info;  // summary lines (mode content, spectral band, propagator grid)
  double modeLoss = 0.0;           // fraction of the energy in azimuthal modes the code does not keep
  double edgeField = 0.0;          // largest |env| at the outer radius / at the ends of the time axis, relative

 private:
  Params P;
  int M = 1, K = 2;                // envelope modes |k| <= K = M + 1 contribute to the code modes 0..M
  double omegaL = 1.0;
  cd px = 1.0, py = 0.0;
  int nt = 0, nr = 0;
  double t0 = 0, dt = 1, r0 = 0, dr = 1;
  std::vector<cd> A;               // envelope modes [(k + K) * nt + l] * nr + j, uniform r_j = r0 + j dr
  std::vector<int> active;         // indices k + K with nonzero data
  double rNeed = 0.0;
  double tau(double zeta) const { return P.t_center + P.zeta_center - zeta; }

  // spectral representation (prepared on first use)
  mutable bool prepared = false;
  mutable int nt2 = 0, qlo = 0, nq = 0, nk = 0;
  mutable double Tstart = 0, dOm = 0, dk = 0;
  mutable std::vector<cd> G;       // [(a * nq + q) * nk + p]: Hankel-FFT coefficients x inverse weights (a: active)
  mutable std::vector<double> beta;// [q * nk + p]: k_z - omega (propagating), NaN: dropped
  void prepare() const;
  void interp(double r, double t, cd* a) const;  // envelope modes at (r, t) by linear interpolation
};

}  // namespace gem
