// ============================================================================
//  config.hpp -- typed run configuration, built from the input deck
//
//  One structure per deck block:
//     grid      GridConfig       (exactly one)   gridding, steps, modes
//     laser     LaserConfig      (0 ... N)       laser pulses
//     plasma    PlasmaConfig     (0 or 1)        background plasma profile
//     ion       IonConfig        (0 ... N)       mobile ion species of the plasma
//     bunch     BunchConfig      (0 ... N)       particle bunches
//     undulator UndulatorConfig  (0 ... N)       static undulator fields
//     output    OutputConfig     (0 or 1)        diagnostics
//     numerics  NumericsConfig   (0 or 1)        algorithm parameters
//
//  All values are stored in CODE UNITS (see units.hpp); the deck may use
//  physical units. The full key reference is docs/input_deck.md.
// ============================================================================
#pragma once
#include <string>
#include <algorithm>
#include <cmath>
#include <memory>
#include <vector>

#include "deck.hpp"
#include "laser_file.hpp"
#include "expr.hpp"
#include "mapping.hpp"
#include "units.hpp"

namespace gem {

// ---------------------------------------------------------------- grid
struct GridConfig {
  int M = 1;  // highest azimuthal mode, modes 0..M

  // radial grid, faces at r = R(i), centres at R(i + 1/2):
  //   table:     dR/dxi = dr(r), piecewise linear between reference points
  //              (free-hand grid, see mapping.hpp)
  //   stretched: dR/dxi = dr_min + (dr_max - dr_min) (1 + tanh((xi - xi0)/w)) / 2
  //   uniform:   dR/dxi = dr
  enum class Radial { Uniform, Stretched, Table };
  Radial radial = Radial::Stretched;
  bool stretched = true;              // Radial::Stretched (kept for existing code)
  int Nr = 100;
  double dr_min = 0.3, dr_max = 3.0;  // smallest / largest radial step
  double xi0 = 40.0;                  // stretched: transition position [cells]
  double wcells = 8.0;                // stretched: transition width [cells]
  StepMapping rtable;                 // table: reference points and steps
  double r_max = 0.0;                 // resulting outer radius (derived)

  // longitudinal window zeta = 0 (rear) ... zeta_Nz (front), slices zeta_n:
  //   uniform: zeta_n = n dzeta
  //   table:   d zeta/dn = dzeta(zeta), piecewise linear between reference
  //            points (free-hand grid, the same mapping as the radial table)
  int Nz = 1200;
  bool zetaTable = false;
  double dzeta = 0.25;      // uniform: the step; table: the smallest step
  double dzeta_max = 0.25;  // table: the largest step
  StepMapping ztable;       // table: reference points and steps

  // propagation: steps s^0 -> s^1 = s^0 + Delta(s^0) along s = z
  //   uniform: Delta = dz
  //   table:   the levels follow ds/dk = dz(s), piecewise linear between
  //            reference points; the step from s is X(xi(s) + 1) - s
  bool zTable = false;
  double dz = 40.0;     // uniform: the step; table: the smallest step
  double dz_max = 40.0; // table: the largest step
  StepMapping stable;   // table: reference points and steps
  double z_end = 4000;  // total propagation distance

  // radial absorber for R, active for r > r_absorber
  double r_absorber = 160.0;
  double absorber_strength = 0.2;

  double window() const { return zetaNode(Nz); }
  // a run ends at the first level s >= z_end - endTol() (relative: rounding of many steps)
  double endTol() const { return 1e-9 * std::max(1.0, std::fabs(z_end)); }
  double mapping(double xi) const;  // R(xi)
  double zetaNode(int n) const;     // zeta_n
  double zetaStep(double zeta) const;  // local step dzeta at zeta
  double stepAt(double s) const;    // Delta of the step that starts at s
  int stepsFrom(double s) const;    // number of steps from s to z_end
};

// ---------------------------------------------------------------- laser
struct LaserConfig {
  enum class Pol { Linear, Circular };
  double a0 = 1.0;
  double omega = 1.0;      // laser frequency / omega0
  double waist = 10.0;     // w0: 1/e field radius at focus
  double length_1e = 10.0; // 1/e half length of the field envelope
  double zeta_center = 0;  // envelope centre in the window (envelope_function: phase reference)
  double z_focus = 0;      // focal plane (z = 0 at simulation start)
  double cep = 0;          // carrier-envelope phase [rad]
  Expr envelope_fn;        // optional envelope f(zeta) of the field (replaces the Gaussian in zeta)
  Expr profile_fn;         // optional transverse profile g(r) at z = 0 (replaces the Gaussian beam)
  Pol pol = Pol::Linear;
  double pol_angle = 0;    // linear: angle of E to the x axis [rad]
  int helicity = 1;        // circular: +1 or -1
  double zR() const { return 0.5 * omega * waist * waist; }
  // laser from a LASY / openPMD envelope file (laser.file; see laser_file.hpp): the envelope at the plane
  // z = file_z, its time t_center at zeta_center; inject: box (initial field only) or boundary (also fed through
  // the radial absorber layer at every step, total-field / scattered-field, like the flying focus)
  std::shared_ptr<const LaserFileData> file;
  double file_scale = 1.0;   // factor on the file's envelope (laser.a0 given: peak |E| = a0 omega_L)
  double t_center = 0.0;     // [1/omega0], lab time at the file's plane
  double file_z = 0.0;
  bool boundary = false;
  double background = 0.0;   // omega_p^2 [n_c] of the uniform medium of the propagation
  int n_theta = 0;           // cartesian files: azimuthal samples (0: automatic)
  int time_padding = 2;
  double file_mode_loss = 0.0;  // energy fraction in envelope modes beyond m_max + 1 (summary)
  LaserFileField::Params fileParams() const {
    LaserFileField::Params p;
    p.scale = file_scale; p.zeta_center = zeta_center; p.t_center = t_center; p.file_z = file_z;
    p.background = background; p.n_theta = n_theta; p.padding = time_padding;
    return p;
  }
};

// ---------------------------------------------------------------- flying focus
// Luminal Bessel pulse (focus-wave mode): a converging conical wave of cone angle theta whose cone angle
// depends on the frequency, k_perp(omega)^2 = 2 omega kappa - kappa^2 with kappa = omega0 (1 - cos theta), so
// that every frequency has k_z = omega - kappa: the pattern moves with the window at c (no dephasing),
//     E_x(r, zeta, s) = Re[ e^{-i kappa s} sum_omega a(omega) J0(k_perp(omega) r) e^{i omega (zeta - zeta_c) + i cep} ],
// In a uniform background plasma (omega_p^2 = background) k_perp^2 = 2 omega kappa - kappa^2 - omega_p^2.
// a(omega) Gaussian: on the axis E_x = a0 omega0 exp(-(zeta - zeta_c)^2 / L^2) cos(...). It is fed through the
// radial absorber layer (total-field / scattered-field): there R relaxes to the incident field, and only the
// difference (what the plasma radiates or scatters) is absorbed. Default: linear polarization along x, mode m = 1.
// Several components (beat_* keys): U = sum_c U_c, U_c = A_c e^{i psi_c} e^{-i kappa_c s} sum_omega a(omega)
// J_|l_c|(k_c(omega) r) e^{i omega (zeta - zeta_c) + i cep}, Bessel order l_c (azimuthal phase e^{i l_c phi}), cone
// angle theta_c; every component moves at c, two of them beat with the period 2 pi / |kappa_1 - kappa_0| in s.
// Physical field: linear E_x = Re sum_c U_c e^{i l_c phi}; circular with helicity h E_x + i h E_y = sum_c U_c e^{i l_c phi}.
// With W = e^{-i phi} (E_x + i E_y) = sum_n w_n e^{i n phi}: R_r,m = w_m + conj(w_-m), R_phi,m = -i (w_m - conj(w_-m))
// (as envelopeToR in laser_file.cpp). Helicity +1 gives w_{l-1} = U_c, helicity -1 gives w_{-(l+1)} = conj(U_c);
// linear polarization is the mean of both. terms() lists these contributions.
struct FlyingFocusConfig {
  struct Component {
    int order = 0;         // Bessel order l
    double theta = 0.05;   // cone angle at omega0 [rad]
    double amp = 1.0;      // amplitude relative to a0
    double phase = 0.0;    // psi_c
    double kappa() const { return 1.0 - std::cos(theta); }  // omega0 (1 - cos theta), code units
  };
  // contribution of component comp to mode m: Y = U_c or conj(U_c) (conj), U_c = G_c e^{-i kappa_c s} at level s;
  // R_r,m += w Y, R_phi,m += i rps w Y
  struct Term { int comp, m; bool conj; double w, rps; };
  bool present = false;
  double a0 = 1.0;
  double length_1e = 10.0; // 1/e half length of the on-axis field envelope
  double zeta_center = 0;
  double cep = 0;
  double background = 0;   // density [n_c] of a uniform plasma in which the incident field propagates (0: vacuum)
  bool circular = false;
  int helicity = 1;
  std::vector<Component> comps;  // comps[0]: cone_angle, order (amplitude 1); then the beat_* components
  double thetaMax() const {
    double t = 0.0;
    for (const auto& c : comps) t = std::max(t, c.theta);
    return t;
  }
  // all terms, also those of modes m > m_max (the caller drops them)
  std::vector<Term> terms() const {
    std::vector<Term> T;
    const double w = circular ? 1.0 : 0.5;
    for (int c = 0; c < int(comps.size()); ++c)
      for (int h : {1, -1}) {
        if (circular && h != helicity) continue;
        const int n = h > 0 ? comps[c].order - 1 : -(comps[c].order + 1);
        const bool cu = h < 0;                            // w_n = conj(U_c)
        if (n >= 0) T.push_back({c, n, cu, w, -1.0});     // w_m, m = n
        if (n <= 0) T.push_back({c, -n, !cu, w, 1.0});    // conj(w_-m), m = -n
      }
    return T;
  }
};

// ---------------------------------------------------------------- plasma
struct PlasmaConfig {
  enum class Radial { Uniform, Parabolic, Table };
  bool present = false;
  double n0 = 0.0;                 // reference density [n_c]
  bool free_electrons = true;      // load free electrons with density n0 f(z) g(r)
  std::vector<double> z_points;    // longitudinal profile (empty: constant 1)
  std::vector<double> z_values;
  Radial radial = Radial::Uniform;
  double r_max = 1e300;            // plasma column radius
  double smooth_length = 0.0;      // radial smoothing length a of the plasma sources (0: off)
  double smooth_length_m0 = -1.0;  // a of the mode m = 0 if >= 0 (else smooth_length)
  double channel_depth = 0.0;      // parabolic: n ~ 1 + depth r^2 / r_channel^2
  double r_channel = 1.0;
  std::vector<double> r_points, r_values;  // radial table
  // macroparticles
  int ppc_r = 2, ppc_phi = 8;              // per radial cell / per ring
  int ppc_r_axis = 4, ppc_phi_axis = 16;   // for cells with rho_{i+1} <= r_axis_refine
  double r_axis_refine = 10.0;

  Expr fn;                         // density_function(x, y, z, r): further factor (empty: 1)
  double longitudinal(double z) const;  // relative profile f(z)
  double radialFactor(double r) const;  // relative profile g(r)
  double density(double r, double z) const;  // n0 f(z) g(r), 0 outside r_max (axisymmetric part)
  double density(double x, double y, double z) const;  // ... times density_function(x, y, z, r)
};

// ---------------------------------------------------------------- ions
// An ion species of the plasma, represented by macroparticles. Its atom
// density follows the plasma profile,  n_atom(r, z) = atoms * n0 f(z) g(r),
// and it is loaded at the positions of the plasma electrons. At the start it
// carries the share  fraction() = Z * atoms  of the positive charge that
// neutralizes the free electrons; the charge of all ion blocks together must
// not exceed that of the free electrons, the rest is an immobile background.
// With field ionization every macroparticle has its own charge state
// Z ... Zmax; each ionization event creates a free electron.
struct IonConfig {
  std::string name;
  std::string element;         // chemical symbol, if given
  int Znuc = 0;                // atomic number (0: not known, ion given by its mass)
  int Z = 1;                   // initial charge state (0 = neutral atoms, needs ionization)
  double mass = 1.0;           // ion mass [m_e]
  double atoms = 1.0;          // atom density / plasma reference density n0
  bool mobile = true;          // pushed by the fields (false: fixed, e.g. heavy dopant)
  bool ionization = false;     // ADK field ionization of the levels Z ... Zmax-1
  int Zmax = 1;                // highest charge state reachable by ionization
  std::vector<double> ip_eV;   // ionization energies, index = charge state k (k -> k+1)
  // own profile (deck: profile = own): atom density prof.density(x, y, z) instead of
  // atoms x the plasma profile; loaded with Z electrons per atom (neutral)
  bool own = false;
  PlasmaConfig prof;
  double fraction() const { return own ? 0.0 : Z * atoms; }  // share of the plasma's positive charge
  double qm() const { return double(Z) / mass; }  // initial charge-to-mass ratio [e/m_e]
};

// ---------------------------------------------------------------- bunch
struct BunchConfig {
  enum class Shape { Gaussian, Flattop, Function, File };
  std::string name;
  double qm = -1.0;         // charge-to-mass ratio: -1 electron, +1 positron
  double charge = 0.0;      // total charge [code units], signed
  double charge_pC = 0.0;   // |charge| in pC, for printing
  double gamma = 1000.0;
  double energy_spread = 0; // relative rms of p_z
  double sigma_x = 1, sigma_y = 1;
  double emittance_x = 0, emittance_y = 0;  // normalized rms emittance
  Shape shape = Shape::Gaussian;
  double sigma_z = 1;       // Gaussian: rms length
  double length = 1;        // flat-top: full length
  double zeta_center = 0;
  Expr profile_z;           // Function: current profile f(zeta) on [zeta_min, zeta_max]
  double zeta_min = 0, zeta_max = 0;
  Expr gamma_fn;            // optional gamma(zeta) (energy chirp); empty: gamma
  // optional transverse profile f(x, y, r, zeta) >= 0, sampled by rejection in the disk r < r_extent
  // (replaces the Gaussian of sigma_x, sigma_y; those then hold the rms of the profile, for messages)
  Expr profile_r;
  double r_extent = 0.0, profile_r_max = 0.0;  // sampling disk; bound of f used for the rejection
  // optional mean transverse momenta u_x, u_y (x, y, r, zeta) [m c], added to the emittance spread; the energy
  // is kept: u_z = sqrt(u_z,sampled^2 - u_x^2 - u_y^2)
  Expr ux_fn, uy_fn;
  double x0 = 0, y0 = 0;    // transverse offset
  double z_waist = 0;       // position of the beta waist (z = 0: at start)
  double z_start = 0;       // the bunch is loaded (with its co-moving self-field) when the run reaches z_start
  int macroparticles = 10000;
  // File: particles read from an openPMD / native file (bunch.file), after selection and
  // placement: x, y, zeta, ux, uy, uz [code units] and |charge| per macroparticle [code units]
  struct FileParticle { double x, y, zeta, ux, uy, uz, q; };
  std::vector<FileParticle> file_particles;
  std::string file_info;              // for the summary
  double file_center = 0.0;           // centroid (zeta_center of the deck)
  double resample_smoothing = 0.0;    // jitter of resampled particles, in units of the rms per coordinate
  // shot noise: random positions (noise of the macroparticles), or Fawley's beamlets
  // with the bunching noise of the real number of electrons at noise_wavelength; quiet = the
  // beamlets without their displacements (quiet start: no bunching at harmonics h < beamlet)
  bool fawley = false;
  bool quiet = false;
  double noise_wavelength = 0.0;
  int beamlet = 4;
  // more macroparticles near the axis (lower noise of the modes m >= 1 there): centres with r < axis_radius are
  // sampled axis_refine times as often, each with 1 / axis_refine of the charge (same density)
  int axis_refine = 1;
  double axis_radius = 0.0;  // 0: the first four radial cells
};

// ---------------------------------------------------------------- undulator
// static undulator field acting on the particles (see external.hpp)
struct UndulatorConfig {
  enum class Type { Planar, Helical };
  Type type = Type::Planar;
  double period = 1.0;     // lambda_u
  double K = 1.0;          // deflection parameter (peak value)
  double z_start = 0.0;    // entrance (lab z)
  double length = 1.0;     // magnetic length, including the ramps
  double ramp = 0.0;       // linear ramps of the field amplitude at both ends
  bool focusing = true;    // planar: natural focusing (off-axis field)
  Expr K_fn;               // optional K(z) (taper); empty: K
  double Kat(double z) const { return K_fn.empty() ? K : K_fn({z}); }
  double ku() const;
  double B0() const { return K * ku(); }  // code units
  // on-axis resonant wavelength for an electron of energy gamma
  double resonance(double gamma) const;
};

// ---------------------------------------------------------------- transport elements
// static fields acting on all particles over [z_start, z_start + length]
struct MagnetConfig {
  enum class Kind { Quadrupole, PlasmaLens, Solenoid };
  Kind kind = Kind::Quadrupole;
  std::string label;       // "quadrupole[0]", ...
  double z_start = 0.0, length = 0.0;
  double gradient = 0.0;   // quadrupole, plasma lens: dB/dx [code units]; > 0 focuses electrons (in x for a quadrupole)
  double radius = 1e300;   // plasma lens: radius of the current (B_phi ~ 1/r outside)
  double field = 0.0;      // solenoid: B_z [code units]
  double ramp = 0.0;       // solenoid: linear fringe at both ends (inside the length)
};
// perfectly conducting lateral wall at r = r_w(z) (block 'wall'): a waveguide or beam pipe of variable radius.
// With s = z a level is the plane z = s, so the wall has one radius r_w(s) for the whole window at every step.
struct WallConfig {
  bool present = false;
  double radius = 0.0;                       // constant radius (no table, no function)
  std::vector<double> z_points, r_points;    // piecewise linear r_w(z), constant beyond the ends
  Expr fn;                                   // radius_function(z)
  double at(double z) const;                 // r_w(z), code units
};

// thin magnetic chicane at z: the beam particles move by r56 * delta - delay in zeta
struct ChicaneConfig {
  double z = 0.0, r56 = 0.0, delay = 0.0, gamma_ref = 1.0;
};

// ---------------------------------------------------------------- output / numerics
struct OutputConfig {
  std::string dir = "diag";
  int every = 25;          // write files every N steps (0: only at start and end)
  int map_stride = 4;      // every map_stride-th slice in 2D maps
  bool axis = true, maps = true, particles = true;
  // text: axis/map/beam ASCII files; hdf5: one file data_NNNNN.h5 per output
  // step with all modes, densities and the beam list; both: both
  bool text = true, hdf5 = false;
  // openpmd: openPMD 1.1.0 files <dir>/openpmd/openpmd_%T.h5 (T = step): fields in
  // thetaMode, interpolated to a uniform r-z grid; the beam list per species
  bool openpmd = false;
  double openpmd_dr = 0, openpmd_dzeta = 0;  // output grid (0: smallest steps of the grid)
  // openpmd_grid: uniform (interpolated, for any openPMD reader), native (the code's own nodes,
  // rectilinear-mesh records "native_*"), both
  bool openpmd_uniform = true, openpmd_native = true;
  // fel.txt every step: radiation energy, peak power, bunching (default: on with undulators)
  bool fel = false;
  double fel_wavelength = 0;  // bunching wavelength (default: resonance of undulator[0] for bunch[0])
  int fel_mode = -1;          // deck: -1 auto, 0 off, 1 on
  double fel_slice = 0;       // slices for the largest bunching (default 20 fel_wavelength)
  int fel_harmonics = 1;      // > 1: fel_harmonics.txt, power per harmonic 1..H and mode m = 0..M, bunching
  bool hdf5_single = false;  // store floating-point data in single precision
  int hdf5_compression = 0;  // gzip level 0..9 (0: none)
  // checkpoints for restarts (checkpoint_NNNNN.h5, complete state in double precision)
  int checkpoint_every = 0;     // every N steps, and at the end of the run (0: none)
  int checkpoint_keep = 2;      // keep the newest N checkpoint files (0: all)
  double walltime_hours = 0.0;  // stop with a checkpoint after this wall-clock time (0: no limit)
  // LASY envelope files of the laser (<dir>/lasy/laser_NNNNN.h5, laser_file.hpp): the forward field R/2 at the
  // plane z = s, demodulated at lasy_omega, for LASY and other codes
  bool lasy = false;
  int lasy_every = 0;           // every N steps (0: only at the end of the run)
  double lasy_dr = 0.0, lasy_dt = 0.0, lasy_r_max = 0.0;  // output grid (0: smallest steps; r_absorber)
  double lasy_omega = 0.0;      // carrier frequency of the envelope (0: laser[0], else omega0)
  double lasy_band = 0.5;       // kept spectrum |omega - omega_L| < lasy_band omega_L
};

enum { BEAM_FIELD_LEVEL = 0, BEAM_FIELD_MEAN = 1 };

struct NumericsConfig {
  // serial: reference implementation (simulation.cpp); openmp / cuda: parallel
  // engine (engine.cpp) on CPU threads / an NVIDIA GPU. Default: cuda in a
  // CUDA build, serial otherwise.
  enum class Backend { Serial, OpenMP, Cuda };
#if defined(GEMPIC_HAVE_CUDA)
  Backend backend = Backend::Cuda;
#else
  Backend backend = Backend::Serial;
#endif
  int iterations = 2;        // passes per slice (1 = predictor only)
  double beam_smooth_length = 0.0;  // radial smoothing of the beam current (as plasma.smooth_length); 0: off
  int beam_filter = 1;       // (1,2,1) filter passes along zeta on the beam current (Nyquist Cherenkov)
  double trap_gamma = 0.0;       // ... and only if its gamma is at least this (0: any gamma)
  double trap_advance = 1.0;     // a streaming particle leaves its level after advancing this many steps in z
  double theta = 0.5;            // implicitness of the R update in s (0.5: trapezoid)
  double theta_stiff = 0.0;      // per-mode implicitness theta_k -> 1 for c h dz |lambda_k| >> 1 (FieldSolver); 0: off
  // beam sources consistent with the slice equations: (1/4,1/2,1/4) of the beam currents in the half-cell equations
  // (E_z, L_r, L_phi) and the charge-conserving radial flux (mode 0) in the R update
  bool beam_consistent = false;
  int beam_field = -1;       // fields of the beam push: BEAM_FIELD_LEVEL (s^0), BEAM_FIELD_MEAN; -1 = auto
  bool beam_field_auto = true;  // beam_field chosen by the code (mean with undulators, else level)
  bool beam_vay = true;      // beam list pushed with the Vay pusher (false: Boris)
  unsigned long seed = 12345;
  double beam_dz_max = 10.0; // max sub-step of the beam push along s
  double vz_min = 0.8;       // beam particles slower than this return to streaming
  double vz_beam = 0.9;      // particles reaching the next s-level join the beam list if v_z >= vz_beam
  // tridiagonal solves of the field update in the parallel engine: Thomas (one thread per mode) or the
  // Thomas/PCR partition method (all threads of the team); auto: partition on the GPU, Thomas on CPU threads
  enum class Tridiag { Auto, Thomas, Partition };
  Tridiag tridiag = Tridiag::Auto;
};

// ---------------------------------------------------------------- all
struct Config {
  Units units;
  GridConfig grid;
  std::vector<LaserConfig> lasers;
  FlyingFocusConfig flying_focus;
  PlasmaConfig plasma;
  std::vector<IonConfig> ions;
  std::vector<BunchConfig> bunches;
  std::vector<UndulatorConfig> undulators;
  std::vector<MagnetConfig> magnets;    // quadrupole, plasma_lens, solenoid blocks
  std::vector<ChicaneConfig> chicanes;  // chicane blocks, sorted by z
  WallConfig wall;                      // block 'wall'
  std::map<std::string, double> constants;  // block 'constants' (code units)
  OutputConfig output;
  NumericsConfig numerics;

  // Build from a parsed deck: converts units, applies defaults, validates.
  // Throws DeckError with file:line information. Warnings go to `warnings`.
  static Config fromDeck(const Deck& deck, std::vector<std::string>* warnings = nullptr);

  // Human-readable summary in physical and code units.
  std::string summary() const;

  // smallest box for the fields of every bunch with v < c (see config.cpp)
  struct BoxAdvice {
    std::string name;
    double rmax_need = 0, rmax_now = 0;  // radial boundary
    std::string rmax_why;
    double gap_need = 0;                 // gap at both ends of the window
    double front_now = 0, rear_end = 0;  // gap ahead at the start, behind at the end of the run
    double slip = 0;                     // slippage of the bunch over the run
    double window_need = 0;              // smallest window length with both gaps
  };
  std::vector<BoxAdvice> boxAdvice() const;
};

}  // namespace gem
