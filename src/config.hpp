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
#include <vector>

#include "deck.hpp"
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
  double fraction() const { return Z * atoms; }
  double qm() const { return double(Z) / mass; }  // initial charge-to-mass ratio [e/m_e]
};

// ---------------------------------------------------------------- bunch
struct BunchConfig {
  enum class Shape { Gaussian, Flattop, Function };
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
  double x0 = 0, y0 = 0;    // transverse offset
  double z_waist = 0;       // position of the beta waist (z = 0: at start)
  int macroparticles = 10000;
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
  bool hdf5_single = false;  // store floating-point data in single precision
  int hdf5_compression = 0;  // gzip level 0..9 (0: none)
  // checkpoints for restarts (checkpoint_NNNNN.h5, complete state in double precision)
  int checkpoint_every = 0;     // every N steps, and at the end of the run (0: none)
  int checkpoint_keep = 2;      // keep the newest N checkpoint files (0: all)
  double walltime_hours = 0.0;  // stop with a checkpoint after this wall-clock time (0: no limit)
};

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
  int beam_filter = 1;       // (1,2,1) filter passes along zeta on the beam current (Nyquist Cherenkov)
  bool beam_vay = true;      // beam list pushed with the Vay pusher (false: Boris)
  unsigned long seed = 12345;
  double beam_dz_max = 10.0; // max sub-step of the beam push along s
  double vz_min = 0.8;       // beam particles slower than this return to streaming
  double vz_beam = 0.9;      // particles reaching the next s-level join the beam list if v_z >= vz_beam
};

// ---------------------------------------------------------------- all
struct Config {
  Units units;
  GridConfig grid;
  std::vector<LaserConfig> lasers;
  PlasmaConfig plasma;
  std::vector<IonConfig> ions;
  std::vector<BunchConfig> bunches;
  std::vector<UndulatorConfig> undulators;
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
