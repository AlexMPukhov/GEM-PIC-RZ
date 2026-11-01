// ============================================================================
//  simulation.hpp -- GEM-PIC-RZ time stepping (Galilean transformation GT2:
//  s = z slow, zeta = z - t fast, c = 1)
//
//  Per step s^0 -> s^1 = s^0 + dz:
//   1. push the beam list (bunches + trapped electrons) along s, deposit j^1
//   2. load fresh plasma rings (electrons + mobile ions) at the front of the window
//   3. march zeta from the front to the rear; per slice predictor-corrector:
//      push streaming particles, deposit, solve the slice for all modes
//   4. streaming particles that reach z = s^1 join the beam list if v_z >= vz_beam;
//      slower ones continue streaming in the next step, inserted at their zeta
//      (beam particles with v_z < vz_min return to streaming the same way)
//   5. radial absorber, swap levels
// ============================================================================
#pragma once
#include <random>
#include <string>
#include <vector>

#include <cstdint>
#include <map>
#include <memory>

#include "config.hpp"
#include "engine.hpp"
#include "external.hpp"
#include "fields.hpp"
#include "grid.hpp"
#include "ionization.hpp"
#include "particles.hpp"
#include "smoothing.hpp"
#include "solver.hpp"
#include "zgrid.hpp"

namespace gem {

// moments of the particles of one species (combinable over MPI ranks)
struct Moments {
  long long n = 0;
  double Q = 0, x1 = 0, x2 = 0, u1 = 0, u2 = 0, xu = 0, y1 = 0, y2 = 0, v1 = 0, v2 = 0, yv = 0, z1 = 0, z2 = 0;
  // gamma: charge-weighted mean and sum of |q| (gamma - mean)^2, updated and combined centred (West / Chan et al.):
  // sum gamma^2 - (sum gamma)^2 / Q cancels to ~(dE/E)^2 of the sums (1e-11 at gamma = 1000, dE/E = 3e-6)
  double gmean = 0, gM2 = 0;
  double br = 0, bi = 0;  // sum |q| exp(i k zeta) at the FEL wavelength (bunching factor)
  void add(const Particle& p);
  void add(const Moments& o);
};

// everything the per-step console line and tables need, as partial sums /
// extrema that can be combined over ranks
struct DiagSums {
  double maxEx = 0, minEz = 0, maxEz = 0, nimax = 0;  // on-axis fields, mobile-ion density
  std::map<int, Moments> moments;                      // per species id of the beam list
  std::vector<double> freed;                           // per ion block: charge freed in the step
  // FEL diagnostics (output.fel): field energy of the modes m >= 1 (radiation of a centred beam)
  // and m = 0, largest power of the modes m >= 1 through a slice, and its zeta
  double Wrad = 0, W0 = 0, Pmax = 0, zPmax = 0;
  // bunching per bunch in slices of output.fel_slice: {Q, sum |q| cos, sum |q| sin} per slice
  std::map<int, std::vector<double>> slices;
  // power of the modes m >= 1 in slices of output.fel_slice: {sum P h, sum h} per slice
  std::vector<double> pslices;
  // resonant part: per slice, mode m >= 1 and radial point (faces 1..N, then centres 0..N-1) the sums
  // sum_n F_n exp(-+ i k zeta_n) h_n of R and L (8 doubles: R+, R-, L+, L- as re, im)
  std::vector<double> rslices;
  // harmonics (output.fel_harmonics = H > 1): per slice, harmonic h = 1..H, mode m = 0..M and radial point the
  // same 8 sums at +-h k; and the bunching of all bunches per slice and harmonic {sum |q| cos, sum |q| sin}
  std::vector<double> hslices, hbunch;
  void add(const DiagSums& o);
};

class PipelineRun;

// MPI pipeline: the part of the window a rank owns and stores (see pipeline.cpp). The state is built only there.
struct StorageBlock {
  int P = 1, r = 0;  // ranks, own rank
  int a = 0, b = 0;  // own slices / half-cells [a, b)
  int sa = 0, sb = 0;  // stored slices [sa, sb] and half-cells [sa, sb) (block and ghosts; not clipped)
  int g = 0;           // ghost width
};

class Simulation {
 public:
  // restartFile empty: initial state from the deck; else: state of a checkpoint
  // useEngine = false: no parallel engine even for backend openmp/cuda (MPI pipeline, threads in its loops)
  // blk: build and store only that part of the window (MPI pipeline); null: the whole window
  explicit Simulation(const Config& cfg, const std::string& restartFile = "", bool useEngine = true,
                      const StorageBlock* blk = nullptr);
  ~Simulation();

  void step();
  void diagnostics(bool writeFiles);
  bool finite() const;  // no NaN/Inf in the fields

  double z() const { return s0; }
  int stepNumber() const { return step_no; }
  long restreamedParticles() const { return n_restreamed; }
  const RadialGrid& grid() const { return g; }
  std::string backendInfo() const;  // backend and device
  void setDeckText(const std::string& t) { deckText = t; }  // stored in the HDF5 files
  // checkpoint of the complete state: <output.dir>/checkpoint_NNNNN.h5 (returns the path)
  std::string writeCheckpoint();
  // newest checkpoint_NNNNN.h5 in dir ("" if none)
  static std::string latestCheckpoint(const std::string& dir);

 private:
  friend class Engine;  // parallel backends (engine.cpp) share the state
  friend class PipelineRun;  // MPI pipeline along zeta (pipeline.cpp)
  Config cfg;
  RadialGrid g;
  Fields fold, fnew;            // levels s^0 and s^1
  Currents c0b, c0s, c1b, c1s;  // beam / streaming currents at s^0 and s^1
  FieldSolver solver;
  std::vector<Particle> stream, backup, beam;
  std::vector<cd> predR, predP, predRprev, predPprev, Jr, Jp, jr1, jp1, jz1;
  std::mt19937_64 rng;
  uint64_t ionKey = 0;                       // key of the ionization random numbers of the current pass
  std::vector<AdkRates> adk;                 // per ion block (empty levels without ionization)
  std::vector<Particle> born;                // electrons freed in the current slice
  std::vector<std::vector<double>> exitStates;  // per ion block: atoms per charge state leaving the window
  std::vector<double> freedCharge;           // per ion block: charge freed in the last step [code]
  double sliceFieldBound = 0.0;              // upper bound of |E| on the two slices of the current push
  double s0 = 0.0;
  int step_no = 0;
  long n_restreamed = 0, n_trapped_total = 0;
  std::vector<Particle> pending;      // particles joining the current march at their zeta (they wait
                                      // in the streaming list with status STATUS_WAITING)
  std::vector<Particle> pendingNext;  // ... joining the march of the next step
  double n_ref = 1.0;  // density for the n_e/n0 output
  std::unique_ptr<Engine> engine;  // parallel backend (null: serial code below)
  std::string deckText;            // deck after overrides (provenance of the output)

  // shortcuts to often used parameters
  int M = 0, Nz = 0;
  StorageBlock blk;   // own and stored part of the window (whole window without MPI)
  bool partial = false;  // only a block is stored
  // stored slices with data, [sLo, sHi) within 0..Nz-1
  int sLo = 0, sHi = 0;
  double Delta = 0;  // step of the current s^0 -> s^1 (set at the start of every step)
  ZetaGrid zg;       // longitudinal grid of the window
  ZetaView zv;       // view of zg
  std::vector<UndulatorV> undV;  // external fields (undulators)
  std::vector<double> undK;      // tapered undulators: tables of B0(z)
  std::vector<MagnetV> magV;     // transport elements
  ExternalV ext;                 // view of undV

  void beginStep();  // Delta of this step, matrices of the field solver

  // metal wall (block 'wall'): the solver boundary of level s, particles are lost at r >= rWall (r_max without a
  // wall); maskWall zeroes the stored fields outside the wall
  double rWall = 0.0;
  void setWall(double s);
  void maskWall(Fields& F) const;
  // the wall moved outwards past faces this step: wallFrom = its last face K0 at s^0 (-1: it did not). The field at
  // the old wall continues into the opened annulus (zeros there would make a jump of E_r, B_phi, B_z and radiate
  // high TM modes): faces K0+1..K1 take E_r, B_phi, B_z of face K0, the centres K0+1..K1 the wall value of E_z,
  // R_phi, L_phi stay 0 (E_phi is tangential).
  int wallFrom = -1;
  void extendWall(Fields& F) const;

  void initLasers();
  // flying focus (cfg.flying_focus): incident field G_c(r, zeta) of every component c at s = 0 on the faces
  // (c ffStrideF + n (N+1) + i) and centres (c ffStrideC + n N + i); at level s U_c = G_c e^{-i kappa_c s}, and
  // the terms ffTerms (modes m <= M) give R_r,m and R_phi,m (FlyingFocusConfig)
  std::vector<cd> ffGf, ffGc;
  std::vector<FlyingFocusConfig::Term> ffTerms;
  size_t ffStrideF = 0, ffStrideC = 0;
  int ffN0 = 0;  // first slice of ffGf, ffGc (stored slices only)
  void setupFlyingFocus();
  // incident R_r,m (face) or R_phi,m (centre) of the flying focus at level s
  cd ffInc(int m, bool face, int n, int i, double s) const {
    cd R = 0.0;
    for (const auto& t : ffTerms) {
      if (t.m != m) continue;
      const cd G = face ? ffGf[t.comp * ffStrideF + size_t(n - ffN0) * (g.N + 1) + i]
                        : ffGc[t.comp * ffStrideC + size_t(n - ffN0) * g.N + i];
      cd Y = G * std::polar(1.0, -cfg.flying_focus.comps[t.comp].kappa() * s);
      if (t.conj) Y = std::conj(Y);
      R += face ? t.w * Y : cd(0.0, t.rps * t.w) * Y;
    }
    return R;
  }
  // lasers from LASY / openPMD envelope files (laser.file): field evaluators (null for the other lasers)
  std::vector<std::unique_ptr<LaserFileField>> fileFields;
  // inject = boundary: incident field of the absorber layer at level incS, prepared for the slices [incN0, incN1):
  // incF[(m * Nz + n) * nLF + i - iaF] (faces i = iaF..N, R_r), incC[(m * Nz + n) * nLC + i - iaC] (centres, R_phi)
  bool hasIncident = false;
  int iaF = 0, iaC = 0, nLF = 0, nLC = 0, incN0 = 0, incN1 = 0;
  std::vector<cd> incF, incC;
  double incS = 0.0;
  void prepareIncident(double s, int n0, int n1);
  // (slices of the stored range: index n - nLo, nInc slices)
  int nLo = 0, nInc = 0;
  cd incidentF(int m, int n, int i) const {
    return (hasIncident && i >= iaF) ? incF[(size_t(m) * nInc + n - nLo) * nLF + i - iaF] : cd(0.0);
  }
  cd incidentC(int m, int n, int i) const {
    return (hasIncident && i >= iaC) ? incC[(size_t(m) * nInc + n - nLo) * nLC + i - iaC] : cd(0.0);
  }
  // LASY output (output.lasy): lab time of slice zeta at level s is s - zeta + lasyTref (the time base of the
  // first file laser, else 0)
  double lasyTref = 0.0;
 public:
  bool lasyDue(int step, bool last) const {
    return cfg.output.lasy && ((cfg.output.lasy_every > 0 && step % cfg.output.lasy_every == 0) || last);
  }
  // ---- output of the window, block by block (no rank ever needs the whole window in memory)
  // A step's output files are opened by beginOutput, filled by outputBlock for the blocks [n0, n1) in increasing
  // order (fold must hold the slices [n0, n1], c0b / c0s the half-cells [n0, n1] (n1 < Nz) - one slice / half-cell
  // beyond the block for the interpolations), and closed by endOutput, which returns the console messages.
  // Without MPI the whole window is one block; the MPI pipeline streams the blocks of the ranks through rank 0.
  struct OutputPlan {
    bool files = false;       // the formats of output.format (text, hdf5, openpmd)
    bool lasy = false;        // output.lasy
    bool checkpoint = false;  // checkpoint_NNNNN.h5 (beam, pendingNext and rng must be complete)
    bool any() const { return files || lasy || checkpoint; }
  };
  struct Sink {
    virtual ~Sink() = default;
    virtual void block(int n0, int n1) = 0;
    virtual std::string finish() = 0;  // message for the console ("" for none)
  };
  struct OutputSet {
    std::vector<std::unique_ptr<Sink>> sinks;
  };
  std::unique_ptr<OutputSet> beginOutput(const OutputPlan& plan);
  void outputBlock(OutputSet& o, int n0, int n1);
  std::vector<std::string> endOutput(OutputSet& o);
  bool externalOutput = false;  // MPI pipeline: report() writes no files, rank 0 streams the blocks itself
 private:
  struct TextSink;        // diagnostics.cpp: axis_, map_, beam_NNNNN.txt
  struct Hdf5Sink;        // output_hdf5.cpp: data_NNNNN.h5
  struct OpenPMDSink;     // output_hdf5.cpp: openpmd/openpmd_<step>.h5
  struct LasySink;        // output_hdf5.cpp: lasy/laser_NNNNN.h5
  struct CheckpointSink;  // output_hdf5.cpp: checkpoint_NNNNN.h5
  static std::unique_ptr<Sink> makeHdf5Sink(Simulation& S, int kind);  // output_hdf5.cpp (0 data, 1 openPMD, 2 LASY,
                                                                       // 3 checkpoint)
  // the bunches of the deck (only < 0: all with z_start = 0), or only bunch `only`
  void initBunches(int only = -1);
  void initBunchFields();
  void addBeamFields(double sign, bool deckGamma, int only = -1);  // only >= 0: the field of that bunch only
  void loadBunches(double za, double zb);  // bunches with z_start in (za, zb]: loaded with their self-field
  void addAzimuthalFields(const Currents& C, double sign, double beta);  // co-moving self-fields of the beam list
  void applyChicanes(double za, double zb);          // thin chicanes passed in (za, zb]

  EB gather(const Fields& F, int n1, double w1, int n2, double w2, double x, double y) const;
  void accumulate(const Fields& F, int m, int n, double w, double r, const RadialGrid::W& wc,
                  const RadialGrid::W& wf, cd& er, cd& ep, cd& ez, cd& br, cd& bp, cd& bz) const;
  // kc >= 0: store half-cell k at index kc of C (a buffer holding fewer half-cells)
  void deposit(Currents& C, int k, double wk, double x, double y, double vx, double vy, double vz, double W,
               bool ion = false, int kc = -1);
  void depositBeamParticle(Currents& C, const Particle& p);
  // (1,2,1)/4 filter along zeta of the beam current (numerics.beam_filter passes):
  // removes the zeta-Nyquist component, which the scheme radiates resonantly
  // (it has omega = 0 there for every k_perp: numerical Cherenkov)
  void filterBeamCurrent(Currents& C) const;
  void smoothBeamRadially(Currents& C) const;
  // numerics.beam_consistent: charge-conserving radial flux (mode 0) of the beam list from the positions `start`
  // (index i of the list at s^0) to the positions of the first push (beamMap), into bflux.jr, filtered as the beam
  // current; ccStart: the first-push states of the last step (beam_field = mean)
  void beamFlux(const std::vector<Particle>& start);
  Currents bflux;
  std::vector<Particle> ccStart;

  void pushBeam();
  void march(double s1);
  void remeanBeam();
  void pushBeamParticle(Particle& p);
  // numerics.beam_field = mean: weight of F^1 in the second push, the beam at s^0, index map of the first push
  double wNew = 0.0;
  std::vector<Particle> beamSave;
  std::vector<int> beamMap;
  void loadFresh();
  void pushStreamParticle(Particle& p, int n, double s1);
  // plasma.smooth_length: radial smoothing of the streaming-plasma sources of half-cell n
  void smoothStream(Currents& C, int n);
  SourceFilter sfilter;
  mutable SourceFilter bfilter;  // numerics.beam_smooth_length: radial smoothing of the beam current
  // the same, depositing into half-cell kc of C and collecting freed electrons in bornOut
  void pushStreamParticle(Particle& p, int n, double s1, Currents& C, int kc, std::vector<Particle>& bornOut);
  // all streaming particles through half-cell n (threads > 1: OpenMP, per-thread buffers)
  void pushStreamSlice(int n, double s1, int threads);
  void ionize(Particle& p, double E1, double E0, double dt, std::vector<Particle>& bornOut);
  double fieldBound(const Fields& F, int n) const;
  void absorber(double s);  // at level s (flying focus: the incident field of that level)

  DiagSums diagSums(int n0, int n1) const;          // slices / half-cells [n0, n1) and the beam list
  // the same for a given beam list and freed charge (the MPI pipeline computes them while the next step runs)
  DiagSums diagSums(int n0, int n1, const std::vector<Particle>& beamList, const std::vector<double>& freedQ) const;
  void report(const DiagSums& D, bool files);       // console, tables, files
  void readCheckpoint(const std::string& path);
  void trimTables();  // restart: remove lines of later steps from bunches.txt, ionization.txt
  void writeBunchTable(const DiagSums& D);
  void writeFelTable(const DiagSums& D);
  void writeFelHarmonics(const DiagSums& D);
  void writeIonizationTable(const DiagSums& D);
};

}  // namespace gem
