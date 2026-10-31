// ============================================================================
//  pipeline.cpp -- MPI parallelization by domain decomposition along zeta
//
//  Rank q owns the slices and half-cells [lo_q, lo_{q+1}) (rank P-1 is at the
//  front of the window, rank 0 at the rear). Within one step s^0 -> s^1 the
//  zeta-march is sequential from the front to the rear, so the ranks work as a
//  pipeline in the steps: rank q computes step k while rank q+1 computes step
//  k+1. Per step and rank (a = lo_q, b = lo_{q+1}):
//
//   front -> rear, end of the step ("F" message, rank q+1 -> q):
//     slice b of the new level (needed to solve slice b-1), the prediction
//     predR(b) used there (the correction of R at slice b needs it), the beam
//     current that rank q+1 deposited into half-cells b-1, b-2, the streaming
//     particles crossing zeta_b, and beam particles that drift into rank q's
//     part of the window (with the sub-step at which their push continues).
//
//   rear -> front, after the first 3 slices ("R" message, rank q -> q+1):
//     the corrected R of slice b, the final slices b-2, b-1 and currents of
//     half-cell b-1 (ghosts for the next step), the current that rank q's
//     trapped particles deposited into half-cell b, and copies of rank q's beam
//     particles in the half-cell b-1 ("halo"): rank q+1 pushes these copies as
//     well, so that its prediction of j^1 in half-cell b-1 is complete.
//
//  A beam particle belongs to the rank of the upper half-cell of its
//  deposition (the slice nearest to it, ZetaView::nearest). Rank 0 finishes every step last; it
//  collects the diagnostics of all ranks and writes the output and the
//  checkpoints (same files as a serial run). Storage: every rank keeps only the slices and
//  half-cells [a - G, b + G) (G = beam_filter + 4); rank 0 assembles the whole window for output steps.
//  With one rank the pipelined step is identical to Simulation::step().
//
//  numerics.beam_field = mean (Simulation::remeanBeam): the second push of the beam of step k, in
//  (F^0 + F^1)/2, needs F^1 below the block, which the rank behind computes later. It is done at the start
//  of step k+1 ("mean phase", meanPhase), when the R message of step k has brought F^1 below the block
//  and the F message of step k+1 the results of the mean phase of the rank in front. The mean phase runs
//  from the front to the rear like the first push: particles below lowZ continue on the rank behind, the
//  halo copies complete the current below the block, and the raw current of the rank in front completes
//  the filter. The beam particles of s^0 carry a tag (Particle::id, unused on the beam list); the rank where
//  the second push of a particle ends decides with its first-push outcome, as Simulation::remeanBeam. The
//  output of step k follows its mean phase; after the last step one more exchange does its mean phase.
// ============================================================================
#include "pipeline.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include "comm.hpp"
#if defined(_OPENMP)
#include <omp.h>
#endif
#include "rng.hpp"
#include "simulation.hpp"

namespace gem {

namespace {
enum { TAG_F = 101, TAG_R = 102, TAG_DIAG = 103, TAG_BLOCK = 104 };
constexpr int MIN_SLICES = 4;  // per rank without beam filter (the R message needs 3 processed slices)
// beam_field = mean: tag of a beam particle of s^0 in Particle::id (streaming ids stay below 2^63)
constexpr unsigned long long TAG_BIT = 1ULL << 63;
inline bool tagged(const Particle& p) { return (p.id & TAG_BIT) != 0; }
// first-push outcome of a tag (status of a record in fpKnow)
enum { FP_BELOW = 0, FP_LOST = 2, FP_RESTREAMED = 3, FP_UNKNOWN = 4, FP_CARRIED = 5 };
}  // namespace

class PipelineRun {
 public:
  PipelineRun(Simulation& s, const RunControl& c) : S(s), rc(c) {}
  int run();

 private:
  Simulation& S;
  RunControl rc;
  int P = 1, r = 0, a = 0, b = 0;
  int sa = 0, sb = 0;  // stored slices [sa, sb] and half-cells [sa, sb) (block and ghosts, clipped)
  std::vector<int> lo;
  int lastStep = 0;
  bool wallStop = false;
  // particles from the neighbours for the current step
  std::vector<Particle> halo, inContinue, inFinished, inTransfer;
  std::vector<cd> predRfront, predPfront;  // prediction of the front rank at slice b
  Currents H;                              // currents of the halo copies (half-cells below a)
  Currents Tr;                             // trapped deposits of the front rank into b-1, b-2
  Currents Rw;                             // raw beam current of the front rank, half-cells b ... b+pf-1
  Currents Snap;                           // own raw beam current, half-cells a-2 ... a+pf-1
  std::vector<cd> predBelowR, predBelowP;  // filtered beam current of half-cell a-1 (prediction at slice a)
  int pf = 0;                              // passes of the beam-current filter
  int threads = 1;                         // OpenMP threads per rank (numerics.backend = openmp)
  long long myRestreamed = 0, myTrapped = 0;
  long baseRestreamed = 0, baseTrapped = 0;
  std::chrono::steady_clock::time_point t0;
  int stepsDone = 0;

  // ---- numerics.beam_field = mean (see the header): state kept from step k for its mean phase
  bool mean = false;
  std::vector<Particle> beamSave;          // own beam list at s^0 (tagged)
  std::vector<Particle> haloSave, haloFp;  // halo copies at s^0 and their first push (status FP_UNKNOWN: continued)
  std::unordered_map<unsigned long long, Particle> fpKnow;  // first-push outcomes by tag (status FP_*), except
                                                            // those on the own beam list
  std::vector<Particle> fpBackOut;  // first-push outcomes of the particles continued here from the front rank
  std::vector<Particle> trapped;    // particles trapped here during the march, in order
  Buffer trStash;                   // trapped deposits of the front rank into b-1, b-2
  Buffer rearT;                     // trapped deposits of the rank behind into a-1, a
  Currents Tq;                      // scratch: trapped deposits of this rank
  Buffer meanOut;                   // results of the mean phase for the rank behind (sent with the F message)
  bool meanReady = false;           // meanOut holds the results of a mean phase
  bool frontMean = false;           // the received F message has a mean section
  long long meanGuessed = 0;  // second pushes decided without the first-push outcome (it ended two ranks behind)
  // F^1 of the absorber layer before the radial absorber (the second push uses it, as Simulation::step)
  int iF = 0, iC = 0;               // first absorbed face / centre
  std::vector<std::vector<cd>> preAbs;  // per stored slice n - sa: R_r of faces >= iF, R_phi of centres >= iC
  std::vector<char> preHave;
  void absorb(int n, double s);  // absorbSlice of fnew, keeping the values before it (mean)
  void putPreAbs(Buffer& B, int n0, int n1) const;
  void getPreAbs(Buffer& B, int n0, int n1);
  void swapPreAbs();  // fold <-> values before the absorber (twice: unchanged)
  // the beam push (Simulation::pushBeamParticle) from sub-step p.sub; F1 != null: the second push in
  // (F0 + F1)/2. Returns true if the push must continue on the rank behind (zeta < lowZ).
  int nsub = 1;
  double pdz = 0, lowZ = -1e300;
  void beamPushSetup();
  bool pushBeamParticle(Particle& p, const Fields& F0, const Fields* F1, bool count);
  Buffer recvFront(int k);
  void meanPhase(Buffer* F);
  // numerics.backend = cuda (or GEMPIC_PIPELINE_ENGINE=1, testing): the march of the block on the device
  std::unique_ptr<Engine> eng;
  void deviceMarch(int k, double s1, std::vector<Particle>& newForFront, std::vector<Particle>& transferDown);
  void deviceBeamCurrent(Engine::Curr C, const std::vector<Particle>& depB, const std::vector<Particle>& depH,
                         bool meanPhase);

  int kh(const Particle& p) const { return S.zv.nearest(p.zeta); }
  // half-cells of the rank behind whose beam current this rank needs (filter reach)
  bool inHalo(int k) const { return r > 0 && k >= a - pf - 1 && k <= a - 1; }
  int ownerHalf(int k) const {
    k = std::max(0, std::min(S.Nz - 1, k));
    return int(std::upper_bound(lo.begin(), lo.end(), k) - lo.begin()) - 1;
  }
  int beamOwner(const Particle& p) const { return ownerHalf(kh(p)); }
  int joinOwner(const Particle& p) const { return ownerHalf(S.zv.joinCell(p.zeta)); }
  int firstOut = -1;  // step of the initial output (not after a restart)
  bool filesAt(int k) const {
    const int e = S.cfg.output.every;
    return (e > 0 && k % e == 0) || k == lastStep || k == firstOut;
  }
  bool lasyAt(int k) const { return S.lasyDue(k, k == lastStep); }
  bool checkpointAt(int k) const {
    const int c = S.cfg.output.checkpoint_every;
    return (c > 0 && (k % c == 0 || k == lastStep)) || (wallStop && k == lastStep);
  }

  void partition();
  void step(Buffer& F);
  void filterRaw(Currents& C, bool belowIntoC);
  void filterBeam();
  void sendFront(int k, const std::vector<Particle>& newForFront);
  void recvRear();
  void sendDiag(int k);
  void reportStep(int k);
  bool finiteOwn() const;
  void absorbSlice(Fields& F, int n, double s) const;

  // ---- packing of slices / half-cells
  void putSlice(Buffer& B, const Fields& F, int n) const {
    for (const auto* v : {&F.Rr, &F.Lr, &F.Bz})
      for (int m = 0; m <= S.M; ++m) B.putRaw(&(*v)[F.f(m, n, 0)], sizeof(cd) * (S.g.N + 1));
    for (const auto* v : {&F.Rp, &F.Lp, &F.Ez})
      for (int m = 0; m <= S.M; ++m) B.putRaw(&(*v)[F.c(m, n, 0)], sizeof(cd) * S.g.N);
  }
  void getSlice(Buffer& B, Fields& F, int n) {
    for (auto* v : {&F.Rr, &F.Lr, &F.Bz})
      for (int m = 0; m <= S.M; ++m) B.getRaw(&(*v)[F.f(m, n, 0)], sizeof(cd) * (S.g.N + 1));
    for (auto* v : {&F.Rp, &F.Lp, &F.Ez})
      for (int m = 0; m <= S.M; ++m) B.getRaw(&(*v)[F.c(m, n, 0)], sizeof(cd) * S.g.N);
  }
  void putSliceR(Buffer& B, const Fields& F, int n) const {
    for (int m = 0; m <= S.M; ++m) B.putRaw(&F.Rr[F.f(m, n, 0)], sizeof(cd) * (S.g.N + 1));
    for (int m = 0; m <= S.M; ++m) B.putRaw(&F.Rp[F.c(m, n, 0)], sizeof(cd) * S.g.N);
  }
  void getSliceR(Buffer& B, Fields& F, int n) {
    for (int m = 0; m <= S.M; ++m) B.getRaw(&F.Rr[F.f(m, n, 0)], sizeof(cd) * (S.g.N + 1));
    for (int m = 0; m <= S.M; ++m) B.getRaw(&F.Rp[F.c(m, n, 0)], sizeof(cd) * S.g.N);
  }
  void putHalf(Buffer& B, const Currents& C, int k) const {
    const int N = S.g.N;
    for (int m = 0; m <= S.M; ++m) B.putRaw(&C.jr[C.f(m, k, 0)], sizeof(cd) * (N + 1));
    for (int m = 0; m <= S.M; ++m) B.putRaw(&C.jp[C.c(m, k, 0)], sizeof(cd) * N);
    for (int m = 0; m <= S.M; ++m) B.putRaw(&C.jz[C.c(m, k, 0)], sizeof(cd) * N);
    B.putRaw(&C.rho[C.d(k, 0)], sizeof(double) * N);
    B.putRaw(&C.rhoi[C.d(k, 0)], sizeof(double) * N);
  }
  void getHalf(Buffer& B, Currents& C, int k, bool add) {
    const int N = S.g.N;
    std::vector<cd> t(N + 1);
    std::vector<double> d(N);
    auto take = [&](cd* dst, int n) {
      B.getRaw(t.data(), sizeof(cd) * n);
      for (int i = 0; i < n; ++i) dst[i] = add ? dst[i] + t[i] : t[i];
    };
    auto taked = [&](double* dst) {
      B.getRaw(d.data(), sizeof(double) * N);
      for (int i = 0; i < N; ++i) dst[i] = add ? dst[i] + d[i] : d[i];
    };
    for (int m = 0; m <= S.M; ++m) take(&C.jr[C.f(m, k, 0)], N + 1);
    for (int m = 0; m <= S.M; ++m) take(&C.jp[C.c(m, k, 0)], N);
    for (int m = 0; m <= S.M; ++m) take(&C.jz[C.c(m, k, 0)], N);
    taked(&C.rho[C.d(k, 0)]);
    taked(&C.rhoi[C.d(k, 0)]);
  }
  // slices [n0, n1) of the fields, half-cells [n0, n1) of the currents
  void putRange(Buffer& B, const Fields& F, int n0, int n1) const {
    const size_t nf = size_t(n1 - n0) * (S.g.N + 1), nc = size_t(n1 - n0) * S.g.N;
    for (const auto* v : {&F.Rr, &F.Lr, &F.Bz})
      for (int m = 0; m <= S.M; ++m) B.putRaw(&(*v)[F.f(m, n0, 0)], sizeof(cd) * nf);
    for (const auto* v : {&F.Rp, &F.Lp, &F.Ez})
      for (int m = 0; m <= S.M; ++m) B.putRaw(&(*v)[F.c(m, n0, 0)], sizeof(cd) * nc);
  }
  void getRange(Buffer& B, Fields& F, int n0, int n1) {
    const size_t nf = size_t(n1 - n0) * (S.g.N + 1), nc = size_t(n1 - n0) * S.g.N;
    for (auto* v : {&F.Rr, &F.Lr, &F.Bz})
      for (int m = 0; m <= S.M; ++m) B.getRaw(&(*v)[F.f(m, n0, 0)], sizeof(cd) * nf);
    for (auto* v : {&F.Rp, &F.Lp, &F.Ez})
      for (int m = 0; m <= S.M; ++m) B.getRaw(&(*v)[F.c(m, n0, 0)], sizeof(cd) * nc);
  }
  void putDens(Buffer& B, const Currents& C, int k0, int k1) const {
    B.putRaw(&C.rho[C.d(k0, 0)], sizeof(double) * size_t(k1 - k0) * S.g.N);
    B.putRaw(&C.rhoi[C.d(k0, 0)], sizeof(double) * size_t(k1 - k0) * S.g.N);
  }
  void getDens(Buffer& B, Currents& C, int k0, int k1) {
    B.getRaw(&C.rho[C.d(k0, 0)], sizeof(double) * size_t(k1 - k0) * S.g.N);
    B.getRaw(&C.rhoi[C.d(k0, 0)], sizeof(double) * size_t(k1 - k0) * S.g.N);
  }
  void putJ(Buffer& B, const Currents& C, int k0, int k1) const {
    const size_t nf = size_t(k1 - k0) * (S.g.N + 1), nc = size_t(k1 - k0) * S.g.N;
    for (int m = 0; m <= S.M; ++m) B.putRaw(&C.jr[C.f(m, k0, 0)], sizeof(cd) * nf);
    for (int m = 0; m <= S.M; ++m) B.putRaw(&C.jp[C.c(m, k0, 0)], sizeof(cd) * nc);
    for (int m = 0; m <= S.M; ++m) B.putRaw(&C.jz[C.c(m, k0, 0)], sizeof(cd) * nc);
  }
  void getJ(Buffer& B, Currents& C, int k0, int k1) {
    const size_t nf = size_t(k1 - k0) * (S.g.N + 1), nc = size_t(k1 - k0) * S.g.N;
    for (int m = 0; m <= S.M; ++m) B.getRaw(&C.jr[C.f(m, k0, 0)], sizeof(cd) * nf);
    for (int m = 0; m <= S.M; ++m) B.getRaw(&C.jp[C.c(m, k0, 0)], sizeof(cd) * nc);
    for (int m = 0; m <= S.M; ++m) B.getRaw(&C.jz[C.c(m, k0, 0)], sizeof(cd) * nc);
  }
  static void putDiag(Buffer& B, const DiagSums& D) {
    B.put(D.maxEx); B.put(D.minEz); B.put(D.maxEz); B.put(D.nimax);
    B.put<unsigned long long>(D.moments.size());
    for (const auto& kv : D.moments) { B.put(kv.first); B.put(kv.second); }
    B.putVec(D.freed);
    B.put(D.Wrad); B.put(D.W0); B.put(D.Pmax); B.put(D.zPmax);
    B.putVec(D.pslices);
    B.putVec(D.rslices);
    B.putVec(D.hslices);
    B.putVec(D.hbunch);
    B.put<unsigned long long>(D.slices.size());
    for (const auto& kv : D.slices) { B.put(kv.first); B.putVec(kv.second); }
  }
  static DiagSums getDiag(Buffer& B) {
    DiagSums D;
    D.maxEx = B.get<double>(); D.minEz = B.get<double>(); D.maxEz = B.get<double>(); D.nimax = B.get<double>();
    const size_t n = size_t(B.get<unsigned long long>());
    for (size_t i = 0; i < n; ++i) {
      const int sp = B.get<int>();
      D.moments[sp] = B.get<Moments>();
    }
    D.freed = B.getVec<double>();
    D.Wrad = B.get<double>(); D.W0 = B.get<double>(); D.Pmax = B.get<double>(); D.zPmax = B.get<double>();
    D.pslices = B.getVec<double>();
    D.rslices = B.getVec<double>();
    D.hslices = B.getVec<double>();
    D.hbunch = B.getVec<double>();
    const size_t ns = size_t(B.get<unsigned long long>());
    for (size_t i = 0; i < ns; ++i) {
      const int sp = B.get<int>();
      D.slices[sp] = B.getVec<double>();
    }
    return D;
  }
};

// ----------------------------------------------------------------------------
// decomposition of the window and of the (global) initial state
// ----------------------------------------------------------------------------
StorageBlock pipelineBlock(const Config& cfg, int P, int r) {
  const int Nz = cfg.grid.Nz, pf = cfg.numerics.beam_filter;
  // ghost width: the beam push reaches down to half-cell a - pf - 2, the filter and the exchanges up to b + pf - 1;
  // the bunch fields of the start are built on the stored range, where the filtered current is valid pf away
  // from its ends
  const int G = 2 * pf + 6;
  const int minSlices = std::max(MIN_SLICES + pf, G);
  if (P > 1 && cfg.numerics.beam_smooth_length > 0.0)
    throw std::runtime_error("MPI: numerics.beam_smooth_length is not available with more than one rank yet");
  if (cfg.numerics.beam_consistent)
    throw std::runtime_error("MPI pipeline: numerics.beam_consistent is not available in the pipelined step yet");
  for (const auto& B : cfg.bunches)
    if (B.z_start > 0)
      throw std::runtime_error("MPI pipeline: bunch.z_start > 0 (bunch '" + B.name + "') is not available yet");
  for (int q = 0; q < P; ++q) {
    const int lq = int((long long)Nz * q / P), hq = int((long long)Nz * (q + 1) / P);
    if (hq - lq < minSlices && P > 1)
      throw std::runtime_error("MPI: " + std::to_string(P) + " ranks for " + std::to_string(Nz) +
                               " slices; every rank needs at least " + std::to_string(minSlices) +
                               " slices (use at most " + std::to_string(Nz / minSlices) + " ranks)");
  }
  StorageBlock B;
  B.P = P;
  B.r = r;
  B.a = int((long long)Nz * r / P);
  B.b = int((long long)Nz * (r + 1) / P);
  B.g = G;
  B.sa = B.a - G;
  B.sb = B.b + G;
  if (P == 1 || std::getenv("GEMPIC_FULL_STORAGE")) { B.sa = 0; B.sb = Nz; }  // the whole window (testing)
  return B;
}

void PipelineRun::partition() {
  P = comm::size();
  r = comm::rank();
  const int Nz = S.Nz;
  lo.resize(P + 1);
  for (int q = 0; q <= P; ++q) lo[q] = int((long long)Nz * q / P);
  pf = S.cfg.numerics.beam_filter;
  const StorageBlock blk = pipelineBlock(S.cfg, P, r);
  a = blk.a;
  b = blk.b;
  // beam list: own particles, and the copies of the particles of rank r-1 in half-cell a-1
  std::vector<Particle> own;
  halo.clear();
  for (const auto& p : S.beam) {
    if (beamOwner(p) == r) own.push_back(p);
    else if (inHalo(kh(p))) halo.push_back(p);
  }
  S.beam.swap(own);
  std::vector<Particle> pend;
  for (const auto& p : S.pendingNext)
    if (joinOwner(p) == r) pend.push_back(p);
  S.pendingNext.swap(pend);
  baseRestreamed = (r == 0) ? S.n_restreamed : 0;
  baseTrapped = (r == 0) ? S.n_trapped_total : 0;
  // distributed storage: every rank keeps (and has built) only the slices and half-cells [sa, sb) of its block
  // and ghosts; all current arrays of a rank share one layout (filterBeam indexes the halo / front arrays like c1b)
  const int M = S.M, N = S.g.N;
  sa = blk.sa;
  sb = blk.sb;
  if (S.fold.n0 != std::max(0, sa) || S.c0b.k0 != std::max(0, sa))
    throw std::runtime_error("MPI pipeline: the simulation was not built for the block of this rank");
  for (auto* C : {&H, &Tr, &Rw, &Snap}) C->initRange(M, Nz, N, sa, sb);
  mean = S.cfg.numerics.beam_field == BEAM_FIELD_MEAN;
  if (mean) {
    Tq.initRange(M, Nz, N, sa, sb);
    const double ra = S.cfg.grid.r_absorber;
    iF = S.g.N + 1;
    iC = S.g.N;
    for (int i = S.g.N; i >= 0 && S.g.rho[i] > ra; --i) iF = i;
    for (int i = S.g.N - 1; i >= 0 && S.g.c[i] > ra; --i) iC = i;
    preAbs.assign(size_t(sb - sa + 1), {});
    preHave.assign(size_t(sb - sa + 1), 0);
  }
  predBelowR.assign(S.predR.size(), cd(0.0));
  predBelowP.assign(S.predP.size(), cd(0.0));
  predRfront.assign(S.predR.size(), cd(0.0));
  predPfront.assign(S.predP.size(), cd(0.0));
  // number of the last step (same arithmetic as the serial loop)
  lastStep = S.step_no + S.cfg.grid.stepsFrom(S.s0);
  // the march on the device: backend = cuda; GEMPIC_PIPELINE_ENGINE=1 also with backend = openmp (the GPU
  // emulation build runs the GPU code paths on CPU threads: testing)
  const char* pe = std::getenv("GEMPIC_PIPELINE_ENGINE");
  if (S.cfg.numerics.backend == NumericsConfig::Backend::Cuda ||
      (S.cfg.numerics.backend == NumericsConfig::Backend::OpenMP && pe && pe[0] == '1'))
    eng = std::make_unique<Engine>(S);
}

void PipelineRun::absorbSlice(Fields& F, int n, double s) const {
  const double ra = S.cfg.grid.r_absorber, rm = S.g.rho[S.g.N];
  if (ra >= rm) return;
  auto damp = [&](double rr) {
    if (rr <= ra) return 1.0;
    const double x = (rr - ra) / (rm - ra);
    return std::max(0.0, 1.0 - S.cfg.grid.absorber_strength * x * x);
  };
  const bool ff = S.cfg.flying_focus.present;  // as Simulation::absorber
  for (int m = 0; m <= S.M; ++m) {
    for (int i = 0; i <= S.g.N; ++i) {
      const double d = damp(S.g.rho[i]);
      if (d == 1.0) continue;
      const cd inc = ((ff && m == 1) ? cd(S.ffInc(true, n, i, s)) : cd(0.0)) + S.incidentF(m, n, i);
      cd& R = F.Rr[F.f(m, n, i)];
      R = inc + d * (R - inc);
    }
    for (int i = 0; i < S.g.N; ++i) {
      const double d = damp(S.g.c[i]);
      if (d == 1.0) continue;
      const cd inc = ((ff && m == 1) ? cd(0.0, 1.0) * S.ffInc(false, n, i, s) : cd(0.0)) + S.incidentC(m, n, i);
      cd& R = F.Rp[F.c(m, n, i)];
      R = inc + d * (R - inc);
    }
  }
}

void PipelineRun::absorb(int n, double s) {
  if (mean && iF + iC <= 2 * S.g.N) {
    const int nf = S.g.N + 1 - iF, nc = S.g.N - iC;
    std::vector<cd>& v = preAbs[size_t(n - sa)];
    v.resize(size_t(S.M + 1) * (nf + nc));
    for (int m = 0; m <= S.M; ++m) {
      std::copy_n(&S.fnew.Rr[S.fnew.f(m, n, iF)], nf, &v[size_t(m) * (nf + nc)]);
      std::copy_n(&S.fnew.Rp[S.fnew.c(m, n, iC)], nc, &v[size_t(m) * (nf + nc) + nf]);
    }
    preHave[size_t(n - sa)] = 1;
  }
  absorbSlice(S.fnew, n, s);
}

void PipelineRun::putPreAbs(Buffer& B, int n0, int n1) const {
  for (int n = n0; n <= n1; ++n) {
    const bool h = preHave[size_t(n - sa)] != 0;
    B.put<int>(h ? 1 : 0);
    if (h) B.putVec(preAbs[size_t(n - sa)]);
  }
}

void PipelineRun::getPreAbs(Buffer& B, int n0, int n1) {
  for (int n = n0; n <= n1; ++n) {
    const bool h = B.get<int>() != 0;
    preHave[size_t(n - sa)] = h ? 1 : 0;
    if (h) preAbs[size_t(n - sa)] = B.getVec<cd>();
  }
}

void PipelineRun::swapPreAbs() {
  const int nf = S.g.N + 1 - iF, nc = S.g.N - iC;
  for (int n = std::max(0, sa); n <= std::min(sb, S.Nz); ++n) {
    if (!preHave[size_t(n - sa)]) continue;
    std::vector<cd>& v = preAbs[size_t(n - sa)];
    for (int m = 0; m <= S.M; ++m) {
      std::swap_ranges(&v[size_t(m) * (nf + nc)], &v[size_t(m) * (nf + nc)] + nf, &S.fold.Rr[S.fold.f(m, n, iF)]);
      std::swap_ranges(&v[size_t(m) * (nf + nc) + nf], &v[size_t(m) * (nf + nc) + nf] + nc,
                       &S.fold.Rp[S.fold.c(m, n, iC)]);
    }
  }
}

bool PipelineRun::finiteOwn() const {
  for (int m = 0; m <= S.M; ++m)
    for (int n = a; n < b; ++n) {
      for (int i = 0; i <= S.g.N; ++i) {
        const cd v = S.fold.Rr[S.fold.f(m, n, i)];
        if (!std::isfinite(v.real()) || !std::isfinite(v.imag())) return false;
      }
      for (int i = 0; i < S.g.N; ++i) {
        const cd v = S.fold.Ez[S.fold.c(m, n, i)];
        if (!std::isfinite(v.real()) || !std::isfinite(v.imag())) return false;
      }
    }
  return true;
}

// ----------------------------------------------------------------------------
// beam push of one particle (Simulation::pushBeamParticle)
// ----------------------------------------------------------------------------
void PipelineRun::beamPushSetup() {
  nsub = std::max(1, int(std::ceil(S.Delta / S.cfg.numerics.beam_dz_max)));
  pdz = S.Delta / nsub;
  // lowest zeta whose current this rank needs: centre of half-cell a-pf-2
  lowZ = (r > 0 && a - pf - 2 >= -1) ? S.zv.centre(a - pf - 2) : -1e300;
}

bool PipelineRun::pushBeamParticle(Particle& p, const Fields& F0, const Fields* F1, bool count) {
  const auto& num = S.cfg.numerics;
  const double rmax = S.rWall;  // r_max, or the wall
  for (int sub = p.sub; sub < nsub && p.status == 0; ++sub) {
    if (p.zeta < lowZ) { p.sub = sub; return true; }
    double gm = gammaOf(p);
    const double vz = p.uz / gm;
    if (vz < num.vz_min) { p.status = 3; if (count) ++myRestreamed; break; }
    const double dt = pdz / vz;
    double w;
    const int na = S.zv.cell(p.zeta, w);
    if (na < 0 || na >= S.Nz) { p.status = 2; break; }
    EB f = S.gather(F0, na, 1.0 - w, na + 1, w, p.x, p.y);
    if (F1) {  // numerics.beam_field = mean, second push: (F^0 + F^1)/2 (as Simulation::pushBeamParticle)
      const double wNew = 0.5;
      const EB f1 = S.gather(*F1, na, 1.0 - w, na + 1, w, p.x, p.y);
      const double bo = 1.0 - wNew;
      f = {bo * f.Ex + wNew * f1.Ex, bo * f.Ey + wNew * f1.Ey, bo * f.Ez + wNew * f1.Ez,
           bo * f.Bx + wNew * f1.Bx, bo * f.By + wNew * f1.By, bo * f.Bz + wNew * f1.Bz};
    }
    addExternal(S.ext, p.x, p.y, p.z, f, pdz);
    if (num.beam_vay) vay(p.ux, p.uy, p.uz, f, p.qm, dt);
    else boris(p.ux, p.uy, p.uz, f, p.qm, dt);
    gm = gammaOf(p);
    p.x += p.ux / gm * dt; p.y += p.uy / gm * dt;
    p.zeta += (p.uz / gm - 1.0) * dt;
    p.z += pdz;
    if (p.zeta < 0.0 || std::hypot(p.x, p.y) >= rmax) p.status = 2;
  }
  p.sub = 0;
  return false;
}

// F message of step k from the front rank: header (step, last step, wall-time stop); the rest is read by
// meanPhase and step
Buffer PipelineRun::recvFront(int k) {
  Buffer B = comm::recv(r + 1, TAG_F);
  if (B.get<int>() != k) throw std::runtime_error("MPI pipeline: step mismatch (F)");
  lastStep = std::min(lastStep, B.get<int>());
  wallStop = wallStop || B.get<int>() != 0;
  frontMean = mean && B.get<int>() != 0;  // the message has a mean section (not in the first step)
  return B;
}

// ----------------------------------------------------------------------------
// one step s^0 -> s^1 on the slices [a, b) (port of Simulation::step); F: the F message of the front rank
// (header and mean section read)
// ----------------------------------------------------------------------------
void PipelineRun::step(Buffer& F) {
  const int N = S.g.N, M = S.M, Nz = S.Nz;
  S.beginStep();
  S.extendWall(S.fold);  // metal wall moved outwards: the opened annulus (the device march does it on the device)
  const double Delta = S.Delta;
  const ZetaView& zv = S.zv;
  const double s1 = S.s0 + Delta;
  const int k = S.step_no + 1;
  const auto& num = S.cfg.numerics;
  S.prepareIncident(s1, a - 1, b + 1);  // inject = boundary: the layer field of the own slices
  S.c1b.zero();
  S.c1s.zero();
  std::fill(S.freedCharge.begin(), S.freedCharge.end(), 0.0);
  S.pending.swap(S.pendingNext);
  S.pendingNext.clear();
  if (r > 0)
    for (int k = std::max(0, a - pf - 2); k <= a - 1; ++k) H.zeroHalf(k);
  if (mean) {  // tags of the beam particles of s^0; the state for the mean phase of this step
    for (size_t i = 0; i < S.beam.size(); ++i) S.beam[i].id = TAG_BIT | ((unsigned long long)r << 40) | i;
    beamSave = S.beam;
    fpKnow.clear();
    fpBackOut.clear();
    trapped.clear();
    std::fill(preHave.begin(), preHave.end(), 0);
  }

  // ---- from the front rank: slice b, its prediction, currents, particles
  S.stream.clear();
  inContinue.clear(); inFinished.clear(); inTransfer.clear();
  if (r < P - 1) {
    Buffer& B = F;
    getSlice(B, S.fnew, b);
    predRfront = B.getVec<cd>();
    predPfront = B.getVec<cd>();
    getHalf(B, S.c1b, b - 1, true);  // raw beam deposits of the front rank (filtered below)
    getHalf(B, S.c1b, b - 2, true);
    getHalf(B, Tr, b - 1, false);     // its trapped deposits (added after the filter)
    getHalf(B, Tr, b - 2, false);
    if (mean) {  // kept for the mean phase of this step
      trStash = Buffer();
      putHalf(trStash, Tr, b - 1);
      putHalf(trStash, Tr, b - 2);
    }
    for (int k = b; k < b + pf; ++k) getHalf(B, Rw, k, false);  // its raw beam current above b
    S.stream = B.getVec<Particle>();
    inContinue = B.getVec<Particle>();
    inFinished = B.getVec<Particle>();
    inTransfer = B.getVec<Particle>();
  }

  // ---- beam push (fields of s^0; slices a-2 ... b are available)
  beamPushSetup();
  if (eng) eng->putFields(0, sa, sb + 1);  // level s^0 on the device: the own block and the ghosts
  // returns true if the push must continue on the rank behind (zeta < lowZ)
  auto push = [&](Particle& p, bool count) { return pushBeamParticle(p, S.fold, nullptr, count); };
  std::vector<Particle> newBeam, handCont, handFin, transferDown;
  // beam deposits into c1b and H (the device: collected, deposited there after the bookkeeping)
  std::vector<Particle> depB, depH;
  auto toC1b = [&](const Particle& p) { if (eng) depB.push_back(p); else S.depositBeamParticle(S.c1b, p); };
  auto toH = [&](const Particle& p) { if (eng) depH.push_back(p); else S.depositBeamParticle(H, p); };
  // beam_field = mean: the first-push outcome of a tag that does not stay on this rank's list (fpKnow), and
  // for particles continued here from the front rank also in the report to it (fpBackOut)
  auto record = [&](const Particle& p, int st, bool fromFront) {
    if (!mean) return;
    Particle o = p;
    o.status = st;
    fpKnow[p.id] = o;
    if (fromFront) fpBackOut.push_back(o);
  };
  auto afterPush = [&](Particle& p, bool fromFront) {
    if (p.status == 3) {
      Particle q = p; q.status = 0; S.pending.push_back(q);
      record(p, FP_RESTREAMED, fromFront);
      return;
    }
    if (p.status != 0) { record(p, FP_LOST, fromFront); return; }
    if (beamOwner(p) == r) {
      newBeam.push_back(p);
      if (mean && fromFront) { Particle o = p; o.status = FP_BELOW; fpBackOut.push_back(o); }
      return;
    }
    if (inHalo(kh(p))) toH(p);  // completes the current below the block
    handFin.push_back(p);
    record(p, FP_BELOW, fromFront);
  };
  auto continued = [&](const Particle& p, bool fromFront) {
    handCont.push_back(p);
    if (mean && fromFront) { Particle o = p; o.status = FP_UNKNOWN; fpBackOut.push_back(o); }
  };
  if (eng) {
    // on the device (numerics.backend = cuda): the own beam and the pushes continued from the front rank
    std::vector<Particle> L(S.beam);
    L.insert(L.end(), inContinue.begin(), inContinue.end());
    std::vector<char> cont;
    eng->pushBeams(L, lowZ, 0.0, cont);
    const size_t nOwn = S.beam.size();
    for (size_t i = 0; i < L.size(); ++i) {
      if (L[i].status == 3) ++myRestreamed;
      if (cont[i]) continued(L[i], i >= nOwn);
      else afterPush(L[i], i >= nOwn);
    }
  } else if (threads > 1) {
    // the push of every particle is independent; the bookkeeping follows in the serial order
    std::vector<char> cont(S.beam.size());
    long long restreamed = 0;
#if defined(_OPENMP)
#pragma omp parallel for num_threads(threads) schedule(static) reduction(+ : restreamed)
#endif
    for (long long i = 0; i < (long long)S.beam.size(); ++i) {
      Particle& p = S.beam[size_t(i)];
      cont[size_t(i)] = push(p, false) ? 1 : 0;
      if (p.status == 3) ++restreamed;
    }
    myRestreamed += restreamed;
    for (size_t i = 0; i < S.beam.size(); ++i) {
      if (cont[i]) continued(S.beam[i], false);
      else afterPush(S.beam[i], false);
    }
    for (auto& p : inContinue) {
      if (push(p, true)) continued(p, true);
      else afterPush(p, true);
    }
  } else {
    for (auto& p : S.beam) {
      if (push(p, true)) continued(p, false);
      else afterPush(p, false);
    }
    for (auto& p : inContinue) {
      if (push(p, true)) continued(p, true);
      else afterPush(p, true);
    }
  }
  for (const auto& p : newBeam) toC1b(p);
  for (const auto& p : inFinished) {  // pushed by a rank in front
    if (beamOwner(p) == r) { toC1b(p); newBeam.push_back(p); }
    else {
      if (inHalo(kh(p))) toH(p);
      handFin.push_back(p);
      record(p, FP_BELOW, false);
    }
  }
  for (const auto& p : inTransfer) {  // trapped elsewhere, already deposited
    if (beamOwner(p) == r) newBeam.push_back(p);
    else transferDown.push_back(p);
  }
  S.beam.swap(newBeam);
  if (mean) { haloSave = halo; haloFp.clear(); }
  {  // copies of rank r-1's particles near the boundary
    std::vector<Particle> hc(halo);
    std::vector<char> cont(hc.size(), 0);
    if (eng) eng->pushBeams(hc, lowZ, 0.0, cont);
    else
      for (size_t i = 0; i < hc.size(); ++i) cont[i] = push(hc[i], false) ? 1 : 0;
    for (size_t i = 0; i < hc.size(); ++i) {
      const Particle& q = hc[i];
      if (mean) { haloFp.push_back(q); if (cont[i]) haloFp.back().status = FP_UNKNOWN; }
      if (!cont[i] && q.status == 0 && inHalo(kh(q))) toH(q);
    }
  }
  if (eng) deviceBeamCurrent(Engine::C1B, depB, depH, false);
  else filterBeam();

  // ---- streaming particles: fresh plasma (front rank), then the waiting ones
  if (r == P - 1) S.loadFresh();
  for (const auto& p : S.pending) {
    Particle q = p;
    q.status = STATUS_WAITING;
    S.stream.push_back(q);
  }
  S.pending.clear();
  const uint64_t keyStep = ionKeyStep(num.seed, S.step_no);
  if (r < P - 1) { S.predRprev = predRfront; S.predPprev = predPfront; }
  else {
    std::fill(S.predRprev.begin(), S.predRprev.end(), 0.0);
    std::fill(S.predPprev.begin(), S.predPprev.end(), 0.0);
  }

  // ---- zeta march over the own slices (numerics.backend = cuda: on the device, deviceMarch)
  std::vector<Particle> newForFront;
  int done = 0;
  if (eng) deviceMarch(k, s1, newForFront, transferDown);
  else for (int n = b - 1; n >= a; --n) {
    for (int m = 0; m <= M; ++m) {
      for (int i = 0; i <= N; ++i) {
        const cd c1 = (n == a && r > 0) ? predBelowR[m * (N + 1) + i]
                      : (n > 0)         ? S.c1b.jr[S.c1b.f(m, n - 1, i)]
                                        : cd(0.0);
        S.predR[m * (N + 1) + i] = (n > 0) ? c1 + S.c0s.jr[S.c0s.f(m, n - 1, i)] : 0.0;
      }
      for (int i = 0; i < N; ++i) {
        const cd c1 = (n == a && r > 0) ? predBelowP[m * N + i]
                      : (n > 0)         ? S.c1b.jp[S.c1b.c(m, n - 1, i)]
                                        : cd(0.0);
        S.predP[m * N + i] = (n > 0) ? c1 + S.c0s.jp[S.c0s.c(m, n - 1, i)] : 0.0;
      }
    }
    if (num.iterations > 1) S.backup = S.stream;
    for (int it = 0; it < num.iterations; ++it) {
      if (it == 0) S.fnew.copySlice(S.fold, n);
      else S.stream = S.backup;
      S.born.clear();
      if (!S.cfg.ions.empty()) S.sliceFieldBound = std::max(S.fieldBound(S.fnew, n + 1), S.fieldBound(S.fnew, n));
      S.c1s.zeroHalf(n);
      S.ionKey = ionKeySlice(keyStep, n, it);
      if (threads > 1) {
        S.pushStreamSlice(n, s1, threads);
      } else {
        for (auto& p : S.stream) {
          if (p.status == STATUS_WAITING) {
            if (p.zeta > zv.node(n)) p.status = 0;
            else continue;
          }
          S.pushStreamParticle(p, n, s1);
        }
      }
      S.smoothStream(S.c1s, n);
      for (int m = 0; m <= M; ++m) {
        const int nm = (n > 0) ? n - 1 : 0;  // last slice: half-cell -1 := half-cell 0 (Simulation::step)
        const double wl = zv.wLo(n), wh = zv.wHi(n);
        const double th1 = 2.0 * S.cfg.numerics.theta, th0 = 2.0 - th1;  // 1, 1 for the trapezoid
        for (int i = 0; i <= N; ++i) {
          const cd a0 = S.c0b.jr[S.c0b.f(m, nm, i)] + S.c0s.jr[S.c0s.f(m, nm, i)];
          const cd pr = (n > 0) ? S.predR[m * (N + 1) + i] : S.c1b.jr[S.c1b.f(m, 0, i)] + S.c1s.jr[S.c1s.f(m, 0, i)];
          S.Jr[i] = th0 * (wl * a0 + wh * S.c0b.jr[S.c0b.f(m, n, i)] + wh * S.c0s.jr[S.c0s.f(m, n, i)]) +
                    th1 * (wl * pr + wh * S.c1b.jr[S.c1b.f(m, n, i)] + wh * S.c1s.jr[S.c1s.f(m, n, i)]);
        }
        for (int i = 0; i < N; ++i) {
          const cd a0 = S.c0b.jp[S.c0b.c(m, nm, i)] + S.c0s.jp[S.c0s.c(m, nm, i)];
          const cd pr = (n > 0) ? S.predP[m * N + i] : S.c1b.jp[S.c1b.c(m, 0, i)] + S.c1s.jp[S.c1s.c(m, 0, i)];
          S.Jp[i] = th0 * (wl * a0 + wh * S.c0b.jp[S.c0b.c(m, n, i)] + wh * S.c0s.jp[S.c0s.c(m, n, i)]) +
                    th1 * (wl * pr + wh * S.c1b.jp[S.c1b.c(m, n, i)] + wh * S.c1s.jp[S.c1s.c(m, n, i)]);
        }
        for (int i = 0; i <= N; ++i) S.jr1[i] = S.c1b.jr[S.c1b.f(m, n, i)] + S.c1s.jr[S.c1s.f(m, n, i)];
        for (int i = 0; i < N; ++i) {
          S.jp1[i] = S.c1b.jp[S.c1b.c(m, n, i)] + S.c1s.jp[S.c1s.c(m, n, i)];
          S.jz1[i] = S.c1b.jz[S.c1b.c(m, n, i)] + S.c1s.jz[S.c1s.c(m, n, i)];
        }
        Fields& fo = S.fold;
        Fields& fn = S.fnew;
        S.solver.solveSlice(m, n, &fo.Rr[fo.f(m, n, 0)], &fo.Lr[fo.f(m, n, 0)], &fo.Rp[fo.c(m, n, 0)], &fo.Ez[fo.c(m, n, 0)],
                            &fo.Bz[fo.f(m, n, 0)], &fn.Rr[fn.f(m, n + 1, 0)], &fn.Lr[fn.f(m, n + 1, 0)],
                            &fn.Rp[fn.c(m, n + 1, 0)], &fn.Lp[fn.c(m, n + 1, 0)], &fn.Ez[fn.c(m, n + 1, 0)],
                            &fn.Bz[fn.f(m, n + 1, 0)], S.jr1.data(), S.jp1.data(), S.jz1.data(), S.Jr.data(),
                            S.Jp.data(), &fn.Rr[fn.f(m, n, 0)], &fn.Lr[fn.f(m, n, 0)], &fn.Rp[fn.c(m, n, 0)],
                            &fn.Lp[fn.c(m, n, 0)], &fn.Ez[fn.c(m, n, 0)], &fn.Bz[fn.f(m, n, 0)]);
      }
    }
    for (const auto& e : S.born) S.freedCharge[ionizedIndex(e.sp)] += e.w;
    S.stream.insert(S.stream.end(), S.born.begin(), S.born.end());
    S.born.clear();
    for (auto& p : S.stream) {
      if (p.status != 1) continue;
      Particle t = p;
      t.status = 0;
      t.z = s1;
      if (t.uz / gammaOf(t) >= num.vz_beam) {
        if (t.sp == SPECIES_PLASMA) t.sp = SPECIES_TRAPPED;
        else if (isIonizationElectron(t.sp)) t.sp = SPECIES_TRAPPED_IONIZED_BASE - ionizedIndex(t.sp);
        else if (isIon(t.sp)) t.sp = SPECIES_TRAPPED_ION;
        S.depositBeamParticle(S.c1b, t);
        ++myTrapped;
        if (mean) trapped.push_back(t);  // deposited again after the second push (Simulation::remeanBeam)
        const int o = beamOwner(t);
        if (o == r) S.beam.push_back(t);
        else if (o > r) newForFront.push_back(t);
        else transferDown.push_back(t);
      } else {
        S.pendingNext.push_back(t);
      }
    }
    S.stream.erase(std::remove_if(S.stream.begin(), S.stream.end(),
                                  [](const Particle& p) { return p.status != 0 && p.status != STATUS_WAITING; }),
                   S.stream.end());
    if (n + 1 <= Nz - 1) {
      const double q = 0.5 * S.cfg.numerics.theta * Delta * zv.wLo(n + 1);
      const int K = S.solver.wallK();  // inside the wall: faces <= K, centres < K
      for (int m = 0; m <= M; ++m) {
        for (int i = 1; i <= K; ++i) {
          cd jt = S.c1b.jr[S.c1b.f(m, n, i)] + S.c1s.jr[S.c1s.f(m, n, i)];
          S.fnew.Rr[S.fnew.f(m, n + 1, i)] -= q * (jt - S.predRprev[m * (N + 1) + i]);
        }
        for (int i = 0; i < K; ++i) {
          cd jt = S.c1b.jp[S.c1b.c(m, n, i)] + S.c1s.jp[S.c1s.c(m, n, i)];
          S.fnew.Rp[S.fnew.c(m, n + 1, i)] -= q * (jt - S.predPprev[m * N + i]);
        }
      }
      absorb(n + 1, s1);  // slice n+1 is final
    }
    S.predRprev.swap(S.predR);
    S.predPprev.swap(S.predP);
    if (++done == pf + 3 && r < P - 1) sendFront(k, newForFront);
  }
  if (r == 0 && !eng) absorb(0, s1);
  for (auto& v : S.exitStates) std::fill(v.begin(), v.end(), 0.0);
  if (r == 0)
    for (const auto& p : S.stream)
      if (p.status == 0 && isIon(p.sp) && p.sp != SPECIES_TRAPPED_ION) S.exitStates[ionIndex(p.sp)][p.Zi] += p.w;

  // ---- wall-time limit: decided by the front rank, passed down the pipeline
  ++stepsDone;
  if (r == P - 1 && S.cfg.output.walltime_hours > 0 && k < lastStep) {
    const double hrs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count() / 3600.0;
    if (hrs + (P + 1.5) * hrs / stepsDone > S.cfg.output.walltime_hours) {
      lastStep = k;
      wallStop = true;
    }
  }

  // ---- to the rear rank: slice a, prediction, currents, particles
  if (r > 0) {
    Buffer B;
    B.put<int>(k);
    B.put<int>(lastStep);
    B.put<int>(wallStop ? 1 : 0);
    if (mean) {  // results of the mean phase of the previous step
      B.put<int>(meanReady ? 1 : 0);
      if (meanReady) B.putRaw(meanOut.data.data(), meanOut.data.size());
    }
    putSlice(B, S.fnew, a);
    B.putVec(S.predRprev);  // = prediction made at slice a
    B.putVec(S.predPprev);
    putHalf(B, Snap, a - 1);  // raw beam deposits below the block
    putHalf(B, Snap, a - 2);
    putHalf(B, S.c1b, a - 1);  // deposits of the particles trapped here (c1b below a was reset after the filter)
    putHalf(B, S.c1b, a - 2);
    for (int kk = a; kk < a + pf; ++kk) putHalf(B, Snap, kk);  // raw beam current for the rear rank's filter
    B.putVec(S.stream);
    B.putVec(handCont);
    B.putVec(handFin);
    B.putVec(transferDown);
    comm::send(r - 1, TAG_F, std::move(B));
  }
  std::swap(S.fold, S.fnew);
  std::swap(S.c0b, S.c1b);
  std::swap(S.c0s, S.c1s);
  if (eng) eng->swapLevels();
  S.s0 = s1;
  ++S.step_no;
}

// numerics.backend = cuda: the beam current of s^1 on the device (filterBeam): the deposits of the own beam
// (depB) and of the particles below the block (depH, the halo H), the front rank's raw deposits below its
// block and its raw current above b (Rw); the raw current of [a-2, a+pf-1] for the rank behind (Snap) and
// the filtered current of half-cell a-1 (predBelowR / P) come back to the host. C = C1B (first push), or
// C0B in the mean phase (filterRaw with belowIntoC: half-cell a-1 keeps its filtered value)
void PipelineRun::deviceBeamCurrent(Engine::Curr C, const std::vector<Particle>& depB,
                                    const std::vector<Particle>& depH, bool meanPhase) {
  const int Nz = S.Nz;
  eng->zeroCurrents(C, sa, sb);
  if (r < P - 1) eng->putCurrents(C, b - 2, b);  // the front rank's raw deposits (in S.c1b / S.c0b)
  eng->depositBeams(depB, C);
  if (r > 0) {
    eng->zeroCurrents(Engine::AUX_H, sa, a);
    eng->depositBeams(depH, Engine::AUX_H);
  }
  const int k0 = (r > 0) ? a - 2 : a, k1 = std::min(b - 1, a + pf - 1);
  eng->getCurrents(C, k0, k1 + 1, &Snap);
  // the raw current over [klo, khi]: own + halo below a, the front rank's above b; then the filter
  const int klo = std::max(0, (r > 0) ? a - pf - 1 : a), khi = std::min(Nz - 1, b + pf - 1);
  if (r > 0) eng->combineCurrents(C, Engine::AUX_H, klo, a, true);
  if (b <= khi) {
    eng->putCurrents(Engine::AUX_RW, b, khi + 1, &Rw);
    eng->combineCurrents(C, Engine::AUX_RW, b, khi + 1, false);
  }
  eng->filterCurrents(C, klo, khi + 1, pf);
  if (r > 0 && !meanPhase) eng->getHalfJ(C, a - 1, predBelowR, predBelowP);
  // below a: from now on only the deposits of trapped particles; above b: nothing of this rank
  eng->zeroCurrents(C, klo, meanPhase ? a - 1 : a);
  eng->zeroCurrents(C, b, khi + 1);
  if (r < P - 1) eng->addCurrents(C, Tr, b - 2, b);  // trapped at the lowest slices of the front rank
}

// numerics.backend = cuda: the march of the own slices on the device (Engine::marchBlock). The host uploads
// the state of the step (level s^0 with the ghosts from the rank behind and the mean phase, slice b of s^1
// from the front rank, the currents of s^0 and the filtered beam current of s^1), sends the R message from
// device data after pf+3 slices, routes the trapped particles as in step(), and takes back level s^1.
void PipelineRun::deviceMarch(int k, double s1, std::vector<Particle>& newForFront,
                              std::vector<Particle>& transferDown) {
  const int N = S.g.N;
  // (level s^0 and the beam current c1b are on the device since the beam push)
  eng->putFields(1, b, b + 1);
  eng->putCurrents(Engine::C0B, sa, sb);
  eng->putCurrents(Engine::C0S, sa, sb);
  size_t routed = 0;
  auto route = [&]() {  // trapped particles (fast ones, already deposited into c1b on the device)
    for (const Particle& t : eng->trapped(routed)) {
      ++routed;
      ++myTrapped;
      if (mean) trapped.push_back(t);
      const int o = beamOwner(t);
      if (o == r) S.beam.push_back(t);
      else if (o > r) newForFront.push_back(t);
      else transferDown.push_back(t);
    }
  };
  auto fetchPreAbs = [&](int n0, int n1) {  // values of the absorber layer before the absorber (mean phase)
    for (int n = n0; n < n1; ++n)
      if (eng->preAbsorber(n, preAbs[size_t(n - sa)])) preHave[size_t(n - sa)] = 1;
  };
  Engine::BlockMarch bm;
  bm.a = a;
  bm.b = b;
  bm.s1 = s1;
  if (r > 0) { bm.predBelowR = &predBelowR; bm.predBelowP = &predBelowP; }
  if (r < P - 1) { bm.predFrontR = &predRfront; bm.predFrontP = &predPfront; }
  bm.absorbFirst = r == 0;
  bm.iF = mean ? iF : N + 1;
  bm.iC = mean ? iC : N;
  if (r < P - 1) {
    bm.callAfter = pf + 3;
    bm.after = [&]() {
      eng->getFields(1, b - pf - 2, b + 1);
      eng->getCurrents(Engine::C1B, b - 1, b + 1);
      eng->getCurrents(Engine::C1S, b - 1, b);
      if (mean) fetchPreAbs(b - pf - 2, b + 1);
      route();
      sendFront(k, newForFront);
    };
  }
  eng->marchBlock(bm);
  route();
  eng->getFields(1, a, b + 1);
  eng->getCurrents(Engine::C1B, sa, sb);
  eng->getCurrents(Engine::C1S, sa, sb);
  if (mean) fetchPreAbs(a, b + 1);
}

// Filter of the beam current (Simulation::filterBeamCurrent) on the own
// half-cells and on half-cell a-1 (prediction at slice a). The raw current is
// assembled over [a-pf-1, b+pf-1]: own deposits, the halo deposits below a,
// the front rank's raw current above b. Values within pf of the ends of this
// range are not needed. The raw own deposits below a and the raw current of
// [a, a+pf-1] are kept for the rear rank (Snap); c1b below a is reset and then
// collects only the deposits of particles trapped during the march.
void PipelineRun::filterBeam() {
  const int N = S.g.N, M = S.M;
  filterRaw(S.c1b, false);
  if (r > 0) {  // from now on: only the deposits of trapped particles below the block
    S.c1b.zeroHalf(a - 1);
    S.c1b.zeroHalf(a - 2);
  }
  // deposits of the particles trapped at the lowest slices of the front rank
  if (r < P - 1)
    for (int k = b - 2; k <= b - 1; ++k) {
      for (int m = 0; m <= M; ++m) {
        for (int i = 0; i <= N; ++i) S.c1b.jr[S.c1b.f(m, k, i)] += Tr.jr[Tr.f(m, k, i)];
        for (int i = 0; i < N; ++i) {
          S.c1b.jp[S.c1b.c(m, k, i)] += Tr.jp[Tr.c(m, k, i)];
          S.c1b.jz[S.c1b.c(m, k, i)] += Tr.jz[Tr.c(m, k, i)];
        }
      }
      for (int i = 0; i < N; ++i) {
        S.c1b.rho[S.c1b.d(k, i)] += Tr.rho[Tr.d(k, i)];
        S.c1b.rhoi[S.c1b.d(k, i)] += Tr.rhoi[Tr.d(k, i)];
      }
    }
}

// The filter itself on the raw current C (own deposits): the raw current of [a-2, a+pf-1] is kept in Snap;
// the filtered current of half-cell a-1 goes to predBelowR / predBelowP (the prediction at slice a), or with
// belowIntoC into C itself (all components: the mean phase builds the final current of s^1 there).
void PipelineRun::filterRaw(Currents& C, bool belowIntoC) {
  const int N = S.g.N, M = S.M, Nz = S.Nz;
  const int k0 = (r > 0) ? a - 2 : a, k1 = std::min(b - 1, a + pf - 1);
  for (int k = k0; k <= k1; ++k) {
    Snap.zeroHalf(k);
    for (int m = 0; m <= M; ++m) {
      for (int i = 0; i <= N; ++i) Snap.jr[Snap.f(m, k, i)] = C.jr[C.f(m, k, i)];
      for (int i = 0; i < N; ++i) {
        Snap.jp[Snap.c(m, k, i)] = C.jp[C.c(m, k, i)];
        Snap.jz[Snap.c(m, k, i)] = C.jz[C.c(m, k, i)];
      }
    }
    for (int i = 0; i < N; ++i) {
      Snap.rho[Snap.d(k, i)] = C.rho[C.d(k, i)];
      Snap.rhoi[Snap.d(k, i)] = C.rhoi[C.d(k, i)];
    }
  }
  const int klo = std::max(0, (r > 0) ? a - pf - 1 : a), khi = std::min(Nz - 1, b + pf - 1);
  std::vector<cd> col(size_t(khi - klo + 1));
  // raw current of half-cell k, entry i of mode m, for one of the arrays
  auto run = [&](std::vector<cd>& own, const std::vector<cd>& halo, const std::vector<cd>& front,
                 std::vector<cd>* below, int L, auto idx) {
    for (int m = 0; m <= M; ++m)
      for (int i = 0; i < L; ++i) {
        for (int k = klo; k <= khi; ++k) {
          const size_t j = idx(m, k, i);
          col[k - klo] = (k < a) ? own[j] + halo[j] : (k >= b) ? front[j] : own[j];
        }
        for (int p = 0; p < pf; ++p) {
          cd prev(0.0);
          for (int k = klo; k <= khi; ++k) {
            cd& x = col[k - klo];
            const cd cur = x;
            const cd next = (k + 1 <= khi) ? col[k + 1 - klo] : cd(0.0);
            x = 0.25 * prev + 0.5 * cur + 0.25 * next;
            prev = cur;
          }
        }
        for (int k = a; k < b; ++k) own[idx(m, k, i)] = col[k - klo];
        if (r > 0) {
          if (belowIntoC) own[idx(m, a - 1, i)] = col[a - 1 - klo];
          else if (below) (*below)[size_t(m) * L + i] = col[a - 1 - klo];
        }
      }
  };
  auto fi = [&](int m, int k, int i) { return C.f(m, k, i); };
  auto ci = [&](int m, int k, int i) { return C.c(m, k, i); };
  run(C.jr, H.jr, Rw.jr, &predBelowR, N + 1, fi);
  run(C.jp, H.jp, Rw.jp, &predBelowP, N, ci);
  run(C.jz, H.jz, Rw.jz, nullptr, N, ci);
}

// ----------------------------------------------------------------------------
// numerics.beam_field = mean: the second push of the step that just ended (fold = F^1, fnew = F^0, the beam
// list = the first-push results and the trapped particles; see the header). F: the F message of the next
// step with its mean section (null: the front rank). The results for the rank behind go to meanOut.
// ----------------------------------------------------------------------------
void PipelineRun::meanPhase(Buffer* F) {
  using u64 = unsigned long long;
  swapPreAbs();  // F^1 before the absorber, as Simulation::step pushes before absorber()
  beamPushSetup();
  const Fields& F0 = S.fnew;
  const Fields& F1 = S.fold;
  // numerics.backend = cuda: the pushes, deposits and the filter on the device; F^1 (before the absorber) to
  // the device level 0 (F^0 is there as level 1; level 0 is uploaded again, after the absorber, in step())
  if (eng) eng->putFields(0, sa, sb + 1);
  std::vector<Particle> depH;
  auto toH = [&](const Particle& p) { if (eng) depH.push_back(p); else S.depositBeamParticle(H, p); };
  if (r > 0)
    for (int k = std::max(0, a - pf - 2); k <= a - 1; ++k) H.zeroHalf(k);
  S.c0b.zero();  // the beam current of s^1, deposited anew
  std::vector<Particle> inFin, inCont, inCarry;
  std::vector<u64> inDrop;
  if (F && frontMean) {
    getHalf(*F, S.c0b, b - 1, true);  // raw deposits of the front rank's particles below its block
    getHalf(*F, S.c0b, b - 2, true);
    for (int kk = b; kk < b + pf; ++kk) getHalf(*F, Rw, kk, false);  // its raw current above b
    inFin = F->getVec<Particle>();
    inCont = F->getVec<Particle>();
    inCarry = F->getVec<Particle>();
    inDrop = F->getVec<u64>();
  }
  for (const auto& c : inCarry) fpKnow[c.id] = c;

  // the beam list: first-push results (tagged), particles trapped in the step (no tag)
  std::vector<Particle>& L = S.beam;
  std::unordered_map<u64, size_t> at;
  for (size_t i = 0; i < L.size(); ++i)
    if (tagged(L[i])) at[L[i].id] = i;
  std::vector<char> erased(L.size(), 0), fin(L.size(), 0);
  std::vector<Particle> outFin, outCont, outCarry, newHalo;
  std::vector<u64> outDrop;
  auto held = [&](u64 t) -> long long {
    const auto it = at.find(t);
    return (it != at.end() && !erased[it->second]) ? (long long)it->second : -1;
  };
  auto eraseHere = [&](u64 t) {
    const long long i = held(t);
    if (i < 0 || fin[size_t(i)]) return false;
    erased[size_t(i)] = 1;
    return true;
  };
  // the final state X of a particle: on the own list (in place of its first-push state), or towards its owner
  auto route = [&](Particle X) {
    X.status = 0;
    X.sub = 0;
    if (beamOwner(X) == r) {
      const long long i = held(X.id);
      if (i >= 0) { L[size_t(i)] = X; fin[size_t(i)] = 1; }
      else { at[X.id] = L.size(); L.push_back(X); erased.push_back(0); fin.push_back(1); }
      return;
    }
    eraseHere(X.id);
    outFin.push_back(X);
    if (inHalo(kh(X))) { toH(X); newHalo.push_back(X); }  // completes the current below
  };
  // Simulation::remeanBeam: the second-push state s if it survived, else the first-push state; nothing if the
  // first push lost the particle or returned it to streaming
  auto decide = [&](const Particle& s) {
    const u64 t = s.id;
    int where = FP_UNKNOWN;
    Particle fp{};
    if (held(t) >= 0) where = -1;  // on the own list
    else if (const auto k = fpKnow.find(t); k != fpKnow.end()) { fp = k->second; where = fp.status; }
    if (where == FP_LOST || where == FP_RESTREAMED) return;
    if (s.status == 0) {
      if (where == FP_BELOW || where == FP_UNKNOWN) outDrop.push_back(t);  // the first-push state is elsewhere
      if (where == FP_UNKNOWN) ++meanGuessed;  // (assumed to have survived)
      route(s);
    } else if (where == FP_CARRIED) {
      route(fp);
    }
  };
  // a second push that continues on the rank behind takes the first-push outcome along if it is known here
  auto carry = [&](const Particle& s) {
    const long long i = held(s.id);
    if (i >= 0) {
      Particle c = L[size_t(i)];
      c.status = FP_CARRIED;
      outCarry.push_back(c);
      erased[size_t(i)] = 1;
    } else if (const auto k = fpKnow.find(s.id); k != fpKnow.end() && k->second.status != FP_BELOW &&
                                                   k->second.status != FP_UNKNOWN) {
      outCarry.push_back(k->second);
    }
  };

  // ---- final states decided by the ranks in front
  for (const auto& X : inFin) route(X);
  // ---- second push of the own beam of s^0 and of the pushes continued from the front rank
  std::vector<Particle>* lists[2] = {&beamSave, &inCont};
  for (auto* list : lists) {
    std::vector<char> cont(list->size());
    if (eng) eng->pushBeams(*list, lowZ, 0.5, cont);
    else {
#if defined(_OPENMP)
#pragma omp parallel for num_threads(threads) schedule(static) if (threads > 1)
#endif
      for (long long i = 0; i < (long long)list->size(); ++i)
        cont[size_t(i)] = pushBeamParticle((*list)[size_t(i)], F0, &F1, false) ? 1 : 0;
    }
    for (size_t i = 0; i < list->size(); ++i) {
      const Particle& s = (*list)[i];
      if (cont[i]) { outCont.push_back(s); carry(s); }
      else decide(s);
    }
  }
  // ---- first-push states to remove (their particle was decided elsewhere)
  for (const u64 t : inDrop)
    if (!eraseHere(t) && r > 0) outDrop.push_back(t);
  // ---- the halo copies (particles of the rank behind near the boundary), as decide()
  std::vector<Particle> hs(haloSave);
  std::vector<char> hcont(hs.size(), 0);
  if (eng) eng->pushBeams(hs, lowZ, 0.5, hcont);
  else
    for (size_t i = 0; i < hs.size(); ++i) hcont[i] = pushBeamParticle(hs[i], F0, &F1, false) ? 1 : 0;
  for (size_t i = 0; i < hs.size(); ++i) {
    const Particle& s = hs[i];
    if (hcont[i]) continue;
    const Particle& fp = haloFp[i];
    if (fp.status == FP_LOST || fp.status == FP_RESTREAMED) continue;
    Particle X = (s.status == 0) ? s : fp;
    if (X.status != 0) continue;
    X.sub = 0;
    if (inHalo(kh(X))) { toH(X); newHalo.push_back(X); }
  }
  if (r > 0) halo.insert(halo.end(), newHalo.begin(), newHalo.end());  // + the trapped ones (R message)
  // ---- the new list
  {
    std::vector<Particle> nl;
    nl.reserve(L.size());
    for (size_t i = 0; i < L.size(); ++i)
      if (!erased[i]) nl.push_back(L[i]);
    L.swap(nl);
  }
  // ---- the beam current of s^1: the particles of the first push filtered, the trapped ones added unfiltered
  if (eng) {
    std::vector<Particle> depB;
    for (const auto& p : L)
      if (tagged(p)) depB.push_back(p);
    deviceBeamCurrent(Engine::C0B, depB, depH, true);  // (+ the front rank's trapped deposits, as trStash)
    eng->depositBeams(trapped, Engine::C0B);
    if (r > 0) {
      rearT.pos = 0;
      getHalf(rearT, Tq, a - 1, false);
      getHalf(rearT, Tq, a, false);
      eng->addCurrents(Engine::C0B, Tq, a - 1, a + 1);
    }
    eng->getCurrents(Engine::C0B, sa, sb);  // the host copy (output, the next step)
  } else {
    for (const auto& p : L)
      if (tagged(p)) S.depositBeamParticle(S.c0b, p);
    filterRaw(S.c0b, true);
    if (r < P - 1) {
      trStash.pos = 0;
      getHalf(trStash, S.c0b, b - 1, true);
      getHalf(trStash, S.c0b, b - 2, true);
    }
    for (const auto& t : trapped) S.depositBeamParticle(S.c0b, t);
    if (r > 0) {
      rearT.pos = 0;
      getHalf(rearT, S.c0b, a - 1, true);
      getHalf(rearT, S.c0b, a, true);
    }
  }
  swapPreAbs();  // back to F^1 after the absorber
  // ---- for the rank behind (with the F message of the next step)
  meanOut = Buffer();
  if (r > 0) {
    putHalf(meanOut, Snap, a - 1);  // raw deposits below the block
    putHalf(meanOut, Snap, a - 2);
    for (int kk = a; kk < a + pf; ++kk) putHalf(meanOut, Snap, kk);  // raw current for its filter
    meanOut.putVec(outFin);
    meanOut.putVec(outCont);
    meanOut.putVec(outCarry);
    meanOut.putVec(outDrop);
  }
  meanReady = true;
  fpKnow.clear();
}

// after the first pf+3 slices: what the front rank needs for its next step
void PipelineRun::sendFront(int k, const std::vector<Particle>& newForFront) {
  Buffer B;
  B.put<int>(k);
  for (int n = b - pf - 2; n <= b - 1; ++n) putSlice(B, S.fnew, n);
  putSliceR(B, S.fnew, b);
  putHalf(B, S.c1b, b - 1);
  putHalf(B, S.c1s, b - 1);
  putHalf(B, S.c1b, b);  // deposits of the particles trapped at slice b-1
  // halo copies for the front rank. beam_field = mean: only the particles trapped in this step (final); the
  // front rank pushes its own copies of s^0 a second time (meanPhase)
  std::vector<Particle> window;
  for (const auto& p : S.beam)
    if (kh(p) >= b - pf - 1 && kh(p) <= b - 1 && !(mean && tagged(p))) window.push_back(p);
  B.putVec(window);
  B.putVec(newForFront);
  if (mean) {
    putPreAbs(B, b - pf - 2, b);  // F^1 of the absorber layer before the absorber
    // deposits of the particles trapped here so far into half-cell b-1 (alone)
    Tq.zero();
    for (const auto& t : trapped) S.depositBeamParticle(Tq, t);
    putHalf(B, Tq, b - 1);
    B.putVec(fpBackOut);  // first-push outcomes of the particles the front rank continued here
  }
  comm::send(r + 1, TAG_R, std::move(B));
}

// applied after the swap of step k: ghosts and boundary values of level s^1
void PipelineRun::recvRear() {
  Buffer B = comm::recv(r - 1, TAG_R);
  if (B.get<int>() != S.step_no) throw std::runtime_error("MPI pipeline: step mismatch (R)");
  for (int n = a - pf - 2; n <= a - 1; ++n) getSlice(B, S.fold, n);
  getSliceR(B, S.fold, a);
  getHalf(B, S.c0b, a - 1, false);
  getHalf(B, S.c0s, a - 1, false);
  if (mean) getHalf(B, Tq, a, false);  // the mean phase builds c0b anew (rearT)
  else getHalf(B, S.c0b, a, true);
  halo = B.getVec<Particle>();
  const std::vector<Particle> nb = B.getVec<Particle>();
  S.beam.insert(S.beam.end(), nb.begin(), nb.end());
  if (mean) {
    getPreAbs(B, a - pf - 2, a);
    getHalf(B, Tq, a - 1, false);
    rearT = Buffer();
    putHalf(rearT, Tq, a - 1);
    putHalf(rearT, Tq, a);
    for (const auto& o : B.getVec<Particle>()) fpKnow[o.id] = o;  // outcomes of the particles continued there
  }
}

// ----------------------------------------------------------------------------
// diagnostics: partial sums (and data of output / checkpoint steps) to rank 0
// ----------------------------------------------------------------------------
void PipelineRun::sendDiag(int k) {
  const bool files = filesAt(k), ck = checkpointAt(k), fields = files || ck || lasyAt(k);
  Buffer B;
  B.put<int>(k);
  putDiag(B, S.diagSums(a, b));
  B.put(myRestreamed);
  B.put(myTrapped);
  if (files || ck) B.putVec(S.beam);
  if (ck) {
    B.putVec(S.pendingNext);
    std::string rs;
    if (r == P - 1) {
      std::ostringstream os;
      os << S.rng;
      rs = os.str();
    }
    B.putVec(std::vector<char>(rs.begin(), rs.end()));
  }
  comm::send(0, TAG_DIAG, std::move(B));
  if (fields) {  // the own block for the output, streamed through rank 0 (separate message)
    Buffer F;
    F.put<int>(k);
    putRange(F, S.fold, a, (r == P - 1) ? S.Nz + 1 : b);
    putJ(F, S.c0b, a, b);
    putJ(F, S.c0s, a, b);
    putDens(F, S.c0b, a, b);
    putDens(F, S.c0s, a, b);
    comm::send(0, TAG_BLOCK, std::move(F));
  }
}

// rank 0: console line and tables from the partial sums; the output files are written block by block: the block
// of rank q (fields of its slices, currents of its half-cells) with the first slice / half-cell of rank q+1
// appended, so that rank 0 never holds more than two blocks
void PipelineRun::reportStep(int k) {
  DiagSums D = S.diagSums(a, b);
  long long re = myRestreamed, tr = myTrapped;
  const bool files = filesAt(k), ck = checkpointAt(k), lasy = lasyAt(k), fields = files || ck || lasy;
  std::vector<Particle> gbeam, gpend;
  if (files || ck) gbeam = S.beam;
  if (ck) gpend = S.pendingNext;
  std::string rngState;
  for (int q = 1; q < P; ++q) {
    Buffer B = comm::recv(q, TAG_DIAG);
    if (B.get<int>() != k) throw std::runtime_error("MPI pipeline: step mismatch (diagnostics)");
    D.add(getDiag(B));
    re += B.get<long long>();
    tr += B.get<long long>();
    if (files || ck) {
      const std::vector<Particle> pb = B.getVec<Particle>();
      gbeam.insert(gbeam.end(), pb.begin(), pb.end());
    }
    if (ck) {
      const std::vector<Particle> pp = B.getVec<Particle>();
      gpend.insert(gpend.end(), pp.begin(), pp.end());
      const std::vector<char> rs = B.getVec<char>();
      if (q == P - 1) rngState.assign(rs.begin(), rs.end());
    }
  }
  S.n_restreamed = baseRestreamed + long(re);
  S.n_trapped_total = baseTrapped + long(tr);
  S.report(D, files);  // console and tables (S.externalOutput: no files)
  if (!fields) { std::fflush(stdout); return; }

  const int M = S.M, Nz = S.Nz, N = S.g.N;
  struct Blk { Fields F; Currents Cb, Cs; int lo = 0, hi = 0; };
  auto getBlock = [&](int q) {
    auto bl = std::make_unique<Blk>();
    bl->lo = lo[q];
    bl->hi = lo[q + 1];
    const int fe = (q == P - 1) ? Nz + 1 : bl->hi;
    bl->F.initRange(M, Nz, N, bl->lo, fe);
    bl->Cb.initRange(M, Nz, N, bl->lo, bl->hi);
    bl->Cs.initRange(M, Nz, N, bl->lo, bl->hi);
    if (q == 0) {
      bl->F.copyRange(S.fold, bl->lo, fe);
      bl->Cb.copyRange(S.c0b, bl->lo, bl->hi);
      bl->Cs.copyRange(S.c0s, bl->lo, bl->hi);
    } else {
      Buffer B = comm::recv(q, TAG_BLOCK);
      if (B.get<int>() != k) throw std::runtime_error("MPI pipeline: step mismatch (output block)");
      getRange(B, bl->F, bl->lo, fe);
      getJ(B, bl->Cb, bl->lo, bl->hi);
      getJ(B, bl->Cs, bl->lo, bl->hi);
      getDens(B, bl->Cb, bl->lo, bl->hi);
      getDens(B, bl->Cs, bl->lo, bl->hi);
    }
    return bl;
  };
  Simulation::OutputPlan plan;
  plan.files = files;
  plan.lasy = lasy;
  plan.checkpoint = ck;
  if (files || ck) S.beam.swap(gbeam);
  if (ck) {
    S.pendingNext.swap(gpend);
    if (P > 1) {
      std::istringstream is(rngState);
      is >> S.rng;  // rank 0 does not load plasma itself: its generator is unused
    }
  }
  {
    auto out = S.beginOutput(plan);
    std::unique_ptr<Blk> cur = getBlock(0);
    for (int q = 0; q < P; ++q) {
      std::unique_ptr<Blk> next = (q + 1 < P) ? getBlock(q + 1) : nullptr;
      // the block of rank q, extended by one slice / half-cell of rank q+1
      Fields F;
      Currents Cb, Cs;
      const int fe = (q == P - 1) ? Nz + 1 : cur->hi + 1;
      F.initRange(M, Nz, N, cur->lo, fe);
      F.copyRange(cur->F, cur->lo, fe);
      Cb.initRange(M, Nz, N, cur->lo, cur->hi + 1);
      Cs.initRange(M, Nz, N, cur->lo, cur->hi + 1);
      Cb.copyRange(cur->Cb, cur->lo, cur->hi);
      Cs.copyRange(cur->Cs, cur->lo, cur->hi);
      if (next) {
        F.copyRange(next->F, cur->hi, cur->hi + 1);
        Cb.copyRange(next->Cb, cur->hi, cur->hi + 1);
        Cs.copyRange(next->Cs, cur->hi, cur->hi + 1);
      }
      cur.reset();
      std::swap(S.fold, F);
      std::swap(S.c0b, Cb);
      std::swap(S.c0s, Cs);
      S.outputBlock(*out, lo[q], lo[q + 1]);
      std::swap(S.fold, F);
      std::swap(S.c0b, Cb);
      std::swap(S.c0s, Cs);
      cur = std::move(next);
    }
    for (const auto& msg : S.endOutput(*out)) std::printf("%s\n", msg.c_str());
  }
  if (files || ck) S.beam.swap(gbeam);
  if (ck) S.pendingNext.swap(gpend);
  std::fflush(stdout);
}

// ----------------------------------------------------------------------------
int PipelineRun::run() {
  P = comm::size();
#if defined(_OPENMP)
  if (S.cfg.numerics.backend == NumericsConfig::Backend::OpenMP) threads = std::max(1, omp_get_max_threads());
#endif
  if (!S.cfg.chicanes.empty()) {
    if (comm::rank() == 0)
      std::fprintf(stderr, "error: chicane blocks are not supported in MPI runs yet (the chicane moves the beam "
                           "across the domains of the ranks); run on one rank with numerics.backend = openmp\n");
    return 1;
  }
  r = comm::rank();
  S.externalOutput = true;  // rank 0 streams the output block by block (reportStep)
  t0 = std::chrono::steady_clock::now();
  try {
    partition();
    if (!rc.restarted) {  // initial state: collected from the blocks like every output step
      firstOut = S.step_no;
      if (r > 0) sendDiag(firstOut);
      else reportStep(firstOut);
    }
    if (r == 0) {
      std::printf("MPI pipeline along zeta: %d ranks, slices per rank %d ... %d\n", P, lo[1] - lo[0],
                  lo[P] - lo[P - 1]);
      if (rc.restarted) std::printf("continuing at step %d\n", S.step_no);
      std::fflush(stdout);
    }
    const int k0 = S.step_no;
    for (int k = k0 + 1; k <= lastStep; ++k) {
      if (k > k0 + 1 && r > 0) {
        recvRear();
        if (!mean) sendDiag(k - 1);
      }
      Buffer F;
      if (r < P - 1) F = recvFront(k);
      if (k > k0 + 1 && mean) {  // the second push of step k-1, then its output
        meanPhase(r < P - 1 ? &F : nullptr);
        if (r > 0) sendDiag(k - 1);
        else reportStep(k - 1);
      }
      step(F);
      if (!finiteOwn()) {
        std::fprintf(stderr, "non-finite field at step %d on rank %d\n", S.step_no, r);
        comm::abort(2);
      }
      if (r == 0 && !mean) reportStep(k);
    }
    if (lastStep > k0) {
      if (r > 0) recvRear();
      if (mean) {  // the second push of the last step: one more exchange from the front to the rear
        Buffer F;
        if (r < P - 1) F = recvFront(lastStep + 1);
        meanPhase(r < P - 1 ? &F : nullptr);
        if (r > 0) {
          Buffer B;
          B.put<int>(lastStep + 1);
          B.put<int>(lastStep);
          B.put<int>(wallStop ? 1 : 0);
          B.put<int>(1);
          B.putRaw(meanOut.data.data(), meanOut.data.size());
          comm::send(r - 1, TAG_F, std::move(B));
        }
      }
      if (r > 0) sendDiag(lastStep);
      else if (mean) reportStep(lastStep);
    }
    comm::waitAll();
    if (meanGuessed > 0)
      std::printf("rank %d, beam_field = mean: %lld second pushes ended where the first push of the particle was "
                  "not known (it had continued two ranks behind); taken as survived\n", r, meanGuessed);
  } catch (const std::exception& e) {
    std::fprintf(stderr, "error on rank %d: %s\n", r, e.what());
    std::fflush(stderr);
    comm::abort(2);
  }
  if (r == 0) {
    const double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    std::printf("done: %d steps in %.1f s, particles returned from the beam list to streaming: %ld\n", S.step_no,
                sec, S.n_restreamed);
    if (wallStop)
      std::printf("wall-time limit (%.3g h) reached at step %d; continue with:  mpirun ... gempic_rz %s --restart latest\n",
                  S.cfg.output.walltime_hours, S.step_no, rc.deckPath.c_str());
  }
  return wallStop ? 3 : 0;
}

int runPipeline(Simulation& sim, const RunControl& rc) {
  PipelineRun run(sim, rc);
  return run.run();
}

}  // namespace gem
