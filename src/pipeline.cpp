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
#include <vector>

#include "comm.hpp"
#if defined(_OPENMP)
#include <omp.h>
#endif
#include "rng.hpp"
#include "simulation.hpp"

namespace gem {

namespace {
enum { TAG_F = 101, TAG_R = 102, TAG_DIAG = 103 };
constexpr int MIN_SLICES = 4;  // per rank without beam filter (the R message needs 3 processed slices)
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

  int kh(const Particle& p) const { return S.zv.nearest(p.zeta); }
  // half-cells of the rank behind whose beam current this rank needs (filter reach)
  bool inHalo(int k) const { return r > 0 && k >= a - pf - 1 && k <= a - 1; }
  int ownerHalf(int k) const {
    k = std::max(0, std::min(S.Nz - 1, k));
    return int(std::upper_bound(lo.begin(), lo.end(), k) - lo.begin()) - 1;
  }
  int beamOwner(const Particle& p) const { return ownerHalf(kh(p)); }
  int joinOwner(const Particle& p) const { return ownerHalf(S.zv.joinCell(p.zeta)); }
  bool filesAt(int k) const {
    const int e = S.cfg.output.every;
    return (e > 0 && k % e == 0) || k == lastStep;
  }
  bool lasyAt(int k) const { return S.lasyDue(k, k == lastStep); }
  bool checkpointAt(int k) const {
    const int c = S.cfg.output.checkpoint_every;
    return (c > 0 && (k % c == 0 || k == lastStep)) || (wallStop && k == lastStep);
  }

  void partition();
  void step();
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
void PipelineRun::partition() {
  P = comm::size();
  r = comm::rank();
  const int Nz = S.Nz;
  lo.resize(P + 1);
  for (int q = 0; q <= P; ++q) lo[q] = int((long long)Nz * q / P);
  pf = S.cfg.numerics.beam_filter;
  const int minSlices = MIN_SLICES + pf;
  for (int q = 0; q < P; ++q)
    if (lo[q + 1] - lo[q] < minSlices && P > 1)
      throw std::runtime_error("MPI: " + std::to_string(P) + " ranks for " + std::to_string(Nz) +
                               " slices; every rank needs at least " + std::to_string(minSlices) +
                               " slices (use at most " + std::to_string(Nz / minSlices) + " ranks)");
  a = lo[r];
  b = lo[r + 1];
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
  // distributed storage: every rank keeps the slices and half-cells [a - G, b + G) of its block and ghosts
  // (the beam push reaches down to half-cell a - pf - 2, the filter and the exchanges up to b + pf - 1); all
  // current arrays of a rank share one layout (filterBeam indexes the halo / front arrays like c1b)
  const int G = pf + 4, M = S.M, N = S.g.N;
  sa = a - G;
  sb = b + G;
  if (std::getenv("GEMPIC_FULL_STORAGE")) { sa = 0; sb = Nz; }  // the previous behaviour (testing)
  auto shrinkF = [&](Fields& F) {
    Fields t;
    t.initRange(M, Nz, N, sa, sb + 1);
    t.copyRange(F, sa, sb + 1);
    F = std::move(t);
  };
  auto shrinkC = [&](Currents& C) {
    Currents t;
    t.initRange(M, Nz, N, sa, sb);
    t.copyRange(C, sa, sb);
    C = std::move(t);
  };
  shrinkF(S.fold);
  shrinkF(S.fnew);
  for (auto* C : {&S.c0b, &S.c0s, &S.c1b, &S.c1s}) shrinkC(*C);
  for (auto* C : {&H, &Tr, &Rw, &Snap}) C->initRange(M, Nz, N, sa, sb);
  predBelowR.assign(S.predR.size(), cd(0.0));
  predBelowP.assign(S.predP.size(), cd(0.0));
  predRfront.assign(S.predR.size(), cd(0.0));
  predPfront.assign(S.predP.size(), cd(0.0));
  // number of the last step (same arithmetic as the serial loop)
  lastStep = S.step_no + S.cfg.grid.stepsFrom(S.s0);
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
// one step s^0 -> s^1 on the slices [a, b) (port of Simulation::step)
// ----------------------------------------------------------------------------
void PipelineRun::step() {
  const int N = S.g.N, M = S.M, Nz = S.Nz;
  S.beginStep();
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

  // ---- from the front rank: slice b, its prediction, currents, particles
  S.stream.clear();
  inContinue.clear(); inFinished.clear(); inTransfer.clear();
  if (r < P - 1) {
    Buffer B = comm::recv(r + 1, TAG_F);
    if (B.get<int>() != k) throw std::runtime_error("MPI pipeline: step mismatch (F)");
    lastStep = std::min(lastStep, B.get<int>());
    wallStop = wallStop || B.get<int>() != 0;
    getSlice(B, S.fnew, b);
    predRfront = B.getVec<cd>();
    predPfront = B.getVec<cd>();
    getHalf(B, S.c1b, b - 1, true);  // raw beam deposits of the front rank (filtered below)
    getHalf(B, S.c1b, b - 2, true);
    getHalf(B, Tr, b - 1, false);     // its trapped deposits (added after the filter)
    getHalf(B, Tr, b - 2, false);
    for (int k = b; k < b + pf; ++k) getHalf(B, Rw, k, false);  // its raw beam current above b
    S.stream = B.getVec<Particle>();
    inContinue = B.getVec<Particle>();
    inFinished = B.getVec<Particle>();
    inTransfer = B.getVec<Particle>();
  }

  // ---- beam push (fields of s^0; slices a-2 ... b are available)
  const int nsub = std::max(1, int(std::ceil(Delta / num.beam_dz_max)));
  const double dz = Delta / nsub;
  const double rmax = S.g.rho[N];
  // lowest zeta whose current this rank needs: centre of half-cell a-pf-2
  const double lowZ = (r > 0 && a - pf - 2 >= -1) ? zv.centre(a - pf - 2) : -1e300;
  // returns true if the push must continue on the rank behind (zeta < lowZ)
  auto push = [&](Particle& p, bool count) {
    for (int sub = p.sub; sub < nsub && p.status == 0; ++sub) {
      if (p.zeta < lowZ) { p.sub = sub; return true; }
      double gm = gammaOf(p);
      const double vz = p.uz / gm;
      if (vz < num.vz_min) { p.status = 3; if (count) ++myRestreamed; break; }
      const double dt = dz / vz;
      double w;
      const int na = zv.cell(p.zeta, w);
      if (na < 0 || na >= Nz) { p.status = 2; break; }
      EB f = S.gather(S.fold, na, 1.0 - w, na + 1, w, p.x, p.y);
      addExternal(S.ext, p.x, p.y, p.z, f, dz);
      if (num.beam_vay) vay(p.ux, p.uy, p.uz, f, p.qm, dt);
      else boris(p.ux, p.uy, p.uz, f, p.qm, dt);
      gm = gammaOf(p);
      p.x += p.ux / gm * dt; p.y += p.uy / gm * dt;
      p.zeta += (p.uz / gm - 1.0) * dt;
      p.z += dz;
      if (p.zeta < 0.0 || std::hypot(p.x, p.y) >= rmax) p.status = 2;
    }
    p.sub = 0;
    return false;
  };
  std::vector<Particle> newBeam, handCont, handFin, transferDown;
  auto afterPush = [&](Particle& p) {
    if (p.status == 3) { Particle q = p; q.status = 0; S.pending.push_back(q); return; }
    if (p.status != 0) return;
    if (beamOwner(p) == r) { newBeam.push_back(p); return; }
    if (inHalo(kh(p))) S.depositBeamParticle(H, p);  // completes the current below the block
    handFin.push_back(p);
  };
  if (threads > 1) {
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
      if (cont[i]) handCont.push_back(S.beam[i]);
      else afterPush(S.beam[i]);
    }
    for (auto& p : inContinue) {
      if (push(p, true)) handCont.push_back(p);
      else afterPush(p);
    }
  } else {
    for (auto* list : {&S.beam, &inContinue})
      for (auto& p : *list) {
        if (push(p, true)) handCont.push_back(p);
        else afterPush(p);
      }
  }
  for (const auto& p : newBeam) S.depositBeamParticle(S.c1b, p);
  for (const auto& p : inFinished) {  // pushed by a rank in front
    if (beamOwner(p) == r) { S.depositBeamParticle(S.c1b, p); newBeam.push_back(p); }
    else { if (inHalo(kh(p))) S.depositBeamParticle(H, p); handFin.push_back(p); }
  }
  for (const auto& p : inTransfer) {  // trapped elsewhere, already deposited
    if (beamOwner(p) == r) newBeam.push_back(p);
    else transferDown.push_back(p);
  }
  S.beam.swap(newBeam);
  for (const auto& p0 : halo) {  // copies of rank r-1's particles near the boundary
    Particle q = p0;
    if (!push(q, false) && q.status == 0 && inHalo(kh(q))) S.depositBeamParticle(H, q);
  }
  filterBeam();

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

  // ---- zeta march over the own slices
  std::vector<Particle> newForFront;
  int done = 0;
  for (int n = b - 1; n >= a; --n) {
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
        S.solver.solveSlice(m, n, &fo.Rr[fo.f(m, n, 0)], &fo.Rp[fo.c(m, n, 0)], &fo.Ez[fo.c(m, n, 0)],
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
      for (int m = 0; m <= M; ++m) {
        for (int i = 1; i <= N; ++i) {
          cd jt = S.c1b.jr[S.c1b.f(m, n, i)] + S.c1s.jr[S.c1s.f(m, n, i)];
          S.fnew.Rr[S.fnew.f(m, n + 1, i)] -= q * (jt - S.predRprev[m * (N + 1) + i]);
        }
        for (int i = 0; i < N; ++i) {
          cd jt = S.c1b.jp[S.c1b.c(m, n, i)] + S.c1s.jp[S.c1s.c(m, n, i)];
          S.fnew.Rp[S.fnew.c(m, n + 1, i)] -= q * (jt - S.predPprev[m * N + i]);
        }
      }
      absorbSlice(S.fnew, n + 1, s1);  // slice n+1 is final
    }
    S.predRprev.swap(S.predR);
    S.predPprev.swap(S.predP);
    if (++done == pf + 3 && r < P - 1) sendFront(k, newForFront);
  }
  if (r == 0) absorbSlice(S.fnew, 0, s1);
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
  S.s0 = s1;
  ++S.step_no;
}

// Filter of the beam current (Simulation::filterBeamCurrent) on the own
// half-cells and on half-cell a-1 (prediction at slice a). The raw current is
// assembled over [a-pf-1, b+pf-1]: own deposits, the halo deposits below a,
// the front rank's raw current above b. Values within pf of the ends of this
// range are not needed. The raw own deposits below a and the raw current of
// [a, a+pf-1] are kept for the rear rank (Snap); c1b below a is reset and then
// collects only the deposits of particles trapped during the march.
void PipelineRun::filterBeam() {
  const int N = S.g.N, M = S.M, Nz = S.Nz;
  const int k0 = (r > 0) ? a - 2 : a, k1 = std::min(b - 1, a + pf - 1);
  for (int k = k0; k <= k1; ++k) {
    Snap.zeroHalf(k);
    for (int m = 0; m <= M; ++m) {
      for (int i = 0; i <= N; ++i) Snap.jr[Snap.f(m, k, i)] = S.c1b.jr[S.c1b.f(m, k, i)];
      for (int i = 0; i < N; ++i) {
        Snap.jp[Snap.c(m, k, i)] = S.c1b.jp[S.c1b.c(m, k, i)];
        Snap.jz[Snap.c(m, k, i)] = S.c1b.jz[S.c1b.c(m, k, i)];
      }
    }
    for (int i = 0; i < N; ++i) {
      Snap.rho[Snap.d(k, i)] = S.c1b.rho[S.c1b.d(k, i)];
      Snap.rhoi[Snap.d(k, i)] = S.c1b.rhoi[S.c1b.d(k, i)];
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
        if (below && r > 0) (*below)[size_t(m) * L + i] = col[a - 1 - klo];
      }
  };
  auto fi = [&](int m, int k, int i) { return S.c1b.f(m, k, i); };
  auto ci = [&](int m, int k, int i) { return S.c1b.c(m, k, i); };
  run(S.c1b.jr, H.jr, Rw.jr, &predBelowR, N + 1, fi);
  run(S.c1b.jp, H.jp, Rw.jp, &predBelowP, N, ci);
  run(S.c1b.jz, H.jz, Rw.jz, nullptr, N, ci);
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

// after the first pf+3 slices: what the front rank needs for its next step
void PipelineRun::sendFront(int k, const std::vector<Particle>& newForFront) {
  Buffer B;
  B.put<int>(k);
  for (int n = b - pf - 2; n <= b - 1; ++n) putSlice(B, S.fnew, n);
  putSliceR(B, S.fnew, b);
  putHalf(B, S.c1b, b - 1);
  putHalf(B, S.c1s, b - 1);
  putHalf(B, S.c1b, b);  // deposits of the particles trapped at slice b-1
  std::vector<Particle> window;
  for (const auto& p : S.beam)
    if (kh(p) >= b - pf - 1 && kh(p) <= b - 1) window.push_back(p);
  B.putVec(window);
  B.putVec(newForFront);
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
  getHalf(B, S.c0b, a, true);
  halo = B.getVec<Particle>();
  const std::vector<Particle> nb = B.getVec<Particle>();
  S.beam.insert(S.beam.end(), nb.begin(), nb.end());
}

// ----------------------------------------------------------------------------
// diagnostics: partial sums (and data of output / checkpoint steps) to rank 0
// ----------------------------------------------------------------------------
void PipelineRun::sendDiag(int k) {
  Buffer B;
  B.put<int>(k);
  putDiag(B, S.diagSums(a, b));
  B.put(myRestreamed);
  B.put(myTrapped);
  if (filesAt(k) || checkpointAt(k) || lasyAt(k)) {
    putRange(B, S.fold, a, b);
    putDens(B, S.c0b, a, b);
    putDens(B, S.c0s, a, b);
    B.putVec(S.beam);
  }
  if (checkpointAt(k)) {
    putJ(B, S.c0b, a, b);
    putJ(B, S.c0s, a, b);
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
}

void PipelineRun::reportStep(int k) {
  DiagSums D = S.diagSums(a, b);
  long long re = myRestreamed, tr = myTrapped;
  const bool files = filesAt(k), ck = checkpointAt(k), fields = files || ck || lasyAt(k);
  std::vector<Particle> gbeam, gpend;
  if (files || ck) gbeam = S.beam;
  if (ck) gpend = S.pendingNext;
  // output steps: the whole window, assembled on rank 0 only for the output (the ranks store their blocks)
  Fields gF;
  Currents gCb, gCs;
  if (fields) {
    const int M = S.M, Nz = S.Nz, N = S.g.N;
    gF.init(M, Nz, N);
    gF.copyRange(S.fold, a, (P == 1) ? Nz + 1 : b);
    gCb.init(M, Nz, N);
    gCs.init(M, Nz, N);
    gCb.copyRange(S.c0b, a, b);
    gCs.copyRange(S.c0s, a, b);
  }
  std::string rngState;
  for (int q = 1; q < P; ++q) {
    Buffer B = comm::recv(q, TAG_DIAG);
    if (B.get<int>() != k) throw std::runtime_error("MPI pipeline: step mismatch (diagnostics)");
    D.add(getDiag(B));
    re += B.get<long long>();
    tr += B.get<long long>();
    if (fields) {
      getRange(B, gF, lo[q], lo[q + 1]);
      getDens(B, gCb, lo[q], lo[q + 1]);
      getDens(B, gCs, lo[q], lo[q + 1]);
      const std::vector<Particle> pb = B.getVec<Particle>();
      gbeam.insert(gbeam.end(), pb.begin(), pb.end());
    }
    if (ck) {
      getJ(B, gCb, lo[q], lo[q + 1]);
      getJ(B, gCs, lo[q], lo[q + 1]);
      const std::vector<Particle> pp = B.getVec<Particle>();
      gpend.insert(gpend.end(), pp.begin(), pp.end());
      const std::vector<char> rs = B.getVec<char>();
      if (q == P - 1) rngState.assign(rs.begin(), rs.end());
    }
  }
  S.n_restreamed = baseRestreamed + long(re);
  S.n_trapped_total = baseTrapped + long(tr);
  auto swapWindow = [&] {
    if (!fields) return;
    std::swap(S.fold, gF);
    std::swap(S.c0b, gCb);
    std::swap(S.c0s, gCs);
  };
  swapWindow();
  if (files) S.beam.swap(gbeam);
  S.report(D, files);
  if (files) S.beam.swap(gbeam);
  if (ck) {
    S.beam.swap(gbeam);
    S.pendingNext.swap(gpend);
    if (P > 1) {
      std::istringstream is(rngState);
      is >> S.rng;  // rank 0 does not load plasma itself: its generator is unused
    }
    std::printf("checkpoint: %s\n", S.writeCheckpoint().c_str());
    S.beam.swap(gbeam);
    S.pendingNext.swap(gpend);
  }
  swapWindow();
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
  t0 = std::chrono::steady_clock::now();
  try {
    if (r == 0 && !rc.restarted) S.diagnostics(true);  // initial state (global on every rank)
    partition();
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
        sendDiag(k - 1);
      }
      step();
      if (!finiteOwn()) {
        std::fprintf(stderr, "non-finite field at step %d on rank %d\n", S.step_no, r);
        comm::abort(2);
      }
      if (r == 0) reportStep(k);
    }
    if (r > 0 && lastStep > k0) {
      recvRear();
      sendDiag(lastStep);
    }
    comm::waitAll();
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
