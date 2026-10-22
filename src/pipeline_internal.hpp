// ============================================================================
//  pipeline_internal.hpp -- MPI parallelization by domain decomposition along zeta
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
//  checkpoints (same files as a serial run).
//  With one rank the pipelined step is identical to Simulation::step().
// ============================================================================
#pragma once
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

namespace pipedetail {
enum { TAG_F = 101, TAG_R = 102, TAG_DIAG = 103 };
constexpr int MIN_SLICES = 4;  // per rank without beam filter (the R message needs 3 processed slices)
}  // namespace pipedetail
using namespace pipedetail;

class PipelineRun {
 public:
  PipelineRun(Simulation& s, const RunControl& c) : S(s), rc(c) {}
  int run();

 private:
  Simulation& S;
  RunControl rc;
  int P = 1, r = 0, a = 0, b = 0;
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
    B.putRaw(&C.rho[size_t(k) * N], sizeof(double) * N);
    B.putRaw(&C.rhoi[size_t(k) * N], sizeof(double) * N);
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
    taked(&C.rho[size_t(k) * N]);
    taked(&C.rhoi[size_t(k) * N]);
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
    B.putRaw(&C.rho[size_t(k0) * S.g.N], sizeof(double) * size_t(k1 - k0) * S.g.N);
    B.putRaw(&C.rhoi[size_t(k0) * S.g.N], sizeof(double) * size_t(k1 - k0) * S.g.N);
  }
  void getDens(Buffer& B, Currents& C, int k0, int k1) {
    B.getRaw(&C.rho[size_t(k0) * S.g.N], sizeof(double) * size_t(k1 - k0) * S.g.N);
    B.getRaw(&C.rhoi[size_t(k0) * S.g.N], sizeof(double) * size_t(k1 - k0) * S.g.N);
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

}  // namespace gem
