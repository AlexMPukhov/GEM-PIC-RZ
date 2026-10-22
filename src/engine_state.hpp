// ============================================================================
//  engine_state.hpp -- device state of the parallel engine and its phase timing (src/engine*.cpp)
// ============================================================================
#pragma once
#include "engine_kernels.hpp"

namespace gem {

// ============================================================================
//  engine state and host-side control
// ============================================================================
// optional timing of the phases (environment variable GEMPIC_PROFILE=1);
// every timed phase ends with a device synchronization, which itself costs time
struct Profile {
  enum { BEAM, LOAD, PREP, COPY, PUSH, SOLVE, POST, CORR, END, DOWNLOAD, COUNT };
  bool on = false;
  double t[COUNT] = {};
  std::chrono::steady_clock::time_point last;
  void start() { if (on) { par::fence(); last = std::chrono::steady_clock::now(); } }
  void stop(int k) {
    if (!on) return;
    par::fence();
    const auto now = std::chrono::steady_clock::now();
    t[k] += std::chrono::duration<double>(now - last).count();
    last = now;
  }
  void report() const {
    if (!on) return;
    static const char* names[COUNT] = {"beam push+compaction", "fresh plasma (host)", "slice preparation",
                                       "backup/restore", "streaming push+deposit", "field solve",
                                       "trapped/freed", "R correction", "end of step", "download"};
    double tot = 0;
    for (double x : t) tot += x;
    std::printf("engine profile (s):");
    for (int k = 0; k < COUNT; ++k) std::printf("  %s %.3f", names[k], t[k]);
    std::printf("  total %.3f\n", tot);
  }
};

struct EngineState {
  Profile prof;
  // grid, ions, matrices
  par::Array<double> rho, c, Af, Ac, Vf, Vc, kC, kF;
  par::Array<dev::IonV> ions;
  par::Array<dev::AdkV> adk;
  par::Array<dev::LevelV> levels;
  par::Array<double> triData;
  par::Array<dev::TriV> Tf, Tc;
  int triVersion = -1;  // FieldSolver::version() of the uploaded matrices
  par::Array<double> zz, zhs, zzc, zwlo, zwhi;  // longitudinal grid
  par::Array<int> cls;                          // matrix class of the slices
  par::Array<UndulatorV> und;                   // external fields
  par::Array<double> undK;                      // taper tables
  par::Array<MagnetV> mag;                      // transport elements
  // factorized matrices of the solver -> triData (same layout for every step)
  void uploadMatrices(const FieldSolver& fs, int M) {
    std::vector<double> td;
    const int nT = fs.classes() * (M + 1);
    for (int t = 0; t < nT; ++t)
      for (const Tridiag* T : {&fs.Tf[t], &fs.Tc[t]}) {
        td.insert(td.end(), T->a.begin(), T->a.end());
        td.insert(td.end(), T->cp.begin(), T->cp.end());
        td.insert(td.end(), T->den.begin(), T->den.end());
      }
    if (triVersion < 0) {  // first call: allocate and build the views
      triData.assign(td.data(), td.size());
      std::vector<dev::TriV> tf, tc;
      size_t off = 0;
      for (int t = 0; t < nT; ++t)
        for (int which = 0; which < 2; ++which) {
          const Tridiag& T = which ? fs.Tc[t] : fs.Tf[t];
          const double* base = triData.data() + off;
          (which ? tc : tf).push_back(dev::TriV{T.i0, T.n, base, base + T.n, base + 2 * size_t(T.n)});
          off += 3 * size_t(T.n);
        }
      Tf.assign(tf.data(), tf.size());
      Tc.assign(tc.data(), tc.size());
    } else {
      par::upload(triData.data(), td.data(), td.size());
    }
    triVersion = fs.version();
  }
  dev::Params P{};
  // fields and currents: index 0/1 swap every step
  dev::FieldStore F[2];
  dev::CurrStore Cb[2], Cs[2];
  int iold = 0;  // F[iold] = level s^0 (fold), Cb[iold] = c0b, Cs[iold] = c0s
  par::Array<cplx> pred[2], predP[2];  // predicted j^1 (index ip = current slice)
  int ip = 0;
  par::Array<cplx> solveScratch;
  // plasma.smooth_length
  par::Array<cplx> filtBand, filtScratch;
  par::Array<double> filtS;
  dev::FiltV filt{0, 0, nullptr, nullptr, nullptr, nullptr};
  std::vector<dev::CurrStore> local;  // CPU: per-thread currents of one half-cell
  // particles
  par::Array<Particle> stream, backup, born, beam, beamTmp, pnext, pprev;
  par::Array<int> cnt;
  par::Array<unsigned long long> lcnt;
  par::Array<double> bound, freed, exitSt;
  par::Array<cplx> ffGf, ffGc;  // flying focus: incident field G(r, zeta) (Simulation::ffGf, ffGc)
  size_t nf = 0, nc = 0, nd = 0;
  int zStride = 1, maxSpan = 0;
  // host copies of the counters
  int nBeam = 0, nPprev = 0;
  std::vector<int> hcnt = std::vector<int>(dev::C_COUNT, 0);

  void pushCounters() { par::upload(cnt.data(), hcnt.data(), hcnt.size()); }
  void pullCounters() {
    par::fence();
    par::download(hcnt.data(), cnt.data(), hcnt.size());
    if (hcnt[dev::C_OVERFLOW]) throw std::runtime_error("parallel engine: particle buffer overflow (internal error)");
  }
};

}  // namespace gem
