// ============================================================================
//  pipeline.cpp -- MPI pipeline along zeta: partition, boundaries, messages, the run (class in pipeline_internal.hpp,
//  the step in pipeline_step.cpp)
// ============================================================================
#include "pipeline_internal.hpp"

namespace gem {

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
  H.init(S.M, Nz, S.g.N);
  Tr.init(S.M, Nz, S.g.N);
  Rw.init(S.M, Nz, S.g.N);
  Snap.init(S.M, Nz, S.g.N);
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
      const cd inc = (ff && m == 1) ? cd(S.ffInc(true, n, i, s)) : cd(0.0);
      cd& R = F.Rr[F.f(m, n, i)];
      R = inc + d * (R - inc);
    }
    for (int i = 0; i < S.g.N; ++i) {
      const double d = damp(S.g.c[i]);
      if (d == 1.0) continue;
      const cd inc = (ff && m == 1) ? cd(0.0, 1.0) * S.ffInc(false, n, i, s) : cd(0.0);
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
      Snap.rho[size_t(k) * N + i] = S.c1b.rho[size_t(k) * N + i];
      Snap.rhoi[size_t(k) * N + i] = S.c1b.rhoi[size_t(k) * N + i];
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
        S.c1b.rho[size_t(k) * N + i] += Tr.rho[size_t(k) * N + i];
        S.c1b.rhoi[size_t(k) * N + i] += Tr.rhoi[size_t(k) * N + i];
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
  if (filesAt(k) || checkpointAt(k)) {
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
  const bool files = filesAt(k), ck = checkpointAt(k);
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
      getRange(B, S.fold, lo[q], lo[q + 1]);
      getDens(B, S.c0b, lo[q], lo[q + 1]);
      getDens(B, S.c0s, lo[q], lo[q + 1]);
      const std::vector<Particle> pb = B.getVec<Particle>();
      gbeam.insert(gbeam.end(), pb.begin(), pb.end());
    }
    if (ck) {
      getJ(B, S.c0b, lo[q], lo[q + 1]);
      getJ(B, S.c0s, lo[q], lo[q + 1]);
      const std::vector<Particle> pp = B.getVec<Particle>();
      gpend.insert(gpend.end(), pp.begin(), pp.end());
      const std::vector<char> rs = B.getVec<char>();
      if (q == P - 1) rngState.assign(rs.begin(), rs.end());
    }
  }
  S.n_restreamed = baseRestreamed + long(re);
  S.n_trapped_total = baseTrapped + long(tr);
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
