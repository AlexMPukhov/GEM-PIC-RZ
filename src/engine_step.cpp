// ============================================================================
//  engine_step.cpp -- the time step of the parallel engine
// ============================================================================
#include "engine_state.hpp"

namespace gem {

// ----------------------------------------------------------------------------
//  one step s^0 -> s^1 (Simulation::step)
// ----------------------------------------------------------------------------
void Engine::step() {
  using namespace dev;
  EngineState& E = *st;
  // step of this s-level and its matrices (Simulation::beginStep)
  E.P.Delta = sim.Delta;
  E.P.trapExtra = (sim.cfg.numerics.trap_advance - 1.0) * sim.Delta;
  E.P.theta = sim.cfg.numerics.theta;
  if (sim.solver.version() != E.triVersion) E.uploadMatrices(sim.solver, E.P.M);
  const Params P = E.P;
  const int M = P.M, Nz = P.Nz, N = P.N;
  const double Delta = P.Delta;
  const double s1 = sim.s0 + Delta;
  const int iters = sim.cfg.numerics.iterations;
  const int iold = E.iold, inew = 1 - iold;
  const FieldsV fo = E.F[iold].view(M, Nz, N), fn = E.F[inew].view(M, Nz, N);
  const CurrV c0b = E.Cb[iold].view(M, Nz, N), c0s = E.Cs[iold].view(M, Nz, N);
  const CurrV c1b = E.Cb[inew].view(M, Nz, N), c1s = E.Cs[inew].view(M, Nz, N);
  E.Cb[inew].zero();
  E.Cs[inew].zero();
  par::zero(E.freed.data(), E.freed.capacity());
  int* cnt = E.cnt.data();
  unsigned long long* lcnt = E.lcnt.data();
  Profile& prof = E.prof;
  prof.start();

  // ---- 1. beam push (fields of s^0)
  {
    const int nsub = std::max(1, int(std::ceil(Delta / sim.cfg.numerics.beam_dz_max)));
    const double dz = Delta / nsub;
    Particle* beam = E.beam.data();
    par::forDyn(&cnt[C_BEAM], nullptr, E.nBeam,
                GEM_LAMBDA(long i) { pushBeamParticle(P, fo, beam[i], nsub, dz, lcnt); });
  }
  // compaction: alive -> beamTmp, slower than vz_min -> waiting list (pprev)
  E.beamTmp.reserve(std::max(1, E.nBeam));
  E.pprev.reserve(size_t(E.nPprev) + E.nBeam + 1, true);
  {
    E.hcnt[C_TMP] = 0;
    E.hcnt[C_PPREV] = E.nPprev;
    E.pushCounters();
    Particle *beam = E.beam.data(), *tmp = E.beamTmp.data(), *pprev = E.pprev.data();
    par::forDyn(&cnt[C_BEAM], nullptr, E.nBeam, GEM_LAMBDA(long i) {
      const Particle& p = beam[i];
      if (p.status == 0) tmp[par::atomicFetchAdd(&cnt[C_TMP], 1)] = p;
      else if (p.status == 3) {
        Particle q = p;
        q.status = WAITING;
        pprev[par::atomicFetchAdd(&cnt[C_PPREV], 1)] = q;
      }
    }, /*ordered=*/true);
    E.pullCounters();
    E.beam.swap(E.beamTmp);
    E.nBeam = E.hcnt[C_TMP];
    E.nPprev = E.hcnt[C_PPREV];
    E.hcnt[C_BEAM] = E.nBeam;
  }

  prof.stop(Profile::BEAM);
  // ---- 2. fresh plasma (host, same random numbers as the serial code)
  sim.loadFresh();
  const size_t nFresh = sim.stream.size();
  long bornCap = 0;  // (ids of the ionization RNG are set by loadFresh)
  for (size_t k = 0; k < nFresh; ++k) {
    const Particle& p = sim.stream[k];
    if (isIon(p.sp) && p.sp != SPECIES_TRAPPED_ION) {
      const auto& I = sim.cfg.ions[ionIndex(p.sp)];
      if (I.ionization) bornCap += std::max(0, I.Zmax - p.Zi);
    }
  }
  bornCap += long(E.nPprev) * E.maxSpan;
  const long streamCap = long(nFresh) + E.nPprev + bornCap + 1;
  E.stream.reserve(streamCap);
  if (iters > 1) E.backup.reserve(streamCap);
  E.born.reserve(bornCap + 1);
  E.pnext.reserve(streamCap);
  E.beam.reserve(size_t(E.nBeam) + streamCap, true);
  E.beamTmp.reserve(E.beam.capacity());
  par::upload(E.stream.data(), sim.stream.data(), nFresh);
  {
    Particle *stream = E.stream.data(), *pprev = E.pprev.data();
    const long off = long(nFresh);
    par::forN(E.nPprev, GEM_LAMBDA(long i) {
      Particle q = pprev[i];
      q.status = WAITING;
      stream[off + i] = q;
    });
  }
  E.hcnt[C_STREAM] = int(nFresh) + E.nPprev;
  E.hcnt[C_BORN] = 0;
  E.hcnt[C_PNEXT] = 0;
  E.pushCounters();
  E.nPprev = 0;

  // ---- 3. deposit the beam list at s^1
  {
    Particle* beam = E.beam.data();
    par::forDyn(&cnt[C_BEAM], nullptr, E.nBeam, GEM_LAMBDA(long i) { depositBeam(P, c1b, beam[i]); });
  }
  // (1,2,1) filter of the beam current along zeta (Simulation::filterBeamCurrent), one thread per column
  if (const int passes = sim.cfg.numerics.beam_filter) {
    const struct { cplx* v; int L; } cols[3] = {{c1b.jr, N + 1}, {c1b.jp, N}, {c1b.jz, N}};
    for (const auto& cl : cols) {
      cplx* v = cl.v;
      const int L = cl.L;
      par::forN(long(M + 1) * L, GEM_LAMBDA(long col) {
        const int m = int(col / L), i = int(col % L);
        for (int p = 0; p < passes; ++p) {
          cplx prev;
          for (int k = 0; k < Nz; ++k) {
            cplx& x = v[(size_t(m) * Nz + k) * L + i];
            const cplx cur = x;
            const cplx next = (k + 1 < Nz) ? v[(size_t(m) * Nz + k + 1) * L + i] : cplx();
            x = 0.25 * prev + 0.5 * cur + 0.25 * next;
            prev = cur;
          }
        }
      });
    }
  }
  const size_t npf = size_t(M + 1) * (N + 1), npc = size_t(M + 1) * N;
  par::zero(E.pred[1 - E.ip].data(), npf);
  par::zero(E.predP[1 - E.ip].data(), npc);
  prof.stop(Profile::LOAD);

  // ---- 4. zeta march
  Particle *stream = E.stream.data(), *backup = E.backup.data(), *born = E.born.data();
  Particle *beamL = E.beam.data(), *pnext = E.pnext.data();
  const int beamCap = int(E.beam.capacity()), pnextCap = int(E.pnext.capacity());
  const int bornCapI = int(bornCap + 1), streamCapI = int(streamCap);
  double* boundPtr = E.bound.data();
  double* freed = E.freed.data();
  const TriV *Tf = E.Tf.data(), *Tc = E.Tc.data();
  cplx* scratch = E.solveScratch.data();
  const bool ions = P.nIons > 0;
  const double vz_beam = P.vz_beam;
  const uint64_t keyStep = ionKeyStep(P.seed, sim.step_no);

  for (int n = Nz - 1; n >= 0; --n) {
    cplx* predR = E.pred[E.ip].data();
    cplx* predPc = E.predP[E.ip].data();
    const cplx* predRprev = E.pred[1 - E.ip].data();
    const cplx* predPprev = E.predP[1 - E.ip].data();
    for (int it = 0; it < iters; ++it) {
      // prepare: predicted j^1 at n-1, predictor fields, zero half-cell n, |E| bound
      par::teams(1, GEM_LAMBDA(const Team& t) {
        if (it == 0) {
          par::teamFor(t, int(npf), [&](int idx) {
            const int m = idx / (N + 1), i = idx % (N + 1);
            predR[idx] = (n > 0) ? c1b.jr[c1b.f(m, n - 1, i)] + c0s.jr[c0s.f(m, n - 1, i)] : cplx();
            fn.Rr[fn.f(m, n, i)] = fo.Rr[fo.f(m, n, i)];
            fn.Lr[fn.f(m, n, i)] = fo.Lr[fo.f(m, n, i)];
            fn.Bz[fn.f(m, n, i)] = fo.Bz[fo.f(m, n, i)];
          });
          par::teamFor(t, int(npc), [&](int idx) {
            const int m = idx / N, i = idx % N;
            predPc[idx] = (n > 0) ? c1b.jp[c1b.c(m, n - 1, i)] + c0s.jp[c0s.c(m, n - 1, i)] : cplx();
            fn.Rp[fn.c(m, n, i)] = fo.Rp[fo.c(m, n, i)];
            fn.Lp[fn.c(m, n, i)] = fo.Lp[fo.c(m, n, i)];
            fn.Ez[fn.c(m, n, i)] = fo.Ez[fo.c(m, n, i)];
          });
        }
        par::teamFor(t, int(npf), [&](int idx) {
          const int m = idx / (N + 1), i = idx % (N + 1);
          c1s.jr[c1s.f(m, n, i)] = cplx();
        });
        par::teamFor(t, int(npc), [&](int idx) {
          const int m = idx / N, i = idx % N;
          c1s.jp[c1s.c(m, n, i)] = cplx();
          c1s.jz[c1s.c(m, n, i)] = cplx();
        });
        par::teamFor(t, N, [&](int i) {
          c1s.rho[size_t(n) * N + i] = 0.0;
          c1s.rhoi[size_t(n) * N + i] = 0.0;
        });
        if (t.leader()) cnt[C_BORN] = 0;
        if (ions) {
          t.barrier();
          double bnd = 0.0;
          for (int sl = n + 1; sl >= n; --sl) {
            double b = 0.0;
            for (int m = 0; m <= M; ++m) {
              double er = 0, ep = 0, ez = 0;
              par::teamFor(t, N + 1, [&](int i) {
                er = fmax(er, 0.5 * cabs(fn.Rr[fn.f(m, sl, i)] + fn.Lr[fn.f(m, sl, i)]));
              });
              par::teamFor(t, N, [&](int i) {
                ep = fmax(ep, 0.5 * cabs(fn.Rp[fn.c(m, sl, i)] + fn.Lp[fn.c(m, sl, i)]));
                ez = fmax(ez, cabs(fn.Ez[fn.c(m, sl, i)]));
              });
              er = par::teamMax(t, er);
              ep = par::teamMax(t, ep);
              ez = par::teamMax(t, ez);
              b += (m == 0 ? 1.0 : 2.0) * sqrt(er * er + ep * ep + ez * ez);
            }
            bnd = fmax(bnd, 2.0 * b);
          }
          if (t.leader()) *boundPtr = bnd;
        }
      });
      prof.stop(Profile::PREP);
      // predictor-corrector: every pass starts from the same particles
      if (iters > 1) {
        if (it == 0) par::forDyn(&cnt[C_STREAM], nullptr, streamCapI, GEM_LAMBDA(long i) { backup[i] = stream[i]; });
        else par::forDyn(&cnt[C_STREAM], nullptr, streamCapI, GEM_LAMBDA(long i) { stream[i] = backup[i]; });
      }
      prof.stop(Profile::COPY);
      // push, deposit, ionize
      const Ionizer Z{born, cnt, bornCapI, ionKeySlice(keyStep, n, it)};
#if defined(GEM_BACKEND_CUDA) || defined(GEM_EMULATE_GPU)
      par::forDyn(&cnt[C_STREAM], nullptr, streamCapI, GEM_LAMBDA(long i) {
        streamParticle<true>(P, fn, c1s, n, stream[i], n, s1, ions ? *boundPtr : 0.0, Z);
      });
#else
      // CPU: every thread deposits into a private buffer of half-cell n (no
      // atomics); the buffers are summed in a fixed order afterwards
      {
        const long ns = cnt[C_STREAM];
        const double bnd = ions ? *boundPtr : 0.0;
        const int nt = int(E.local.size());
#if defined(_OPENMP)
#pragma omp parallel num_threads(nt) if (ns > 256)
#endif
        {
#if defined(_OPENMP)
          const int tid = omp_get_thread_num();
#else
          const int tid = 0;
#endif
          E.local[tid].zero();
          const CurrV L = E.local[tid].view(M, 1, N);
#if defined(_OPENMP)
#pragma omp for schedule(static)
#endif
          for (long i = 0; i < ns; ++i) streamParticle<false>(P, fn, L, 0, stream[i], n, s1, bnd, Z);
        }
        const int nused = (ns > 256) ? nt : 1;
        std::vector<CurrV> lv;
        for (int k = 0; k < nused; ++k) lv.push_back(E.local[k].view(M, 1, N));
        const CurrV* L = lv.data();
        par::forN(long(npf), GEM_LAMBDA(long idx) {
          const int m = int(idx / (N + 1)), i = int(idx % (N + 1));
          cplx a;
          for (int k = 0; k < nused; ++k) a += L[k].jr[idx];
          c1s.jr[c1s.f(m, n, i)] += a;
        });
        par::forN(long(npc), GEM_LAMBDA(long idx) {
          const int m = int(idx / N), i = int(idx % N);
          cplx a, b;
          for (int k = 0; k < nused; ++k) { a += L[k].jp[idx]; b += L[k].jz[idx]; }
          c1s.jp[c1s.c(m, n, i)] += a;
          c1s.jz[c1s.c(m, n, i)] += b;
        });
        par::forN(long(N), GEM_LAMBDA(long i) {
          double a = 0, b = 0;
          for (int k = 0; k < nused; ++k) { a += L[k].rho[i]; b += L[k].rhoi[i]; }
          c1s.rho[size_t(n) * N + i] += a;
          c1s.rhoi[size_t(n) * N + i] += b;
        });
      }
#endif
      if (E.filt.on) {  // plasma.smooth_length: one team per mode
        const dev::FiltV FV = E.filt;
        cplx* fs = E.filtScratch.data();
        par::teams(M + 1, GEM_LAMBDA(const Team& t) {
          filterTeam(t, FV, t.league, n, c1s, fs + size_t(t.league) * (4 * size_t(N) + 1));
        });
      }
      prof.stop(Profile::PUSH);
      // field solve, one team per mode
      par::teams(M + 1, GEM_LAMBDA(const Team& t) {
        const int m = t.league;
        solveSliceTeam(t, P, m, n, fo, fn, c0b, c0s, c1b, c1s, predR, predPc, Tf, Tc,
                       scratch + size_t(m) * SOLVE_SCRATCH * (N + 1));
      });
      prof.stop(Profile::SOLVE);
    }
    // freed electrons join the streaming list; trapped particles -> beam list
    // (fast) or next step (slow)
    par::forDyn(&cnt[C_STREAM], &cnt[C_BORN], streamCapI, GEM_LAMBDA(long i) {
      const int ns = cnt[C_STREAM];
      if (i >= ns) {
        const Particle& e = born[i - ns];
        if (i < streamCapI) stream[i] = e; else cnt[C_OVERFLOW] = 1;
        par::atomicAddD(&freed[ionizedIndex(e.sp)], e.w);
        return;
      }
      Particle& p = stream[i];
      if (p.status != 1) return;
      Particle t = p;
      t.status = 0;
      t.z = s1;
      p.status = 2;  // leaves the streaming list
      if (t.uz / gammaOf(t) >= vz_beam) {
        if (t.sp == SPECIES_PLASMA) t.sp = SPECIES_TRAPPED;
        else if (isIonizationElectron(t.sp)) t.sp = SPECIES_TRAPPED_IONIZED_BASE - ionizedIndex(t.sp);
        else if (isIon(t.sp)) t.sp = SPECIES_TRAPPED_ION;
        const int j = par::atomicFetchAdd(&cnt[C_BEAM], 1);
        if (j < beamCap) beamL[j] = t; else cnt[C_OVERFLOW] = 1;
        depositBeam(P, c1b, t);
        par::atomicAddU64(&lcnt[L_TRAPPED], 1ULL);
      } else {
        const int j = par::atomicFetchAdd(&cnt[C_PNEXT], 1);
        if (j < pnextCap) pnext[j] = t; else cnt[C_OVERFLOW] = 1;
      }
    }, /*ordered=*/true);
    prof.stop(Profile::POST);
    // half-cell n is final: correct R at slice n+1 for the predicted j^1_{n+1/2}
    par::teams(1, GEM_LAMBDA(const Team& t) {
      if (n + 1 <= Nz - 1) {
        const double q = 0.5 * P.theta * Delta * P.zv.wLo(n + 1);
        par::teamFor(t, int(npf), [&](int idx) {
          const int m = idx / (N + 1), i = idx % (N + 1);
          if (i == 0) return;
          const cplx jt = c1b.jr[c1b.f(m, n, i)] + c1s.jr[c1s.f(m, n, i)];
          fn.Rr[fn.f(m, n + 1, i)] -= q * (jt - predRprev[idx]);
        });
        par::teamFor(t, int(npc), [&](int idx) {
          const int m = idx / N, i = idx % N;
          const cplx jt = c1b.jp[c1b.c(m, n, i)] + c1s.jp[c1s.c(m, n, i)];
          fn.Rp[fn.c(m, n + 1, i)] -= q * (jt - predPprev[idx]);
        });
      }
      t.barrier();
      if (t.leader()) {
        cnt[C_STREAM] += cnt[C_BORN];
        cnt[C_BORN] = 0;
      }
    });
    prof.stop(Profile::CORR);
    E.ip = 1 - E.ip;
  }

  // ---- 5. charge states of the ions leaving the window, absorber
  double* exitSt = E.exitSt.data();
  const int zs = E.zStride;
  par::zero(exitSt, E.exitSt.capacity());
  par::forDyn(&cnt[C_STREAM], nullptr, streamCapI, GEM_LAMBDA(long i) {
    const Particle& p = stream[i];
    if (p.status == 0 && isIon(p.sp) && p.sp != SPECIES_TRAPPED_ION)
      par::atomicAddD(&exitSt[size_t(ionIndex(p.sp)) * zs + p.Zi], p.w);
  });
  {
    const double ra = sim.cfg.grid.r_absorber, rm = sim.g.rho[N], str = sim.cfg.grid.absorber_strength;
    if (ra < rm) {
      const GridV g = P.g;
      // flying focus: R relaxes to the incident field Re(G e^{-i kappa s1}) (m = 1), see Simulation::absorber
      const bool ff = sim.cfg.flying_focus.present;
      const cplx* Gf = ff ? E.ffGf.data() : nullptr;
      const cplx* Gc = ff ? E.ffGc.data() : nullptr;
      const double kph = ff ? -sim.cfg.flying_focus.kappa() * s1 : 0.0;
      const double cph = cos(kph), sph = sin(kph);
      par::forN(long(M + 1) * Nz * (N + 1), GEM_LAMBDA(long idx) {
        const int i = int(idx % (N + 1));
        const long mn = idx / (N + 1);
        const int m = int(mn / Nz), nn = int(mn % Nz);
        auto dmp = [&](double r) {
          if (r <= ra) return 1.0;
          const double x = (r - ra) / (rm - ra);
          return fmax(0.0, 1.0 - str * x * x);
        };
        const double df = dmp(g.rho[i]);
        if (df < 1.0) {
          double inc = 0.0;
          if (ff && m == 1) { const cplx G = Gf[size_t(nn) * (N + 1) + i]; inc = G.re * cph - G.im * sph; }
          const cplx R = fn.Rr[fn.f(m, nn, i)];
          fn.Rr[fn.f(m, nn, i)] = cplx(inc + df * (R.re - inc), df * R.im);
        }
        if (i < N) {
          const double dc = dmp(g.c[i]);
          if (dc < 1.0) {
            double inc = 0.0;  // incident R_phi,1 = i inc
            if (ff && m == 1) { const cplx G = Gc[size_t(nn) * N + i]; inc = G.re * cph - G.im * sph; }
            const cplx R = fn.Rp[fn.c(m, nn, i)];
            fn.Rp[fn.c(m, nn, i)] = cplx(dc * R.re, inc + dc * (R.im - inc));
          }
        }
      });
    }
  }
  // ---- 6. counters, swap levels
  E.pullCounters();
  E.nBeam = E.hcnt[C_BEAM];
  // slow trapped particles continue streaming in the next step
  E.pprev.swap(E.pnext);
  E.nPprev = E.hcnt[C_PNEXT];
  E.iold = inew;
  sim.s0 = s1;
  ++sim.step_no;
  prof.stop(Profile::END);
}

}  // namespace gem
