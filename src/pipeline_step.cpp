// ============================================================================
//  pipeline_step.cpp -- one time step of the MPI pipeline
// ============================================================================
#include "pipeline_internal.hpp"

namespace gem {

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

}  // namespace gem
