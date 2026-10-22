// ============================================================================
//  simulation_step.cpp -- the time step, chicanes, absorber, flying focus
// ============================================================================
#include "simulation.hpp"

#include "eigen.hpp"
#include "rng.hpp"

#include <algorithm>
#if defined(_OPENMP)
#include <omp.h>
#endif
#include <cmath>

namespace gem {

static const cd I1(0.0, 1.0);

// ============================================================================
//  one step s^0 -> s^1
// ============================================================================
void Simulation::beginStep() {
  Delta = cfg.grid.stepAt(s0);
  solver.prepare(Delta, zg);
}

// thin chicanes passed in the last step (z_a, z_b]: the beam particles move by
// r56 delta - delay in zeta (delta = gamma / gamma_ref - 1). Their co-moving
// self-field moves with them (subtracted at the old, added at the new
// positions; vacuum), and the beam current of the level s^0 is deposited anew.
void Simulation::applyChicanes(double za, double zb) {
  for (const auto& C : cfg.chicanes) {
    if (!(za < C.z && C.z <= zb)) continue;
    if (engine) engine->downloadFull();
    addBeamFields(-1.0, false);
    double sh = 0.0, sh2 = 0.0, w = 0.0;
    const double top = zv.node(Nz);
    for (auto& p : beam) {
      const double d = C.r56 * (gammaOf(p) / C.gamma_ref - 1.0) - C.delay;
      p.zeta += d;
      sh += std::fabs(p.q) * d; sh2 += std::fabs(p.q) * d * d; w += std::fabs(p.q);
      if (p.zeta < 0.0 || p.zeta >= top) p.status = 2;
    }
    const size_t n0 = beam.size();
    beam.erase(std::remove_if(beam.begin(), beam.end(), [](const Particle& p) { return p.status != 0; }), beam.end());
    addBeamFields(1.0, false);
    c0b.zero();
    for (const auto& p : beam) depositBeamParticle(c0b, p);
    filterBeamCurrent(c0b);
    if (w > 0) {
      const double um = cfg.units.k0inv_um(), m = sh / w;
      std::printf("chicane at z = %.6g mm: beam moved by %.4g um (rms spread %.3g um)%s\n", C.z * um * 1e-3, m * um,
                  std::sqrt(std::max(0.0, sh2 / w - m * m)) * um,
                  n0 > beam.size() ? (", " + std::to_string(n0 - beam.size()) + " particles left the window").c_str() : "");
    }
    if (engine) engine = std::make_unique<Engine>(*this);  // upload the changed state
  }
}

void Simulation::step() {
  beginStep();
  const double za = s0;
  if (engine) {  // parallel backend: same algorithm, state on the device
    engine->step();
    engine->download();
    if (!cfg.chicanes.empty()) applyChicanes(za, s0);
    return;
  }
  const int N = g.N;
  const double s1 = s0 + Delta;
  c1b.zero(); c1s.zero();
  std::fill(freedCharge.begin(), freedCharge.end(), 0.0);
  pending.swap(pendingNext);  // particles that reached s^0 in the last step with v_z < vz_beam
  pendingNext.clear();
  pushBeam();                 // adds beam particles that became slower than vz_min
  for (const auto& p : beam) depositBeamParticle(c1b, p);
  filterBeamCurrent(c1b);  // (particles trapped during the march are deposited unfiltered)
  loadFresh();
  // particles joining the march at their zeta wait in the streaming list
  for (const auto& p : pending) {
    Particle q = p;
    q.status = STATUS_WAITING;
    stream.push_back(q);
  }
  pending.clear();
  const uint64_t keyStep = ionKeyStep(cfg.numerics.seed, step_no);
  std::fill(predRprev.begin(), predRprev.end(), 0.0);
  std::fill(predPprev.begin(), predPprev.end(), 0.0);

  for (int n = Nz - 1; n >= 0; --n) {
    // prediction of j^1 at half-cell n-1: beam part (known) + streaming part of s^0
    for (int m = 0; m <= M; ++m) {
      for (int i = 0; i <= N; ++i)
        predR[m * (N + 1) + i] = (n > 0) ? c1b.jr[c1b.f(m, n - 1, i)] + c0s.jr[c0s.f(m, n - 1, i)] : 0.0;
      for (int i = 0; i < N; ++i)
        predP[m * N + i] = (n > 0) ? c1b.jp[c1b.c(m, n - 1, i)] + c0s.jp[c0s.c(m, n - 1, i)] : 0.0;
    }
    if (cfg.numerics.iterations > 1) backup = stream;
    for (int it = 0; it < cfg.numerics.iterations; ++it) {
      if (it == 0) fnew.copySlice(fold, n);  // predictor: slow evolution in s
      else stream = backup;
      born.clear();  // electrons freed by ionization: only the last pass counts
      if (!cfg.ions.empty()) sliceFieldBound = std::max(fieldBound(fnew, n + 1), fieldBound(fnew, n));
      c1s.zeroHalf(n);
      ionKey = ionKeySlice(keyStep, n, it);
      for (auto& p : stream) {
        if (p.status == STATUS_WAITING) {  // joins the march inside half-cell n
          if (p.zeta > zv.node(n)) p.status = 0;
          else continue;
        }
        pushStreamParticle(p, n, s1);
      }
      smoothStream(c1s, n);
      // field solve, all modes
      for (int m = 0; m <= M; ++m) {
        // four-point current sums, weighted by the interpolation to slice n
        // (wl = wh = 1 on a uniform grid); at the last slice (n = 0) the missing
        // half-cell -1 is replaced by half-cell 0, at s^1 with the current of
        // this pass (a zero there would halve the source of R(0): a growing spike)
        const int nm = (n > 0) ? n - 1 : 0;
        const double wl = zv.wLo(n), wh = zv.wHi(n);
        const double th1 = 2.0 * cfg.numerics.theta, th0 = 2.0 - th1;  // 1, 1 for the trapezoid
        for (int i = 0; i <= N; ++i) {
          const cd a = c0b.jr[c0b.f(m, nm, i)] + c0s.jr[c0s.f(m, nm, i)];
          const cd pr = (n > 0) ? predR[m * (N + 1) + i] : c1b.jr[c1b.f(m, 0, i)] + c1s.jr[c1s.f(m, 0, i)];
          Jr[i] = th0 * (wl * a + wh * c0b.jr[c0b.f(m, n, i)] + wh * c0s.jr[c0s.f(m, n, i)]) +
                  th1 * (wl * pr + wh * c1b.jr[c1b.f(m, n, i)] + wh * c1s.jr[c1s.f(m, n, i)]);
        }
        for (int i = 0; i < N; ++i) {
          const cd a = c0b.jp[c0b.c(m, nm, i)] + c0s.jp[c0s.c(m, nm, i)];
          const cd pr = (n > 0) ? predP[m * N + i] : c1b.jp[c1b.c(m, 0, i)] + c1s.jp[c1s.c(m, 0, i)];
          Jp[i] = th0 * (wl * a + wh * c0b.jp[c0b.c(m, n, i)] + wh * c0s.jp[c0s.c(m, n, i)]) +
                  th1 * (wl * pr + wh * c1b.jp[c1b.c(m, n, i)] + wh * c1s.jp[c1s.c(m, n, i)]);
        }
        // j^1 at half-cell n = beam + streaming
        for (int i = 0; i <= N; ++i) jr1[i] = c1b.jr[c1b.f(m, n, i)] + c1s.jr[c1s.f(m, n, i)];
        for (int i = 0; i < N; ++i) {
          jp1[i] = c1b.jp[c1b.c(m, n, i)] + c1s.jp[c1s.c(m, n, i)];
          jz1[i] = c1b.jz[c1b.c(m, n, i)] + c1s.jz[c1s.c(m, n, i)];
        }
        solver.solveSlice(m, n, &fold.Rr[fold.f(m, n, 0)], &fold.Rp[fold.c(m, n, 0)], &fold.Ez[fold.c(m, n, 0)],
                          &fold.Bz[fold.f(m, n, 0)], &fnew.Rr[fnew.f(m, n + 1, 0)], &fnew.Lr[fnew.f(m, n + 1, 0)],
                          &fnew.Rp[fnew.c(m, n + 1, 0)], &fnew.Lp[fnew.c(m, n + 1, 0)],
                          &fnew.Ez[fnew.c(m, n + 1, 0)], &fnew.Bz[fnew.f(m, n + 1, 0)], jr1.data(), jp1.data(),
                          jz1.data(), Jr.data(), Jp.data(), &fnew.Rr[fnew.f(m, n, 0)], &fnew.Lr[fnew.f(m, n, 0)],
                          &fnew.Rp[fnew.c(m, n, 0)], &fnew.Lp[fnew.c(m, n, 0)], &fnew.Ez[fnew.c(m, n, 0)],
                          &fnew.Bz[fnew.f(m, n, 0)]);
      }
    }
    // electrons freed in this slice start streaming from slice n
    for (const auto& e : born) freedCharge[ionizedIndex(e.sp)] += e.w;
    stream.insert(stream.end(), born.begin(), born.end());
    born.clear();
    // trapped particles -> beam list; lost particles removed
    for (auto& p : stream) {
      if (p.status == 1) {
        Particle t = p;
        t.status = 0;
        t.z = s1;
        if (t.uz / gammaOf(t) >= cfg.numerics.vz_beam) {
          // fast: continue on the beam list (streaming species get their trapped id)
          if (t.sp == SPECIES_PLASMA) t.sp = SPECIES_TRAPPED;
          else if (isIonizationElectron(t.sp)) t.sp = SPECIES_TRAPPED_IONIZED_BASE - ionizedIndex(t.sp);
          else if (isIon(t.sp)) t.sp = SPECIES_TRAPPED_ION;
          beam.push_back(t);
          depositBeamParticle(c1b, t);
          ++n_trapped_total;
        } else {
          // slow: continue streaming in the next step from its zeta
          pendingNext.push_back(t);
        }
      }
    }
    stream.erase(std::remove_if(stream.begin(), stream.end(),
                                [](const Particle& p) { return p.status != 0 && p.status != STATUS_WAITING; }),
                 stream.end());
    // half-cell n is now final: correct R at slice n+1 for its predicted j^1_{n+1/2}
    if (n + 1 <= Nz - 1) {
      const double q = 0.5 * cfg.numerics.theta * Delta * zv.wLo(n + 1);
      for (int m = 0; m <= M; ++m) {
        for (int i = 1; i <= N; ++i) {
          cd jt = c1b.jr[c1b.f(m, n, i)] + c1s.jr[c1s.f(m, n, i)];
          fnew.Rr[fnew.f(m, n + 1, i)] -= q * (jt - predRprev[m * (N + 1) + i]);
        }
        for (int i = 0; i < N; ++i) {
          cd jt = c1b.jp[c1b.c(m, n, i)] + c1s.jp[c1s.c(m, n, i)];
          fnew.Rp[fnew.c(m, n + 1, i)] -= q * (jt - predPprev[m * N + i]);
        }
      }
    }
    predRprev.swap(predR);
    predPprev.swap(predP);
  }
  // charge states of the ions leaving the rear of the window
  for (auto& v : exitStates) std::fill(v.begin(), v.end(), 0.0);
  for (const auto& p : stream)
    if (p.status == 0 && isIon(p.sp) && p.sp != SPECIES_TRAPPED_ION) exitStates[ionIndex(p.sp)][p.Zi] += p.w;
  absorber(s1);
  std::swap(fold, fnew);
  std::swap(c0b, c1b);
  std::swap(c0s, c1s);
  s0 = s1;
  ++step_no;
  if (!cfg.chicanes.empty()) applyChicanes(za, s0);
}

// damp R near the outer radial boundary; with a flying focus only the difference to the incident field
// (total-field / scattered-field: the layer feeds the converging wave and absorbs what the plasma adds)
void Simulation::absorber(double s) {
  const double ra = cfg.grid.r_absorber, rm = g.rho[g.N];
  if (ra >= rm) return;
  auto damp = [&](double r) {
    if (r <= ra) return 1.0;
    double x = (r - ra) / (rm - ra);
    return std::max(0.0, 1.0 - cfg.grid.absorber_strength * x * x);
  };
  const bool ff = cfg.flying_focus.present;
  for (int m = 0; m <= M; ++m)
    for (int n = 0; n < Nz; ++n) {
      for (int i = 0; i <= g.N; ++i) {
        const double d = damp(g.rho[i]);
        if (d == 1.0) continue;
        const cd inc = (ff && m == 1) ? cd(ffInc(true, n, i, s)) : cd(0.0);
        cd& R = fnew.Rr[fnew.f(m, n, i)];
        R = inc + d * (R - inc);
      }
      for (int i = 0; i < g.N; ++i) {
        const double d = damp(g.c[i]);
        if (d == 1.0) continue;
        const cd inc = (ff && m == 1) ? I1 * ffInc(false, n, i, s) : cd(0.0);
        cd& R = fnew.Rp[fnew.c(m, n, i)];
        R = inc + d * (R - inc);
      }
    }
}

// incident field of the flying focus at s = 0: G(r, zeta) = sum_j c_j J0(k_j r) e^{i w_j (zeta - zeta_c) + i cep},
// c_j the Gaussian spectrum of the on-axis envelope, k_j = sqrt(2 w_j kappa - kappa^2) (k_z = w_j - kappa)
void Simulation::setupFlyingFocus() {
  const auto& F = cfg.flying_focus;
  const double kap = F.kappa(), Lz = F.length_1e, amp = F.a0;  // a0 omega0, omega0 = 1
  const int nw = 241;
  const double w0 = 1.0 - 6.0 / Lz, w1 = 1.0 + 6.0 / Lz, dw = (w1 - w0) / (nw - 1);
  std::vector<double> wj(nw), kj(nw);
  std::vector<cd> cj(nw);
  for (int j = 0; j < nw; ++j) {
    wj[j] = w0 + dw * j;
    kj[j] = std::sqrt(std::max(0.0, 2.0 * wj[j] * kap - kap * kap - F.background));
    // int dw (L / (2 sqrt(pi))) e^{-(w - 1)^2 L^2 / 4} e^{i w d} = e^{i d} e^{-d^2 / L^2}
    cj[j] = amp * dw * Lz / (2.0 * std::sqrt(PI)) * std::exp(-0.25 * (wj[j] - 1.0) * (wj[j] - 1.0) * Lz * Lz) *
            std::polar(1.0, F.cep);
  }
  auto fill = [&](std::vector<cd>& G, const std::vector<double>& rr, int nr, int i0) {
    G.assign(size_t(Nz) * nr, cd(0.0));
    std::vector<double> J(size_t(nr) * nw);
    for (int i = i0; i < nr; ++i)
      for (int j = 0; j < nw; ++j) J[size_t(i) * nw + j] = std::cyl_bessel_j(0.0, kj[j] * rr[i]);
    std::vector<cd> e(nw);
    for (int n = 0; n < Nz; ++n) {
      const double d = zv.node(n) - F.zeta_center;
      for (int j = 0; j < nw; ++j) e[j] = cj[j] * std::polar(1.0, wj[j] * d);
      for (int i = i0; i < nr; ++i) {
        cd sum = 0.0;
        const double* Ji = &J[size_t(i) * nw];
        for (int j = 0; j < nw; ++j) sum += Ji[j] * e[j];
        G[size_t(n) * nr + i] = sum;
      }
    }
  };
  std::vector<double> rf(g.rho.begin(), g.rho.begin() + g.N + 1), rc(g.c.begin(), g.c.begin() + g.N);
  fill(ffGf, rf, g.N + 1, 1);
  fill(ffGc, rc, g.N, 0);
}

bool Simulation::finite() const {
  for (const auto& v : fold.Rr)
    if (!std::isfinite(v.real()) || !std::isfinite(v.imag())) return false;
  for (const auto& v : fold.Ez)
    if (!std::isfinite(v.real()) || !std::isfinite(v.imag())) return false;
  return true;
}

}  // namespace gem
