// ============================================================================
//  simulation_push.cpp -- beam push, fresh plasma, streaming slices, ionization
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
//  particle pushers
// ============================================================================
// beam list: leapfrog along s with sub-steps; fields of level s^0. A particle
// that becomes slower than vz_min stops at its current (z, zeta) and joins the
// streaming particles of the current march at its zeta.
void Simulation::pushBeam() {
  const int nsub = std::max(1, int(std::ceil(Delta / cfg.numerics.beam_dz_max)));
  const double dz = Delta / nsub;
  const double rmax = g.rho[g.N];
  for (auto& p : beam) {
    for (int sub = 0; sub < nsub && p.status == 0; ++sub) {
      double gm = gammaOf(p);
      double vz = p.uz / gm;
      if (vz < cfg.numerics.vz_min) { p.status = 3; ++n_restreamed; break; }
      double dt = dz / vz;
      double w;
      const int na = zv.cell(p.zeta, w);
      if (na < 0 || na >= Nz) { p.status = 2; break; }
      EB f = gather(fold, na, 1.0 - w, na + 1, w, p.x, p.y);
      addExternal(ext, p.x, p.y, p.z, f, dz);
      if (cfg.numerics.beam_vay) vay(p.ux, p.uy, p.uz, f, p.qm, dt);
      else boris(p.ux, p.uy, p.uz, f, p.qm, dt);
      gm = gammaOf(p);
      p.x += p.ux / gm * dt; p.y += p.uy / gm * dt;
      p.zeta += (p.uz / gm - 1.0) * dt;
      p.z += dz;
      if (p.zeta < 0.0 || std::hypot(p.x, p.y) >= rmax) p.status = 2;
    }
  }
  for (auto& p : beam)
    if (p.status == 3) {
      Particle q = p;
      q.status = 0;
      pending.push_back(q);
    }
  beam.erase(std::remove_if(beam.begin(), beam.end(), [](const Particle& p) { return p.status != 0; }), beam.end());
}

// fresh plasma: rings at the front of the window, quiet start in r^2.
// At every position one free electron (if plasma.free_electrons) and one
// macroparticle per ion block are loaded; an ion macroparticle represents
// `atoms` times as many atoms as the electron represents electrons and has
// the charge Z * w, so the loaded plasma is neutral.
void Simulation::loadFresh() {
  stream.clear();
  const auto& P = cfg.plasma;
  if (!P.present) return;
  std::uniform_real_distribution<double> ud(0.0, 1.0);
  const double zmid = s0 + 0.5 * Delta;
  const unsigned long long idBase = (unsigned long long)(step_no + 1) << 32;  // key of the ionization RNG
  // ion blocks with their own profile are loaded wherever their density is positive
  double rlim = P.r_max;
  for (const auto& I : cfg.ions)
    if (I.own) rlim = std::max(rlim, I.prof.r_max);
  for (int i = 0; i < g.N; ++i) {
    if (g.rho[i] >= rlim) break;
    double r2a = g.rho[i] * g.rho[i], r2b = g.rho[i + 1] * g.rho[i + 1];
    const bool ax = (g.rho[i + 1] <= P.r_axis_refine);
    const int npr = ax ? std::max(P.ppc_r, P.ppc_r_axis) : P.ppc_r;
    const int nph = ax ? std::max(P.ppc_phi, P.ppc_phi_axis) : P.ppc_phi;
    for (int j = 0; j < npr; ++j) {
      double r = std::sqrt(r2a + (j + 0.5) / npr * (r2b - r2a));
      const double nAx = P.density(r, zmid);  // axisymmetric part n0 f(z) g(r)
      if (nAx <= 0.0) {  // nothing to load here unless an own ion profile is positive (axisymmetric part)
        bool need = false;
        for (const auto& I : cfg.ions) need = need || (I.own && I.prof.density(r, zmid) > 0.0);
        if (!need) continue;
      }
      const double Q0 = -PI * (r2b - r2a) * Delta / (npr * nph);  // free-electron macro charge per unit density
      double ph0 = 2.0 * PI * ud(rng);
      for (int k = 0; k < nph; ++k) {
        double ph = ph0 + 2.0 * PI * k / nph;
        Particle p{};
        p.x = r * std::cos(ph); p.y = r * std::sin(ph); p.z = s0;
        p.qm = -1.0; p.zeta = zv.node(Nz); p.status = 0; p.sp = SPECIES_PLASMA;
        // density_function(x, y, z, r): the weight of each particle (may break the axial symmetry)
        const double n = (nAx <= 0.0) ? 0.0 : P.fn.empty() ? nAx : P.density(p.x, p.y, zmid);
        const double Q = n * Q0;
        if (n > 0.0) {
          p.q = Q; p.w = -Q;
          if (P.free_electrons) { p.id = idBase + stream.size(); stream.push_back(p); }
        }
        for (size_t s = 0; s < cfg.ions.size(); ++s) {
          const auto& I = cfg.ions[s];
          Particle ion = p;
          if (I.own) {
            const double na = I.prof.density(p.x, p.y, zmid);  // atom density
            if (na <= 0.0) continue;
            ion.w = -na * Q0;
          } else {
            if (n <= 0.0) continue;
            ion.w = -Q * I.atoms;
          }
          if (I.own && I.Z > 0) {  // the electrons of the pre-ionized atoms (before the ions, as above)
            Particle e = p;
            e.q = -I.Z * ion.w; e.w = I.Z * ion.w;
            e.id = idBase + stream.size();
            stream.push_back(e);
          }
          ion.Zi = I.Z;
          ion.q = I.Z * ion.w;
          ion.qm = I.qm();
          ion.sp = SPECIES_ION_BASE - int(s);
          ion.id = idBase + stream.size();
          stream.push_back(ion);
        }
      }
    }
  }
}

// streaming particle from its zeta (normally zeta_{n+1}; later for particles
// that join the march inside the half-cell) to zeta_n (GT2, dt = hl/(1 - v_z))
void Simulation::pushStreamParticle(Particle& p, int n, double s1) { pushStreamParticle(p, n, s1, c1s, n, born); }

void Simulation::pushStreamSlice(int n, double s1, int threads) {
  auto one = [&](Particle& p, Currents& C, int kc, std::vector<Particle>& bo) {
    if (p.status == STATUS_WAITING) {
      if (p.zeta > zv.node(n)) p.status = 0;
      else return;
    }
    pushStreamParticle(p, n, s1, C, kc, bo);
  };
#if defined(_OPENMP)
  if (threads > 1 && stream.size() >= size_t(4 * threads)) {
    // contiguous chunks in thread order: the freed electrons keep the serial order
    const size_t np = stream.size();
    std::vector<Currents> Cb(static_cast<size_t>(threads));
    std::vector<std::vector<Particle>> bo(static_cast<size_t>(threads));
#pragma omp parallel num_threads(threads)
    {
      const int t = omp_get_thread_num(), nt = omp_get_num_threads();
      Currents& C = Cb[size_t(t)];
      C.init(M, 1, g.N);
      const size_t i0 = np * size_t(t) / size_t(nt), i1 = np * size_t(t + 1) / size_t(nt);
      for (size_t i = i0; i < i1; ++i) one(stream[i], C, 0, bo[size_t(t)]);
    }
    for (int t = 0; t < threads; ++t) {
      const Currents& C = Cb[size_t(t)];
      if (C.jr.empty()) continue;
      for (int m = 0; m <= M; ++m) {
        for (int i = 0; i <= g.N; ++i) c1s.jr[c1s.f(m, n, i)] += C.jr[C.f(m, 0, i)];
        for (int i = 0; i < g.N; ++i) {
          c1s.jp[c1s.c(m, n, i)] += C.jp[C.c(m, 0, i)];
          c1s.jz[c1s.c(m, n, i)] += C.jz[C.c(m, 0, i)];
        }
      }
      for (int i = 0; i < g.N; ++i) {
        c1s.rho[size_t(n) * g.N + i] += C.rho[size_t(i)];
        c1s.rhoi[size_t(n) * g.N + i] += C.rhoi[size_t(i)];
      }
      born.insert(born.end(), bo[size_t(t)].begin(), bo[size_t(t)].end());
    }
    return;
  }
#else
  (void)threads;
#endif
  for (auto& p : stream) one(p, c1s, n, born);
}

void Simulation::pushStreamParticle(Particle& p, int n, double s1, Currents& C, int kc,
                                    std::vector<Particle>& bornOut) {
  const bool ion = isIon(p.sp);
  const IonConfig* I = (ion && p.sp != SPECIES_TRAPPED_ION) ? &cfg.ions[ionIndex(p.sp)] : nullptr;
  // length of the path in zeta within this half-cell
  const double h = zv.len(n), zn = zv.node(n);
  const double hl = (p.zeta >= zv.node(n + 1) - 1e-9 * h) ? h : std::max(0.0, p.zeta - zn);
  const double wm = 0.5 * hl / h;  // weight of slice n+1 at the midpoint of the path
  if (I && !I->mobile) {
    // fixed ion: no motion (dt = h), fields only needed for ionization
    deposit(C, n, 1.0, p.x, p.y, 0.0, 0.0, 0.0, p.q * hl / Delta, true, kc);  // density diagnostic only
    const double zs = p.zeta;
    p.zeta = zn;
    if (I->ionization && p.Zi < I->Zmax && sliceFieldBound >= adk[ionIndex(p.sp)].minField()) {
      const double w1 = (zs - zn) / h;  // field at the start of the path
      const EB f1 = gather(fnew, n + 1, w1, n, 1.0 - w1, p.x, p.y);
      const EB f0 = gather(fnew, n, 1.0, n, 0.0, p.x, p.y);
      ionize(p, std::hypot(f1.Ex, f1.Ey, f1.Ez), std::hypot(f0.Ex, f0.Ey, f0.Ez), hl, bornOut);
    }
    return;
  }
  double g0 = gammaOf(p);
  double omv0 = (g0 - p.uz) / g0;  // 1 - v_z > 0
  double dt0 = hl / omv0;
  double xm = p.x + 0.5 * dt0 * p.ux / g0, ym = p.y + 0.5 * dt0 * p.uy / g0;
  // ionizing ions need the field at both slices (the rate is averaged, see
  // ionize); all other particles use the averaged field directly
  const bool ionizing = I && I->ionization && p.Zi < I->Zmax && sliceFieldBound >= adk[ionIndex(p.sp)].minField();
  EB f, f1{}, f0{};
  if (ionizing) {
    const double w1 = 2.0 * wm;  // start of the path
    f1 = gather(fnew, n + 1, w1, n, 1.0 - w1, xm, ym);
    f0 = gather(fnew, n, 1.0, n, 0.0, xm, ym);
    f = {0.5 * (f1.Ex + f0.Ex), 0.5 * (f1.Ey + f0.Ey), 0.5 * (f1.Ez + f0.Ez),
         0.5 * (f1.Bx + f0.Bx), 0.5 * (f1.By + f0.By), 0.5 * (f1.Bz + f0.Bz)};
  } else {
    f = gather(fnew, n + 1, wm, n, 1.0 - wm, xm, ym);
  }
  if (ext.any()) addExternal(ext, xm, ym, p.z + 0.5 * dt0 * p.uz / g0, f, dt0 * p.uz / g0);
  double ux = p.ux, uy = p.uy, uz = p.uz;
  boris(ux, uy, uz, f, p.qm, dt0);
  double g1 = std::sqrt(1.0 + ux * ux + uy * uy + uz * uz);
  double dt = 2.0 * hl / (omv0 + (g1 - uz) / g1);  // time for |dzeta| = hl
  ux = p.ux; uy = p.uy; uz = p.uz;
  boris(ux, uy, uz, f, p.qm, dt);
  g1 = std::sqrt(1.0 + ux * ux + uy * uy + uz * uz);
  double vx = 0.5 * (p.ux / g0 + ux / g1), vy = 0.5 * (p.uy / g0 + uy / g1), vz = 0.5 * (p.uz / g0 + uz / g1);
  double zeta_end = zn;
  const double sTrap = s1 + (cfg.numerics.trap_advance - 1.0) * Delta;
  if (p.z + vz * dt >= sTrap) {  // advanced trap_advance steps inside this half-cell -> trapped
    double dtp = std::max(0.0, (sTrap - p.z) / std::max(vz, 1e-3));
    ux = p.ux; uy = p.uy; uz = p.uz;
    boris(ux, uy, uz, f, p.qm, dtp);
    g1 = std::sqrt(1.0 + ux * ux + uy * uy + uz * uz);
    vx = 0.5 * (p.ux / g0 + ux / g1); vy = 0.5 * (p.uy / g0 + uy / g1);
    vz = 0.5 * (p.uz / g0 + uz / g1);
    dt = dtp;
    zeta_end = p.zeta - (1.0 - vz) * dtp;
    p.status = 1;
  }
  double xo = p.x, yo = p.y;
  p.x += vx * dt; p.y += vy * dt;
  p.z = (p.status == 1) ? s1 : p.z + vz * dt;
  p.ux = ux; p.uy = uy; p.uz = uz;
  p.zeta = zeta_end;
  // streaming current: W = q dt / Delta, into half-cell n
  deposit(C, n, 1.0, 0.5 * (xo + p.x), 0.5 * (yo + p.y), vx, vy, vz, p.q * dt / Delta, ion, kc);
  if (std::hypot(p.x, p.y) >= g.rho[g.N]) p.status = 2;
  // field ionization with the mid-step field (after the deposit: the new
  // charge state acts from the next half-cell on)
  if (ionizing && p.status == 0)
    ionize(p, std::hypot(f1.Ex, f1.Ey, f1.Ez), std::hypot(f0.Ex, f0.Ey, f0.Ez), dt, bornOut);
}

// Upper bound of |E| anywhere on slice n: sum over the modes of the largest
// nodal amplitudes, times 2 for the linear extrapolation to the axis used by
// the gather. Used to skip ionization where no atom can be ionized.
double Simulation::fieldBound(const Fields& F, int n) const {
  double bound = 0.0;
  for (int m = 0; m <= M; ++m) {
    double er = 0, ep = 0, ez = 0;
    for (int i = 0; i <= g.N; ++i) er = std::max(er, 0.5 * std::abs(F.Rr[F.f(m, n, i)] + F.Lr[F.f(m, n, i)]));
    for (int i = 0; i < g.N; ++i) {
      ep = std::max(ep, 0.5 * std::abs(F.Rp[F.c(m, n, i)] + F.Lp[F.c(m, n, i)]));
      ez = std::max(ez, std::abs(F.Ez[F.c(m, n, i)]));
    }
    bound += (m == 0 ? 1.0 : 2.0) * std::sqrt(er * er + ep * ep + ez * ez);
  }
  return 2.0 * bound;
}

// ----------------------------------------------------------------------------
// Monte Carlo field ionization of one ion macroparticle during the time dt.
// E1, E0 = |E| at the two slices of the step. The rate of each level is the
// average (w(E1) + w(E0))/2, i.e. the trapezoidal rule for int w dt. (The rate
// of the averaged field, w((E1+E0)/2), would underestimate the ionization:
// w is very convex and the carrier peaks span only a few slices; for H near
// threshold at 25 points per wavelength the error is -6 %, with the averaged
// rate it is below 1e-4.) Sequential ionization k -> k+1 -> ... is sampled
// with exponential waiting times, so several levels can be ionized in one
// step; for rates constant over the step this reproduces the rate equations
// dN_k/dt = w_{k-1} N_{k-1} - w_k N_k exactly (in the mean).
// Every event raises the charge state of the whole macroparticle by one and
// creates an electron macroparticle (charge -w) at rest relative to the ion.
// ----------------------------------------------------------------------------
void Simulation::ionize(Particle& p, double E1, double E0, double dt, std::vector<Particle>& bornOut) {
  const int s = ionIndex(p.sp);
  const AdkRates& R = adk[s];
  if (E1 < R.minField() && E0 < R.minField()) return;
  const uint64_t base = ionKeyParticle(ionKey, p.id);  // counter-based (rng.hpp)
  unsigned draw = 0;
  double tleft = dt;
  const int Z0 = p.Zi;
  while (p.Zi < R.zmax()) {
    const double w = 0.5 * (R.rate(p.Zi, E1) + R.rate(p.Zi, E0));
    if (!(w > 0.0)) break;
    const double tw = -std::log(1.0 - uniform01(base + draw++)) / w;
    if (tw > tleft) break;
    tleft -= tw;
    ++p.Zi;
    Particle e = p;
    e.q = -p.w;
    e.qm = -1.0;
    e.sp = SPECIES_IONIZED_BASE - s;
    e.Zi = 0;
    bornOut.push_back(e);
  }
  if (p.Zi != Z0) {
    p.q = p.Zi * p.w;
    p.qm = p.Zi / cfg.ions[s].mass;
  }
}

}  // namespace gem
