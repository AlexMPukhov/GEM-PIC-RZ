// ============================================================================
//  engine.cpp -- kernels of the parallel engine (see engine.hpp)
//
//  Every function below that is marked GEM_HD is a line-by-line port of the
//  corresponding member of Simulation / FieldSolver (simulation.cpp,
//  solver.cpp). Differences to the serial code:
//    * deposition with atomic additions (order of the sums not fixed);
//    * random numbers of the field ionization from a counter-based generator,
//      keyed by (seed, particle id, step, slice, pass, draw): statistically
//      equivalent to the serial code, independent of the thread schedule;
//    * particles that join the march later (re-streaming) wait in the
//      streaming list with status WAITING until the march reaches their zeta;
//    * trapped / lost streaming particles are marked, not erased (the
//      streaming list is rebuilt every step anyway).
// ============================================================================
#include "engine.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <stdexcept>
#include <vector>

#include "portable.hpp"
#include "rng.hpp"
#include "simulation.hpp"

namespace gem {
namespace dev {

using par::Team;

constexpr int WAITING = STATUS_WAITING;

// counters in device memory
enum { C_STREAM, C_BORN, C_BEAM, C_TMP, C_PNEXT, C_PPREV, C_OVERFLOW, C_COUNT };
enum { L_RESTREAM, L_TRAPPED, L_COUNT };

// ---------------------------------------------------------------- views
struct W {
  int a;
  double wa, wb;
  int mode;  // 0: between a and a+1;  1: inside r < c_0 (centres);  2: outside
};

struct GridV {
  int N;
  const double *rho, *c, *Af, *Ac, *Vf, *Vc, *kC, *kF;
  // number of entries of x[0..n) that are <= r (std::upper_bound)
  GEM_HD static int upper(const double* x, int n, double r) {
    int lo = 0, hi = n;
    while (lo < hi) {
      const int mid = (lo + hi) >> 1;
      if (x[mid] > r) hi = mid; else lo = mid + 1;
    }
    return lo;
  }
  GEM_HD W locC(double r) const {
    if (r < c[0]) return {0, 1.0, 0.0, 1};
    if (r >= c[N - 1]) return {N - 1, 1.0, 0.0, 2};
    const int a = upper(c, N, r) - 1;
    const double x = (r * r - c[a] * c[a]) / (c[a + 1] * c[a + 1] - c[a] * c[a]);
    const double wb = x + kC[a] * x * (1.0 - x);
    return {a, 1.0 - wb, wb, 0};
  }
  GEM_HD W locF(double r) const {
    if (r >= rho[N]) return {N, 1.0, 0.0, 2};
    const int a = upper(rho, N + 1, r) - 1;
    const double x = (r * r - rho[a] * rho[a]) / (rho[a + 1] * rho[a + 1] - rho[a] * rho[a]);
    const double wb = x + kF[a] * x * (1.0 - x);
    return {a, 1.0 - wb, wb, 0};
  }
  GEM_HD W locC_lin(double r) const {
    if (r < c[0]) return {0, 1.0, 0.0, 1};
    if (r >= c[N - 1]) return {N - 1, 1.0, 0.0, 2};
    const int a = upper(c, N, r) - 1;
    const double wb = (r - c[a]) / (c[a + 1] - c[a]);
    return {a, 1.0 - wb, wb, 0};
  }
  GEM_HD W locF_lin(double r) const {
    if (r >= rho[N]) return {N, 1.0, 0.0, 2};
    const int a = upper(rho, N + 1, r) - 1;
    const double wb = (r - rho[a]) / (rho[a + 1] - rho[a]);
    return {a, 1.0 - wb, wb, 0};
  }
};

struct FieldsV {
  int M, Nz, N;
  cplx *Rr, *Lr, *Bz, *Rp, *Lp, *Ez;
  GEM_HD size_t f(int m, int n, int i) const { return (size_t(m) * (Nz + 1) + n) * (N + 1) + i; }
  GEM_HD size_t c(int m, int n, int i) const { return (size_t(m) * (Nz + 1) + n) * N + i; }
};

struct CurrV {
  int M, Nz, N;
  cplx *jr, *jp, *jz;
  double *rho, *rhoi;
  GEM_HD size_t f(int m, int k, int i) const { return (size_t(m) * Nz + k) * (N + 1) + i; }
  GEM_HD size_t c(int m, int k, int i) const { return (size_t(m) * Nz + k) * N + i; }
};

struct LevelV {
  double logpre, power, expo, Epeak, Emin;
};
struct AdkV {
  int Z0, Zmax;
  double Emin_all;
  const LevelV* lev;  // index k (charge state before the ionization)
};
struct IonV {
  int mobile, ionization, Zmax;
  double mass;
};
struct TriV {
  int i0, n;
  const double *a, *cp, *den;
};

struct Params {
  int M, Nz, N;
  double Delta;  // step of the current s^0 -> s^1
  ZetaView zv;   // longitudinal grid
  const int* cls;  // matrix class of slice n (FieldSolver::classOf)
  ExternalV ext;   // external fields (undulators)
  GridV g;
  int nIons;
  const IonV* ions;
  const AdkV* adk;
  double vz_min, vz_beam;
  unsigned long long seed;
};

// ---------------------------------------------------------------- ADK rate
GEM_HD inline double adkRaw(const LevelV& L, double E) {
  if (!(E > 0)) return 0.0;
  E = fmin(E, L.Epeak);
  return exp(L.logpre - L.power * log(E) + L.expo / E);
}
GEM_HD inline double adkRate(const AdkV& R, int k, double E) {
  if (k < R.Z0 || k >= R.Zmax) return 0.0;
  const LevelV& L = R.lev[k];
  if (E < L.Emin) return 0.0;
  return adkRaw(L, E);
}

// ---------------------------------------------------------------- gather
// (Simulation::accumulate)
GEM_HD inline void accumulate(const GridV& g, const FieldsV& F, int m, int n, double w, double r, const W& wc,
                              const W& wf, cplx& er, cplx& ep, cplx& ez, cplx& br, cplx& bp, cplx& bz) {
  const int N = g.N;
  const size_t cb = F.c(m, n, 0), fb = F.f(m, n, 0);
  auto EpC = [&](int i) { return 0.5 * (F.Rp[cb + i] + F.Lp[cb + i]); };
  auto BrC = [&](int i) { return 0.5 * (F.Lp[cb + i] - F.Rp[cb + i]); };
  auto EzC = [&](int i) { return F.Ez[cb + i]; };
  const double c0 = g.c[0], c1 = g.c[1];
  auto axis = [&](cplx f0, cplx f1) { return (c1 * f0 - c0 * f1) / (c1 - c0); };
  cplx ep0, br0;
  if (m == 1) { ep0 = axis(EpC(0), EpC(1)); br0 = axis(BrC(0), BrC(1)); }

  cplx vep, vbr, vez;
  if (wc.mode == 0) {
    vep = wc.wa * EpC(wc.a) + wc.wb * EpC(wc.a + 1);
    vbr = wc.wa * BrC(wc.a) + wc.wb * BrC(wc.a + 1);
    vez = wc.wa * EzC(wc.a) + wc.wb * EzC(wc.a + 1);
  } else if (wc.mode == 2) {
    vep = EpC(N - 1); vbr = BrC(N - 1); vez = EzC(N - 1);
  } else {
    const double x = r / c0;
    vez = EzC(0) * ((m == 0) ? 1.0 : pow(x, double(m)));
    if (m == 0) { vep = EpC(0) * x; vbr = BrC(0) * x; }
    else if (m == 1) { vep = ep0 + (EpC(0) - ep0) * x; vbr = br0 + (BrC(0) - br0) * x; }
    else { const double fx = pow(x, double(m - 1)); vep = EpC(0) * fx; vbr = BrC(0) * fx; }
  }
  // faces (face 0 for m = 1 from regularity): E_r(0) = -i ep0, B_phi(0) = i br0
  auto ErF = [&](int i) -> cplx {
    if (i == 0) return (m == 1) ? cplx(ep0.im, -ep0.re) : cplx();
    return 0.5 * (F.Rr[fb + i] + F.Lr[fb + i]);
  };
  auto BpF = [&](int i) -> cplx {
    if (i == 0) return (m == 1) ? cplx(-br0.im, br0.re) : cplx();
    return 0.5 * (F.Rr[fb + i] - F.Lr[fb + i]);
  };
  auto BzF = [&](int i) { return F.Bz[fb + i]; };
  cplx ver, vbp, vbz;
  if (wf.mode == 0) {
    ver = wf.wa * ErF(wf.a) + wf.wb * ErF(wf.a + 1);
    vbp = wf.wa * BpF(wf.a) + wf.wb * BpF(wf.a + 1);
    vbz = wf.wa * BzF(wf.a) + wf.wb * BzF(wf.a + 1);
  } else {
    ver = ErF(N); vbp = BpF(N); vbz = BzF(N);
  }
  er += w * ver; ep += w * vep; ez += w * vez; br += w * vbr; bp += w * vbp; bz += w * vbz;
}

// (Simulation::gather)
GEM_HD inline EB gather(const Params& P, const FieldsV& F, int n1, double w1, int n2, double w2, double x, double y) {
  const double r = sqrt(x * x + y * y);
  double cs = 1.0, sn = 0.0;
  if (r > 0.0) { cs = x / r; sn = y / r; }
  const W wc = P.g.locC_lin(r), wf = P.g.locF_lin(r);
  const cplx eip(cs, sn);
  cplx em(1.0, 0.0);
  double Er = 0, Ep = 0, Ez = 0, Br = 0, Bp = 0, Bz = 0;
  for (int m = 0; m <= P.M; ++m) {
    cplx er, ep, ez, br, bp, bz;
    if (w1 != 0.0) accumulate(P.g, F, m, n1, w1, r, wc, wf, er, ep, ez, br, bp, bz);
    if (w2 != 0.0) accumulate(P.g, F, m, n2, w2, r, wc, wf, er, ep, ez, br, bp, bz);
    const double fac = (m == 0) ? 1.0 : 2.0;
    Er += fac * (er * em).re; Ep += fac * (ep * em).re;
    Ez += fac * (ez * em).re; Br += fac * (br * em).re;
    Bp += fac * (bp * em).re; Bz += fac * (bz * em).re;
    em = em * eip;
  }
  return {Er * cs - Ep * sn, Er * sn + Ep * cs, Ez, Br * cs - Bp * sn, Br * sn + Bp * cs, Bz};
}

// |E| as in the serial code (std::hypot of three arguments on the host)
GEM_HD inline double norm3(const EB& f) {
#if defined(__CUDA_ARCH__)
  return norm3d(f.Ex, f.Ey, f.Ez);
#else
  return std::hypot(f.Ex, f.Ey, f.Ez);
#endif
}

// ---------------------------------------------------------------- deposition
// add v to a current / density entry: atomic for shared arrays, plain for the
// private per-thread buffers of the CPU backend
template <bool Atomic> GEM_HD inline void addTo(cplx* a, const cplx& v) {
  if (Atomic) par::atomicAddC(a, v); else *a += v;
}
template <bool Atomic> GEM_HD inline void addTo(double* a, double v) {
  if (Atomic) par::atomicAddD(a, v); else *a += v;
}

// (Simulation::deposit); k = half-cell index in C (C.Nz = 1 for a buffer of one half-cell),
// kz = the half-cell in the window (its length h_kz)
template <bool Atomic>
GEM_HD inline void deposit(const Params& P, const CurrV& C, int k, int kz, double wk, double x, double y, double vx,
                           double vy, double vz, double Wq, bool ion) {
  if (k < 0 || k >= C.Nz || wk == 0.0 || Wq == 0.0) return;
  const GridV& g = P.g;
  const int N = g.N;
  const double r = sqrt(x * x + y * y);
  double cs = 1.0, sn = 0.0;
  if (r > 0.0) { cs = x / r; sn = y / r; }
  const double vr = cs * vx + sn * vy, vp = -sn * vx + cs * vy;
  const W wc = g.locC(r), wf = g.locF(r);
  const double norm = Wq * wk / (2.0 * PI * P.zv.len(kz));
  const cplx emi(cs, -sn);
  cplx em(1.0, 0.0);
  for (int m = 0; m <= P.M; ++m) {
    auto addC = [&](int i, double sp, double sz) {
      const double f = norm / g.Vc[i];
      addTo<Atomic>(&C.jp[C.c(m, k, i)], f * sp * vp * em);
      addTo<Atomic>(&C.jz[C.c(m, k, i)], f * sz * vz * em);
    };
    if (wc.mode == 0) { addC(wc.a, wc.wa, wc.wa); addC(wc.a + 1, wc.wb, wc.wb); }
    else if (wc.mode == 2) addC(N - 1, 1.0, 1.0);
    else addC(0, 1.0, (m == 0) ? 1.0 : pow(r / g.c[0], double(m)));
    auto addF = [&](int i, double s) {
      if (i == 0) return;
      addTo<Atomic>(&C.jr[C.f(m, k, i)], norm * s / g.Vf[i] * vr * em);
    };
    if (wf.mode == 0) { addF(wf.a, wf.wa); addF(wf.a + 1, wf.wb); }
    else addF(N, 1.0);
    em = em * emi;
  }
  double* dens = ion ? C.rhoi : C.rho;
  auto addR = [&](int i, double s) { addTo<Atomic>(&dens[size_t(k) * N + i], norm * s / g.Vc[i]); };
  if (wc.mode == 0) { addR(wc.a, wc.wa); addR(wc.a + 1, wc.wb); }
  else if (wc.mode == 2) addR(N - 1, 1.0);
  else addR(0, 1.0);
}

// (Simulation::depositBeamParticle)
GEM_HD inline void depositBeam(const Params& P, const CurrV& C, const Particle& p) {
  const double gm = gammaOf(p);
  const double vx = p.ux / gm, vy = p.uy / gm, vz = p.uz / gm;
  if (vz <= 0.0) return;
  double f;
  const int j = P.zv.beamCell(p.zeta, f);
  const double Wq = p.q / vz;
  const bool ion = isIon(p.sp);
  deposit<true>(P, C, j - 1, j - 1, 1.0 - f, p.x, p.y, vx, vy, vz, Wq, ion);
  deposit<true>(P, C, j, j, f, p.x, p.y, vx, vy, vz, Wq, ion);
}

// ---------------------------------------------------------------- ionization
struct Ionizer {
  Particle* born;
  int* cnt;
  int bornCap;
  uint64_t key;  // ionKeySlice(...) of the pass (rng.hpp)
};

// (Simulation::ionize)
GEM_HD inline void ionize(const Params& P, Particle& p, double E1, double E0, double dt, const Ionizer& Z) {
  const int s = ionIndex(p.sp);
  const AdkV& R = P.adk[s];
  if (E1 < R.Emin_all && E0 < R.Emin_all) return;
  const uint64_t base = ionKeyParticle(Z.key, p.id);
  double tleft = dt;
  const int Z0 = p.Zi;
  unsigned draw = 0;
  while (p.Zi < R.Zmax) {
    const double w = 0.5 * (adkRate(R, p.Zi, E1) + adkRate(R, p.Zi, E0));
    if (!(w > 0.0)) break;
    const double tw = -log(1.0 - uniform01(base + draw++)) / w;
    if (tw > tleft) break;
    tleft -= tw;
    ++p.Zi;
    Particle e = p;
    e.q = -p.w;
    e.qm = -1.0;
    e.sp = SPECIES_IONIZED_BASE - s;
    e.Zi = 0;
    const int j = par::atomicFetchAdd(&Z.cnt[C_BORN], 1);
    if (j < Z.bornCap) Z.born[j] = e;
    else Z.cnt[C_OVERFLOW] = 1;
  }
  if (p.Zi != Z0) {
    p.q = p.Zi * p.w;
    p.qm = p.Zi / P.ions[s].mass;
  }
}

// ---------------------------------------------------------------- pushers
// (Simulation::pushStreamParticle); F = fnew, C = c1s (deposit into its
// half-cell kdep: n for the shared array, 0 for a one-half-cell buffer)
template <bool Atomic>
GEM_HD inline void pushStream(const Params& P, const FieldsV& F, const CurrV& C, int kdep, Particle& p, int n,
                              double s1, double bound, const Ionizer& Z) {
  const double h = P.zv.len(n), zn = P.zv.node(n);
  const bool ion = isIon(p.sp);
  const IonV* I = (ion && p.sp != SPECIES_TRAPPED_ION) ? &P.ions[ionIndex(p.sp)] : nullptr;
  const double hl = (p.zeta >= P.zv.node(n + 1) - 1e-9 * h) ? h : fmax(0.0, p.zeta - zn);
  const double wm = 0.5 * hl / h;
  if (I && !I->mobile) {
    deposit<Atomic>(P, C, kdep, n, 1.0, p.x, p.y, 0.0, 0.0, 0.0, p.q * hl / P.Delta, true);
    const double zs = p.zeta;
    p.zeta = zn;
    if (I->ionization && p.Zi < I->Zmax && bound >= P.adk[ionIndex(p.sp)].Emin_all) {
      const double w1 = (zs - zn) / h;
      const EB f1 = gather(P, F, n + 1, w1, n, 1.0 - w1, p.x, p.y);
      const EB f0 = gather(P, F, n, 1.0, n, 0.0, p.x, p.y);
      ionize(P, p, norm3(f1), norm3(f0), hl, Z);
    }
    return;
  }
  const double g0 = gammaOf(p);
  const double omv0 = (g0 - p.uz) / g0;
  double dt0 = hl / omv0;
  const double xm = p.x + 0.5 * dt0 * p.ux / g0, ym = p.y + 0.5 * dt0 * p.uy / g0;
  const bool ionizing = I && I->ionization && p.Zi < I->Zmax && bound >= P.adk[ionIndex(p.sp)].Emin_all;
  EB f, f1{}, f0{};
  if (ionizing) {
    const double w1 = 2.0 * wm;
    f1 = gather(P, F, n + 1, w1, n, 1.0 - w1, xm, ym);
    f0 = gather(P, F, n, 1.0, n, 0.0, xm, ym);
    f = {0.5 * (f1.Ex + f0.Ex), 0.5 * (f1.Ey + f0.Ey), 0.5 * (f1.Ez + f0.Ez),
         0.5 * (f1.Bx + f0.Bx), 0.5 * (f1.By + f0.By), 0.5 * (f1.Bz + f0.Bz)};
  } else {
    f = gather(P, F, n + 1, wm, n, 1.0 - wm, xm, ym);
  }
  if (P.ext.n) addExternal(P.ext, xm, ym, p.z + 0.5 * dt0 * p.uz / g0, f);
  double ux = p.ux, uy = p.uy, uz = p.uz;
  boris(ux, uy, uz, f, p.qm, dt0);
  double g1 = sqrt(1.0 + ux * ux + uy * uy + uz * uz);
  double dt = 2.0 * hl / (omv0 + (g1 - uz) / g1);
  ux = p.ux; uy = p.uy; uz = p.uz;
  boris(ux, uy, uz, f, p.qm, dt);
  g1 = sqrt(1.0 + ux * ux + uy * uy + uz * uz);
  double vx = 0.5 * (p.ux / g0 + ux / g1), vy = 0.5 * (p.uy / g0 + uy / g1), vz = 0.5 * (p.uz / g0 + uz / g1);
  double zeta_end = zn;
  if (p.z + vz * dt >= s1) {
    const double dtp = fmax(0.0, (s1 - p.z) / fmax(vz, 1e-3));
    ux = p.ux; uy = p.uy; uz = p.uz;
    boris(ux, uy, uz, f, p.qm, dtp);
    g1 = sqrt(1.0 + ux * ux + uy * uy + uz * uz);
    vx = 0.5 * (p.ux / g0 + ux / g1); vy = 0.5 * (p.uy / g0 + uy / g1);
    vz = 0.5 * (p.uz / g0 + uz / g1);
    dt = dtp;
    zeta_end = p.zeta - (1.0 - vz) * dtp;
    p.status = 1;
  }
  const double xo = p.x, yo = p.y;
  p.x += vx * dt; p.y += vy * dt;
  p.z = (p.status == 1) ? s1 : p.z + vz * dt;
  p.ux = ux; p.uy = uy; p.uz = uz;
  p.zeta = zeta_end;
  deposit<Atomic>(P, C, kdep, n, 1.0, 0.5 * (xo + p.x), 0.5 * (yo + p.y), vx, vy, vz, p.q * dt / P.Delta, ion);
  if (hypot(p.x, p.y) >= P.g.rho[P.g.N]) p.status = 2;
  if (ionizing && p.status == 0) ionize(P, p, norm3(f1), norm3(f0), dt, Z);
}

// one streaming particle in slice n (joins the march if it waits for its zeta)
template <bool Atomic>
GEM_HD inline void streamParticle(const Params& P, const FieldsV& F, const CurrV& C, int kdep, Particle& p, int n,
                                  double s1, double bound, const Ionizer& Z) {
  if (p.status == WAITING) {
    if (p.zeta > P.zv.node(n)) p.status = 0;
    else return;
  }
  if (p.status != 0) return;
  pushStream<Atomic>(P, F, C, kdep, p, n, s1, bound, Z);
}

// (Simulation::pushBeam, one particle); F = fold
GEM_HD inline void pushBeamParticle(const Params& P, const FieldsV& F, Particle& p, int nsub, double dz,
                                    unsigned long long* lcnt) {
  const double rmax = P.g.rho[P.g.N];
  for (int sub = 0; sub < nsub && p.status == 0; ++sub) {
    double gm = gammaOf(p);
    const double vz = p.uz / gm;
    if (vz < P.vz_min) { p.status = 3; par::atomicAddU64(&lcnt[L_RESTREAM], 1ULL); break; }
    const double dt = dz / vz;
    double w;
    const int na = P.zv.cell(p.zeta, w);
    if (na < 0 || na >= P.Nz) { p.status = 2; break; }
    EB f = gather(P, F, na, 1.0 - w, na + 1, w, p.x, p.y);
    addExternal(P.ext, p.x, p.y, p.z, f);
    boris(p.ux, p.uy, p.uz, f, p.qm, dt);
    gm = gammaOf(p);
    p.x += p.ux / gm * dt; p.y += p.uy / gm * dt;
    p.zeta += (p.uz / gm - 1.0) * dt;
    p.z += dz;
    if (p.zeta < 0.0 || hypot(p.x, p.y) >= rmax) p.status = 2;
  }
}

// ---------------------------------------------------------------- field solve
GEM_HD inline void thomas(const TriV& T, const cplx* rhs, cplx* x, cplx* dp) {
  dp[0] = rhs[T.i0] * T.den[0];
  for (int i = 1; i < T.n; ++i) dp[i] = (rhs[T.i0 + i] - T.a[i] * dp[i - 1]) * T.den[i];
  x[T.i0 + T.n - 1] = dp[T.n - 1];
  for (int i = T.n - 2; i >= 0; --i) x[T.i0 + i] = dp[i] - T.cp[i] * x[T.i0 + i + 1];
}

constexpr int SOLVE_SCRATCH = 13;  // arrays of N+1 per mode

// (FieldSolver::solveSlice for mode m, slice n, with the current sums of
//  Simulation::step); one team per mode
GEM_HD inline void solveSliceTeam(const Team& t, const Params& P, int m, int n, const FieldsV& fo, const FieldsV& fn,
                                  const CurrV& c0b, const CurrV& c0s, const CurrV& c1b, const CurrV& c1s,
                                  const cplx* predR, const cplx* predP, const TriV* Tf, const TriV* Tc,
                                  cplx* scratch) {
  const int N = P.N;
  const GridV& g = P.g;
  const double D = P.Delta, h = P.zv.len(n);
  const int T = P.cls[n] * (P.M + 1) + m;  // matrices of this slice and mode
  const cplx im(0.0, double(m));
  const cplx *Rr0 = &fo.Rr[fo.f(m, n, 0)], *Rp0 = &fo.Rp[fo.c(m, n, 0)];
  const cplx *e0 = &fo.Ez[fo.c(m, n, 0)], *b0 = &fo.Bz[fo.f(m, n, 0)];
  const cplx *Rr1 = &fn.Rr[fn.f(m, n + 1, 0)], *Lr1 = &fn.Lr[fn.f(m, n + 1, 0)];
  const cplx *Rp1 = &fn.Rp[fn.c(m, n + 1, 0)], *Lp1 = &fn.Lp[fn.c(m, n + 1, 0)];
  const cplx *e1 = &fn.Ez[fn.c(m, n + 1, 0)], *b1 = &fn.Bz[fn.f(m, n + 1, 0)];
  cplx *Rr = &fn.Rr[fn.f(m, n, 0)], *Lr = &fn.Lr[fn.f(m, n, 0)];
  cplx *Rp = &fn.Rp[fn.c(m, n, 0)], *Lp = &fn.Lp[fn.c(m, n, 0)];
  cplx *e = &fn.Ez[fn.c(m, n, 0)], *b = &fn.Bz[fn.f(m, n, 0)];
  const int L = N + 1;
  cplx *Er1 = scratch, *Bp1 = scratch + L, *Ep1 = scratch + 2 * L, *Br1 = scratch + 3 * L;
  cplx *Pr = scratch + 4 * L, *Qr = scratch + 5 * L, *Pp = scratch + 6 * L, *Qp = scratch + 7 * L;
  cplx *tf1 = scratch + 8 * L, *tc1 = scratch + 9 * L, *rf = scratch + 10 * L, *rc = scratch + 11 * L;
  cplx* dp = scratch + 12 * L;

  // four-point current sums and j^1 (Simulation::step), weighted by the
  // interpolation to slice n (last slice n = 0: half-cell -1 := half-cell 0)
  const int nm = (n > 0) ? n - 1 : 0;
  const double wl = P.zv.wLo(n), wh = P.zv.wHi(n);
  auto Jr = [&](int i) {
    const cplx a = c0b.jr[c0b.f(m, nm, i)] + c0s.jr[c0s.f(m, nm, i)];
    const cplx pr = (n > 0) ? predR[m * (N + 1) + i] : c1b.jr[c1b.f(m, 0, i)] + c1s.jr[c1s.f(m, 0, i)];
    return wl * a + wh * c0b.jr[c0b.f(m, n, i)] + wh * c0s.jr[c0s.f(m, n, i)] + wl * pr +
           wh * c1b.jr[c1b.f(m, n, i)] + wh * c1s.jr[c1s.f(m, n, i)];
  };
  auto Jp = [&](int i) {
    const cplx a = c0b.jp[c0b.c(m, nm, i)] + c0s.jp[c0s.c(m, nm, i)];
    const cplx pr = (n > 0) ? predP[m * N + i] : c1b.jp[c1b.c(m, 0, i)] + c1s.jp[c1s.c(m, 0, i)];
    return wl * a + wh * c0b.jp[c0b.c(m, n, i)] + wh * c0s.jp[c0s.c(m, n, i)] + wl * pr +
           wh * c1b.jp[c1b.c(m, n, i)] + wh * c1s.jp[c1s.c(m, n, i)];
  };
  auto jr1 = [&](int i) { return c1b.jr[c1b.f(m, n, i)] + c1s.jr[c1s.f(m, n, i)]; };
  auto jp1 = [&](int i) { return c1b.jp[c1b.c(m, n, i)] + c1s.jp[c1s.c(m, n, i)]; };
  auto jz1 = [&](int i) { return c1b.jz[c1b.c(m, n, i)] + c1s.jz[c1s.c(m, n, i)]; };

  // discrete operators (FieldSolver::opG ... opCB), one entry each
  auto opG = [&](const cplx* x, int i) -> cplx {
    if (i == 0) return cplx();
    const cplx xi = (i < N) ? x[i] : cplx();
    return (xi - x[i - 1]) / (g.c[i] - g.c[i - 1]);
  };
  auto opg = [&](const cplx* x, int i) { return (x[i + 1] - x[i]) / (g.rho[i + 1] - g.rho[i]); };
  auto opMf = [&](const cplx* x, int i) -> cplx {
    if (i == 0) return cplx();
    return (im / g.rho[i]) * x[i];
  };
  auto opCE = [&](const cplx* Er, const cplx* Ep, int i) -> cplx {
    if (i == 0) return g.c[0] * Ep[0] / g.Af[0];
    if (i == N) return cplx();
    return (g.c[i] * Ep[i] - g.c[i - 1] * Ep[i - 1]) / g.Af[i] - (im / g.rho[i]) * Er[i];
  };
  auto opCB = [&](const cplx* Bp, const cplx* Br, int i) -> cplx {
    const cplx inner = (i == 0) ? cplx() : g.rho[i] * Bp[i];
    return (g.rho[i + 1] * Bp[i + 1] - inner) / g.Ac[i] - (im / g.c[i]) * Br[i];
  };

  par::teamFor(t, N + 1, [&](int i) { Er1[i] = 0.5 * (Rr1[i] + Lr1[i]); Bp1[i] = 0.5 * (Rr1[i] - Lr1[i]); });
  par::teamFor(t, N, [&](int i) { Ep1[i] = 0.5 * (Rp1[i] + Lp1[i]); Br1[i] = 0.5 * (Lp1[i] - Rp1[i]); });
  // known parts P (level s^0) and Q (from zeta_1)
  par::teamFor(t, N + 1, [&](int i) {
    if (i == 0) { Pr[0] = cplx(); Qr[0] = cplx(); return; }
    Pr[i] = Rr0[i] + 0.5 * D * (opG(e0, i) + opMf(b0, i)) - 0.25 * D * Jr(i);
    Qr[i] = Lr1[i] - 0.25 * h * (opG(e1, i) - opMf(b1, i)) - 0.5 * h * jr1(i);
  });
  par::teamFor(t, N, [&](int i) {
    Pp[i] = Rp0[i] + 0.5 * D * ((im / g.c[i]) * e0[i] - opg(b0, i)) - 0.25 * D * Jp(i);
    Qp[i] = Lp1[i] - 0.25 * h * ((im / g.c[i]) * e1[i] + opg(b1, i)) - 0.5 * h * jp1(i);
  });
  t.barrier();
  // B_z
  par::teamFor(t, N + 1, [&](int i) { tf1[i] = 0.5 * (Pr[i] + Qr[i]); });
  par::teamFor(t, N, [&](int i) { tc1[i] = 0.5 * (Pp[i] + Qp[i]); });
  t.barrier();
  par::teamFor(t, N + 1, [&](int i) { rf[i] = 16.0 * b1[i] - 8.0 * h * (opCE(Er1, Ep1, i) + opCE(tf1, tc1, i)); });
  t.barrier();
  if (t.leader()) {
    for (int i = 0; i <= N; ++i) b[i] = cplx();
    thomas(Tf[T], rf, b, dp);
  }
  t.barrier();
  // E_z
  par::teamFor(t, N + 1, [&](int i) { tf1[i] = 0.5 * (Pr[i] - Qr[i]); });
  par::teamFor(t, N, [&](int i) { tc1[i] = 0.5 * (Qp[i] - Pp[i]); });
  t.barrier();
  par::teamFor(t, N, [&](int i) {
    rc[i] = 16.0 * e1[i] + 8.0 * h * (opCB(Bp1, Br1, i) + opCB(tf1, tc1, i)) - 16.0 * h * jz1(i);
  });
  t.barrier();
  if (t.leader()) thomas(Tc[T], rc, e, dp);
  t.barrier();
  // explicit transverse update
  par::teamFor(t, N + 1, [&](int i) {
    if (i == 0) { Rr[0] = cplx(); Lr[0] = cplx(); return; }
    const cplx G = opG(e, i), Mf = opMf(b, i);
    Rr[i] = Pr[i] + 0.5 * D * (G + Mf);
    Lr[i] = Qr[i] - 0.25 * h * (G - Mf);
  });
  par::teamFor(t, N, [&](int i) {
    const cplx gg = opg(b, i), ie = (im / g.c[i]) * e[i];
    Rp[i] = Pp[i] + 0.5 * D * (ie - gg);
    Lp[i] = Qp[i] - 0.25 * h * (ie + gg);
  });
}

// ---------------------------------------------------------------- storage
struct FieldStore {
  par::Array<cplx> Rr, Lr, Bz, Rp, Lp, Ez;
  void init(size_t nf, size_t nc) {
    for (auto* a : {&Rr, &Lr, &Bz}) a->clear(nf);
    for (auto* a : {&Rp, &Lp, &Ez}) a->clear(nc);
  }
  FieldsV view(int M, int Nz, int N) {
    return {M, Nz, N, Rr.data(), Lr.data(), Bz.data(), Rp.data(), Lp.data(), Ez.data()};
  }
};
struct CurrStore {
  par::Array<cplx> jr, jp, jz;
  par::Array<double> rho, rhoi;
  size_t nf = 0, nc = 0, nd = 0;
  void init(size_t nf_, size_t nc_, size_t nd_) {
    nf = nf_; nc = nc_; nd = nd_;
    jr.clear(nf); jp.clear(nc); jz.clear(nc); rho.clear(nd); rhoi.clear(nd);
  }
  void zero() {
    par::zero(jr.data(), nf); par::zero(jp.data(), nc); par::zero(jz.data(), nc);
    par::zero(rho.data(), nd); par::zero(rhoi.data(), nd);
  }
  CurrV view(int M, int Nz, int N) { return {M, Nz, N, jr.data(), jp.data(), jz.data(), rho.data(), rhoi.data()}; }
};

}  // namespace dev

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
  std::vector<dev::CurrStore> local;  // CPU: per-thread currents of one half-cell
  // particles
  par::Array<Particle> stream, backup, born, beam, beamTmp, pnext, pprev;
  par::Array<int> cnt;
  par::Array<unsigned long long> lcnt;
  par::Array<double> bound, freed, exitSt;
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

bool Engine::gpuBuild() {
#if defined(GEM_BACKEND_CUDA)
  return true;
#else
  return false;
#endif
}

std::string Engine::describe() { return par::describe(); }

Engine::~Engine() { st->prof.report(); }

Engine::Engine(Simulation& s) : sim(s), st((par::initDevice(), new EngineState)) {
  EngineState& E = *st;
  if (const char* e = std::getenv("GEMPIC_PROFILE")) E.prof.on = (e[0] == '1');
  const RadialGrid& g = sim.g;
  const int M = sim.M, Nz = sim.Nz, N = g.N;
  E.nf = size_t(M + 1) * (Nz + 1) * (N + 1);
  E.nc = size_t(M + 1) * (Nz + 1) * N;
  // grid
  E.rho.assign(g.rho.data(), g.rho.size());
  E.c.assign(g.c.data(), g.c.size());
  E.Af.assign(g.Af.data(), g.Af.size());
  E.Ac.assign(g.Ac.data(), g.Ac.size());
  E.Vf.assign(g.Vf.data(), g.Vf.size());
  E.Vc.assign(g.Vc.data(), g.Vc.size());
  E.kC.assign(g.kC.data(), g.kC.size());
  E.kF.assign(g.kF.data(), g.kF.size());
  // ions and ADK tables
  const auto& ionsCfg = sim.cfg.ions;
  std::vector<dev::IonV> iv;
  std::vector<dev::LevelV> lv;
  std::vector<size_t> levOff;
  for (size_t s = 0; s < ionsCfg.size(); ++s) {
    const auto& I = ionsCfg[s];
    iv.push_back({I.mobile ? 1 : 0, I.ionization ? 1 : 0, I.Zmax, I.mass});
    levOff.push_back(lv.size());
    const AdkRates& R = sim.adk[s];
    for (const auto& L : R.lev) lv.push_back({L.logpre, L.power, L.expo, L.Epeak, L.Emin});
    E.zStride = std::max(E.zStride, int(sim.exitStates[s].size()));
    if (I.ionization) E.maxSpan = std::max(E.maxSpan, I.Zmax);
  }
  E.levels.assign(lv.data(), lv.size());
  std::vector<dev::AdkV> av;
  for (size_t s = 0; s < ionsCfg.size(); ++s) {
    const AdkRates& R = sim.adk[s];
    av.push_back({R.Z0, R.Zmax, R.Emin_all, E.levels.data() + levOff[s]});
  }
  E.ions.assign(iv.data(), iv.size());
  E.adk.assign(av.data(), av.size());
  // tridiagonal matrices (factorized on the host, per matrix class and mode)
  E.uploadMatrices(sim.solver, M);
  {
    const ZetaGrid& zg = sim.zg;
    E.zz.assign(zg.z.data(), zg.z.size());
    E.zhs.assign(zg.hs.data(), zg.hs.size());
    E.zzc.assign(zg.zc.data(), zg.zc.size());
    E.zwlo.assign(zg.wlo.data(), zg.wlo.size());
    E.zwhi.assign(zg.whi.data(), zg.whi.size());
    std::vector<int> cl(Nz);
    for (int n = 0; n < Nz; ++n) cl[n] = sim.solver.classOf(n);
    E.cls.assign(cl.data(), cl.size());
  }
  // parameters
  dev::Params& P = E.P;
  P.M = M; P.Nz = Nz; P.N = N; P.Delta = sim.Delta;
  P.zv = ZetaView{Nz, sim.zg.uniform ? 1 : 0, sim.zg.h, E.zz.data(), E.zhs.data(), E.zzc.data(), E.zwlo.data(),
                  E.zwhi.data()};
  P.cls = E.cls.data();
  if (!sim.undV.empty()) E.und.assign(sim.undV.data(), sim.undV.size());
  P.ext = ExternalV{int(sim.undV.size()), sim.undV.empty() ? nullptr : E.und.data()};
  P.g = {N, E.rho.data(), E.c.data(), E.Af.data(), E.Ac.data(), E.Vf.data(), E.Vc.data(), E.kC.data(), E.kF.data()};
  P.nIons = int(ionsCfg.size());
  P.ions = E.ions.data();
  P.adk = E.adk.data();
  P.vz_min = sim.cfg.numerics.vz_min;
  P.vz_beam = sim.cfg.numerics.vz_beam;
  P.seed = sim.cfg.numerics.seed;
  // fields: fold from the host (lasers, bunch fields), fnew = 0
  for (auto& F : E.F) F.init(E.nf, E.nc);
  auto up = [](par::Array<cplx>& d, const std::vector<cd>& h) {
    par::upload(d.data(), reinterpret_cast<const cplx*>(h.data()), h.size());
  };
  dev::FieldStore& Fo = E.F[E.iold];
  up(Fo.Rr, sim.fold.Rr); up(Fo.Lr, sim.fold.Lr); up(Fo.Bz, sim.fold.Bz);
  up(Fo.Rp, sim.fold.Rp); up(Fo.Lp, sim.fold.Lp); up(Fo.Ez, sim.fold.Ez);
  const size_t cf = size_t(M + 1) * Nz * (N + 1), cc = size_t(M + 1) * Nz * N;
  E.nd = size_t(Nz) * N;
  for (int k = 0; k < 2; ++k) { E.Cb[k].init(cf, cc, E.nd); E.Cs[k].init(cf, cc, E.nd); }
  up(E.Cb[E.iold].jr, sim.c0b.jr); up(E.Cb[E.iold].jp, sim.c0b.jp); up(E.Cb[E.iold].jz, sim.c0b.jz);
  par::upload(E.Cb[E.iold].rho.data(), sim.c0b.rho.data(), E.nd);
  par::upload(E.Cb[E.iold].rhoi.data(), sim.c0b.rhoi.data(), E.nd);
  // streaming currents of s^0 (zero at the start, from the checkpoint after a restart)
  up(E.Cs[E.iold].jr, sim.c0s.jr); up(E.Cs[E.iold].jp, sim.c0s.jp); up(E.Cs[E.iold].jz, sim.c0s.jz);
  par::upload(E.Cs[E.iold].rho.data(), sim.c0s.rho.data(), E.nd);
  par::upload(E.Cs[E.iold].rhoi.data(), sim.c0s.rhoi.data(), E.nd);
  for (int k = 0; k < 2; ++k) { E.pred[k].clear(size_t(M + 1) * (N + 1)); E.predP[k].clear(size_t(M + 1) * N); }
  E.solveScratch.clear(size_t(M + 1) * dev::SOLVE_SCRATCH * (N + 1));
#if !defined(GEM_BACKEND_CUDA) && !defined(GEM_EMULATE_GPU)
#if defined(_OPENMP)
  E.local = std::vector<dev::CurrStore>(size_t(omp_get_max_threads()));
#else
  E.local = std::vector<dev::CurrStore>(1);
#endif
  for (auto& L : E.local) L.init(size_t(M + 1) * (N + 1), size_t(M + 1) * N, size_t(N));
#endif
  // particles
  E.nBeam = int(sim.beam.size());
  E.beam.reserve(std::max<size_t>(1, sim.beam.size()));
  par::upload(E.beam.data(), sim.beam.data(), sim.beam.size());
  E.cnt.clear(dev::C_COUNT);
  E.lcnt.clear(dev::L_COUNT);
  E.bound.clear(1);
  E.freed.clear(std::max<size_t>(1, ionsCfg.size()));
  E.exitSt.clear(std::max<size_t>(1, ionsCfg.size() * E.zStride));
  E.hcnt[dev::C_BEAM] = E.nBeam;
  E.pushCounters();
  // particles that continue streaming in the next step (after a restart)
  E.nPprev = int(sim.pendingNext.size());
  E.pprev.reserve(std::max<size_t>(1, sim.pendingNext.size()));
  par::upload(E.pprev.data(), sim.pendingNext.data(), sim.pendingNext.size());
  const unsigned long long lc[dev::L_COUNT] = {(unsigned long long)sim.n_restreamed,
                                               (unsigned long long)sim.n_trapped_total};
  par::upload(E.lcnt.data(), lc, dev::L_COUNT);
  par::fence();
}

// ----------------------------------------------------------------------------
//  one step s^0 -> s^1 (Simulation::step)
// ----------------------------------------------------------------------------
void Engine::step() {
  using namespace dev;
  EngineState& E = *st;
  // step of this s-level and its matrices (Simulation::beginStep)
  E.P.Delta = sim.Delta;
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
        const double q = 0.25 * Delta * P.zv.wLo(n + 1);
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
      par::forN(long(M + 1) * Nz * (N + 1), GEM_LAMBDA(long idx) {
        const int i = int(idx % (N + 1));
        const long mn = idx / (N + 1);
        const int m = int(mn / Nz), nn = int(mn % Nz);
        auto dmp = [&](double r) {
          if (r <= ra) return 1.0;
          const double x = (r - ra) / (rm - ra);
          return fmax(0.0, 1.0 - str * x * x);
        };
        fn.Rr[fn.f(m, nn, i)] = fn.Rr[fn.f(m, nn, i)] * dmp(g.rho[i]);
        if (i < N) fn.Rp[fn.c(m, nn, i)] = fn.Rp[fn.c(m, nn, i)] * dmp(g.c[i]);
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

// ----------------------------------------------------------------------------
void Engine::downloadFull() {
  download();
  EngineState& E = *st;
  auto dn = [](std::vector<cd>& h, par::Array<cplx>& d) {
    par::download(reinterpret_cast<cplx*>(h.data()), d.data(), h.size());
  };
  dn(sim.c0b.jr, E.Cb[E.iold].jr); dn(sim.c0b.jp, E.Cb[E.iold].jp); dn(sim.c0b.jz, E.Cb[E.iold].jz);
  dn(sim.c0s.jr, E.Cs[E.iold].jr); dn(sim.c0s.jp, E.Cs[E.iold].jp); dn(sim.c0s.jz, E.Cs[E.iold].jz);
  sim.pendingNext.resize(size_t(E.nPprev));
  par::download(sim.pendingNext.data(), E.pprev.data(), sim.pendingNext.size());
}

// ----------------------------------------------------------------------------
void Engine::download() {
  EngineState& E = *st;
  E.prof.start();
  par::fence();
  dev::FieldStore& Fo = E.F[E.iold];
  auto dn = [](std::vector<cd>& h, par::Array<cplx>& d) {
    par::download(reinterpret_cast<cplx*>(h.data()), d.data(), h.size());
  };
  dn(sim.fold.Rr, Fo.Rr); dn(sim.fold.Lr, Fo.Lr); dn(sim.fold.Bz, Fo.Bz);
  dn(sim.fold.Rp, Fo.Rp); dn(sim.fold.Lp, Fo.Lp); dn(sim.fold.Ez, Fo.Ez);
  par::download(sim.c0b.rho.data(), E.Cb[E.iold].rho.data(), E.nd);
  par::download(sim.c0b.rhoi.data(), E.Cb[E.iold].rhoi.data(), E.nd);
  par::download(sim.c0s.rho.data(), E.Cs[E.iold].rho.data(), E.nd);
  par::download(sim.c0s.rhoi.data(), E.Cs[E.iold].rhoi.data(), E.nd);
  sim.beam.resize(size_t(E.nBeam));
  par::download(sim.beam.data(), E.beam.data(), sim.beam.size());
  std::vector<double> ex(E.exitSt.capacity());
  par::download(ex.data(), E.exitSt.data(), ex.size());
  for (size_t s = 0; s < sim.exitStates.size(); ++s)
    for (size_t z = 0; z < sim.exitStates[s].size(); ++z) sim.exitStates[s][z] = ex[s * E.zStride + z];
  if (!sim.freedCharge.empty()) par::download(sim.freedCharge.data(), E.freed.data(), sim.freedCharge.size());
  unsigned long long lc[dev::L_COUNT];
  par::download(lc, E.lcnt.data(), dev::L_COUNT);
  sim.n_restreamed = long(lc[dev::L_RESTREAM]);
  sim.n_trapped_total = long(lc[dev::L_TRAPPED]);
  E.prof.stop(Profile::DOWNLOAD);
}

}  // namespace gem
