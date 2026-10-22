// ============================================================================
//  engine_kernels.hpp -- kernels of the parallel engine (see engine.hpp); used by src/engine*.cpp
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
#pragma once
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
  double theta;      // implicitness of the R update (numerics.theta)
  double trapExtra;  // (trap_advance - 1) * Delta: advance in z beyond s^1 before a particle leaves the level
  ZetaView zv;   // longitudinal grid
  const int* cls;  // matrix class of slice n (FieldSolver::classOf)
  ExternalV ext;   // external fields (undulators)
  int vay;         // beam list: Vay pusher (else Boris)
  double beamDzMax;  // numerics.beam_dz_max (sub-steps of the beam push)
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
  const int nsub = P.beamDzMax > 0.0 ? (int)ceil(P.Delta / P.beamDzMax) : 1;
  double ux, uy, uz;
  velocityAtPosition(P.ext, p, P.Delta / (nsub > 1 ? nsub : 1), ux, uy, uz);
  const double gm = sqrt(1.0 + ux * ux + uy * uy + uz * uz);
  const double vx = ux / gm, vy = uy / gm, vz = uz / gm;
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
  if (P.ext.any()) addExternal(P.ext, xm, ym, p.z + 0.5 * dt0 * p.uz / g0, f, dt0 * p.uz / g0);
  double ux = p.ux, uy = p.uy, uz = p.uz;
  boris(ux, uy, uz, f, p.qm, dt0);
  double g1 = sqrt(1.0 + ux * ux + uy * uy + uz * uz);
  double dt = 2.0 * hl / (omv0 + (g1 - uz) / g1);
  ux = p.ux; uy = p.uy; uz = p.uz;
  boris(ux, uy, uz, f, p.qm, dt);
  g1 = sqrt(1.0 + ux * ux + uy * uy + uz * uz);
  double vx = 0.5 * (p.ux / g0 + ux / g1), vy = 0.5 * (p.uy / g0 + uy / g1), vz = 0.5 * (p.uz / g0 + uz / g1);
  double zeta_end = zn;
  const double sTrap = s1 + P.trapExtra;
  if (p.z + vz * dt >= sTrap) {
    const double dtp = fmax(0.0, (sTrap - p.z) / fmax(vz, 1e-3));
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
    addExternal(P.ext, p.x, p.y, p.z, f, dz);
    if (P.vay) vay(p.ux, p.uy, p.uz, f, p.qm, dt);
    else boris(p.ux, p.uy, p.uz, f, p.qm, dt);
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

// plasma.smooth_length (SourceFilter, smoothing.cpp): factorized filters per mode
struct FiltV {
  int on, N;
  const cplx* band;   // per mode 5 (2N+1): LU of 1 - a^2 Lap_vec, offsets -2..+2 per row (1/diagonal)
  const double *sa, *scp, *sden;  // per mode N: Thomas factors of 1 - a^2 Lap_m
};

GEM_HD inline void filterScalarDev(const FiltV& F, int m, cplx* f, cplx* w) {
  const int N = F.N;
  const double *a = F.sa + size_t(m) * N, *cp = F.scp + size_t(m) * N, *den = F.sden + size_t(m) * N;
  w[0] = f[0] * den[0];
  for (int i = 1; i < N; ++i) w[i] = (f[i] - a[i] * w[i - 1]) * den[i];
  f[N - 1] = w[N - 1];
  for (int i = N - 2; i >= 0; --i) f[i] = w[i] - cp[i] * f[i + 1];
}

// smooth the streaming currents of half-cell n, mode m (leader of one team; Simulation::smoothStream)
GEM_HD inline void filterTeam(const Team& t, const FiltV& F, int m, int n, const CurrV& c, cplx* x) {
  if (!t.leader()) return;
  const int N = F.N, nb = 2 * N + 1;
  cplx* w = x + nb;
  cplx* y = w + N;
  for (int i = 0; i <= N; ++i) x[2 * i] = c.jr[c.f(m, n, i)];
  for (int i = 0; i < N; ++i) x[2 * i + 1] = c.jp[c.c(m, n, i)];
  const cplx* A = F.band + size_t(m) * 5 * nb;
  for (int i = 0; i < nb; ++i)
    for (int p = -2; p <= -1; ++p)
      if (i + p >= 0) x[i] -= A[size_t(i) * 5 + size_t(p + 2)] * x[i + p];
  for (int i = nb - 1; i >= 0; --i) {
    for (int k = 1; k <= 2; ++k)
      if (i + k < nb) x[i] -= A[size_t(i) * 5 + size_t(k + 2)] * x[i + k];
    x[i] = x[i] * A[size_t(i) * 5 + 2];  // reciprocal of the diagonal (Banded5::factor)
  }
  for (int i = 1; i <= N; ++i) c.jr[c.f(m, n, i)] = x[2 * i];
  for (int i = 0; i < N; ++i) c.jp[c.c(m, n, i)] = x[2 * i + 1];
  for (int i = 0; i < N; ++i) y[i] = c.jz[c.c(m, n, i)];
  filterScalarDev(F, m, y, w);
  for (int i = 0; i < N; ++i) c.jz[c.c(m, n, i)] = y[i];
  if (m == 0)
    for (int k = 0; k < 2; ++k) {
      double* d = (k == 0 ? c.rho : c.rhoi) + size_t(n) * N;
      for (int i = 0; i < N; ++i) y[i] = cplx(d[i]);
      filterScalarDev(F, 0, y, w);
      for (int i = 0; i < N; ++i) d[i] = y[i].re;
    }
}

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
  const double th1 = 2.0 * P.theta, th0 = 2.0 - th1;  // 1, 1 for the trapezoid
  auto Jr = [&](int i) {
    const cplx a = c0b.jr[c0b.f(m, nm, i)] + c0s.jr[c0s.f(m, nm, i)];
    const cplx pr = (n > 0) ? predR[m * (N + 1) + i] : c1b.jr[c1b.f(m, 0, i)] + c1s.jr[c1s.f(m, 0, i)];
    return th0 * (wl * a + wh * c0b.jr[c0b.f(m, n, i)] + wh * c0s.jr[c0s.f(m, n, i)]) +
           th1 * (wl * pr + wh * c1b.jr[c1b.f(m, n, i)] + wh * c1s.jr[c1s.f(m, n, i)]);
  };
  auto Jp = [&](int i) {
    const cplx a = c0b.jp[c0b.c(m, nm, i)] + c0s.jp[c0s.c(m, nm, i)];
    const cplx pr = (n > 0) ? predP[m * N + i] : c1b.jp[c1b.c(m, 0, i)] + c1s.jp[c1s.c(m, 0, i)];
    return th0 * (wl * a + wh * c0b.jp[c0b.c(m, n, i)] + wh * c0s.jp[c0s.c(m, n, i)]) +
           th1 * (wl * pr + wh * c1b.jp[c1b.c(m, n, i)] + wh * c1s.jp[c1s.c(m, n, i)]);
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
    Pr[i] = Rr0[i] + (1.0 - P.theta) * D * (opG(e0, i) + opMf(b0, i)) - 0.25 * D * Jr(i);
    Qr[i] = Lr1[i] - 0.25 * h * (opG(e1, i) - opMf(b1, i)) - 0.5 * h * jr1(i);
  });
  par::teamFor(t, N, [&](int i) {
    Pp[i] = Rp0[i] + (1.0 - P.theta) * D * ((im / g.c[i]) * e0[i] - opg(b0, i)) - 0.25 * D * Jp(i);
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
    Rr[i] = Pr[i] + P.theta * D * (G + Mf);
    Lr[i] = Qr[i] - 0.25 * h * (G - Mf);
  });
  par::teamFor(t, N, [&](int i) {
    const cplx gg = opg(b, i), ie = (im / g.c[i]) * e[i];
    Rp[i] = Pp[i] + P.theta * D * (ie - gg);
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
}  // namespace gem
