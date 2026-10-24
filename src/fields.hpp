// ============================================================================
//  fields.hpp -- storage of mode fields and currents
//
//  Fields of one s-level: slices n = 0..Nz (slice Nz is an empty ghost ahead
//  of the window), modes m = 0..M.  F(r,phi) = F_0 + 2 Re sum_{m>=1} F_m e^{i m phi}
//  Currents live on half-cells k = 0..Nz-1 (between slices k and k+1).
//
//  Storage may cover only a part of the window (MPI ranks: their block and
//  ghosts): Fields hold the slices [n0, n0 + ns), Currents the half-cells
//  [k0, k0 + nk). The accessors take the GLOBAL slice / half-cell index, so
//  the code is the same for full and partial storage. GEMPIC_RANGE_CHECK
//  makes every access outside the stored range throw (debug builds).
// ============================================================================
#pragma once
#include <algorithm>
#include <complex>
#include <stdexcept>
#include <string>
#include <vector>

namespace gem {

using cd = std::complex<double>;

#if defined(GEMPIC_RANGE_CHECK)
inline void rangeFail(const char* what, int n, int lo, int cnt) {
  throw std::out_of_range(std::string(what) + " " + std::to_string(n) + " outside the stored range [" +
                          std::to_string(lo) + ", " + std::to_string(lo + cnt) + ")");
}
#define GEM_RANGE(what, n, lo, cnt) \
  do { if ((n) < (lo) || (n) >= (lo) + (cnt)) rangeFail(what, n, lo, cnt); } while (0)
#else
#define GEM_RANGE(what, n, lo, cnt) do { } while (0)
#endif

struct Fields {
  int M = 0, Nz = 0, N = 0;
  int n0 = 0, ns = 0;          // stored slices [n0, n0 + ns)
  std::vector<cd> Rr, Lr, Bz;  // faces   (N+1 per slice)
  std::vector<cd> Rp, Lp, Ez;  // centres (N   per slice)

  // all slices 0..Nz (Nz: the ghost ahead of the window)
  void init(int M_, int Nz_, int N_) { initRange(M_, Nz_, N_, 0, Nz_ + 1); }
  // slices [lo, hi) only (clipped to 0..Nz)
  void initRange(int M_, int Nz_, int N_, int lo, int hi) {
    M = M_; Nz = Nz_; N = N_;
    n0 = std::max(0, lo);
    ns = std::min(Nz + 1, hi) - n0;
    size_t nf = size_t(M + 1) * ns * (N + 1), nc = size_t(M + 1) * ns * N;
    Rr.assign(nf, 0.0); Lr.assign(nf, 0.0); Bz.assign(nf, 0.0);
    Rp.assign(nc, 0.0); Lp.assign(nc, 0.0); Ez.assign(nc, 0.0);
  }
  bool full() const { return n0 == 0 && ns == Nz + 1; }
  bool has(int n) const { return n >= n0 && n < n0 + ns; }
  size_t f(int m, int n, int i) const {
    GEM_RANGE("slice", n, n0, ns);
    return (size_t(m) * ns + (n - n0)) * (N + 1) + i;
  }
  size_t c(int m, int n, int i) const {
    GEM_RANGE("slice", n, n0, ns);
    return (size_t(m) * ns + (n - n0)) * N + i;
  }
  // copy the slices [lo, hi) present in both
  void copyRange(const Fields& o, int lo, int hi) {
    lo = std::max({lo, n0, o.n0});
    hi = std::min({hi, n0 + ns, o.n0 + o.ns});
    for (int n = lo; n < hi; ++n) copySlice(o, n);
  }

  void copySlice(const Fields& o, int n) {
    for (int m = 0; m <= M; ++m) {
      std::copy_n(&o.Rr[o.f(m, n, 0)], N + 1, &Rr[f(m, n, 0)]);
      std::copy_n(&o.Lr[o.f(m, n, 0)], N + 1, &Lr[f(m, n, 0)]);
      std::copy_n(&o.Bz[o.f(m, n, 0)], N + 1, &Bz[f(m, n, 0)]);
      std::copy_n(&o.Rp[o.c(m, n, 0)], N, &Rp[c(m, n, 0)]);
      std::copy_n(&o.Lp[o.c(m, n, 0)], N, &Lp[c(m, n, 0)]);
      std::copy_n(&o.Ez[o.c(m, n, 0)], N, &Ez[c(m, n, 0)]);
    }
  }
};

struct Currents {
  int M = 0, Nz = 0, N = 0;
  int k0 = 0, nk = 0;       // stored half-cells [k0, k0 + nk)
  std::vector<cd> jr;       // faces
  std::vector<cd> jp, jz;   // centres
  std::vector<double> rho;  // m = 0 charge density of electrons and bunches (diagnostic)
  std::vector<double> rhoi; // m = 0 charge density of mobile ions (diagnostic)

  void init(int M_, int Nz_, int N_) { initRange(M_, Nz_, N_, 0, Nz_); }
  // half-cells [lo, hi) only (clipped to 0..Nz-1)
  void initRange(int M_, int Nz_, int N_, int lo, int hi) {
    M = M_; Nz = Nz_; N = N_;
    k0 = std::max(0, lo);
    nk = std::max(0, std::min(Nz, hi) - k0);
    jr.assign(size_t(M + 1) * nk * (N + 1), 0.0);
    jp.assign(size_t(M + 1) * nk * N, 0.0);
    jz.assign(size_t(M + 1) * nk * N, 0.0);
    rho.assign(size_t(nk) * N, 0.0);
    rhoi.assign(size_t(nk) * N, 0.0);
  }
  bool full() const { return k0 == 0 && nk == Nz; }
  bool has(int k) const { return k >= k0 && k < k0 + nk; }
  size_t f(int m, int k, int i) const {
    GEM_RANGE("half-cell", k, k0, nk);
    return (size_t(m) * nk + (k - k0)) * (N + 1) + i;
  }
  size_t c(int m, int k, int i) const {
    GEM_RANGE("half-cell", k, k0, nk);
    return (size_t(m) * nk + (k - k0)) * N + i;
  }
  // densities rho, rhoi of half-cell k at centre i
  size_t d(int k, int i) const {
    GEM_RANGE("half-cell", k, k0, nk);
    return size_t(k - k0) * N + i;
  }

  // copy the half-cells [lo, hi) present in both
  void copyRange(const Currents& o, int lo, int hi) {
    lo = std::max({lo, k0, o.k0});
    hi = std::min({hi, k0 + nk, o.k0 + o.nk});
    for (int k = lo; k < hi; ++k) {
      for (int m = 0; m <= M; ++m) {
        std::copy_n(&o.jr[o.f(m, k, 0)], N + 1, &jr[f(m, k, 0)]);
        std::copy_n(&o.jp[o.c(m, k, 0)], N, &jp[c(m, k, 0)]);
        std::copy_n(&o.jz[o.c(m, k, 0)], N, &jz[c(m, k, 0)]);
      }
      std::copy_n(&o.rho[o.d(k, 0)], N, &rho[d(k, 0)]);
      std::copy_n(&o.rhoi[o.d(k, 0)], N, &rhoi[d(k, 0)]);
    }
  }
  void zero() {
    std::fill(jr.begin(), jr.end(), 0.0); std::fill(jp.begin(), jp.end(), 0.0);
    std::fill(jz.begin(), jz.end(), 0.0); std::fill(rho.begin(), rho.end(), 0.0);
    std::fill(rhoi.begin(), rhoi.end(), 0.0);
  }
  void zeroHalf(int k) {
    for (int m = 0; m <= M; ++m) {
      std::fill_n(&jr[f(m, k, 0)], N + 1, 0.0);
      std::fill_n(&jp[c(m, k, 0)], N, 0.0);
      std::fill_n(&jz[c(m, k, 0)], N, 0.0);
    }
    std::fill_n(&rho[d(k, 0)], N, 0.0);
    std::fill_n(&rhoi[d(k, 0)], N, 0.0);
  }
};

}  // namespace gem
