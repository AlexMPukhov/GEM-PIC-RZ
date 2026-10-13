// ============================================================================
//  fields.hpp -- storage of mode fields and currents
//
//  Fields of one s-level: slices n = 0..Nz (slice Nz is an empty ghost ahead
//  of the window), modes m = 0..M.  F(r,phi) = F_0 + 2 Re sum_{m>=1} F_m e^{i m phi}
//  Currents live on half-cells k = 0..Nz-1 (between slices k and k+1).
// ============================================================================
#pragma once
#include <algorithm>
#include <complex>
#include <vector>

namespace gem {

using cd = std::complex<double>;

struct Fields {
  int M = 0, Nz = 0, N = 0;
  std::vector<cd> Rr, Lr, Bz;  // faces   (N+1 per slice)
  std::vector<cd> Rp, Lp, Ez;  // centres (N   per slice)

  void init(int M_, int Nz_, int N_) {
    M = M_; Nz = Nz_; N = N_;
    size_t nf = size_t(M + 1) * (Nz + 1) * (N + 1), nc = size_t(M + 1) * (Nz + 1) * N;
    Rr.assign(nf, 0.0); Lr.assign(nf, 0.0); Bz.assign(nf, 0.0);
    Rp.assign(nc, 0.0); Lp.assign(nc, 0.0); Ez.assign(nc, 0.0);
  }
  size_t f(int m, int n, int i) const { return (size_t(m) * (Nz + 1) + n) * (N + 1) + i; }
  size_t c(int m, int n, int i) const { return (size_t(m) * (Nz + 1) + n) * N + i; }

  void copySlice(const Fields& o, int n) {
    for (int m = 0; m <= M; ++m) {
      std::copy_n(&o.Rr[f(m, n, 0)], N + 1, &Rr[f(m, n, 0)]);
      std::copy_n(&o.Lr[f(m, n, 0)], N + 1, &Lr[f(m, n, 0)]);
      std::copy_n(&o.Bz[f(m, n, 0)], N + 1, &Bz[f(m, n, 0)]);
      std::copy_n(&o.Rp[c(m, n, 0)], N, &Rp[c(m, n, 0)]);
      std::copy_n(&o.Lp[c(m, n, 0)], N, &Lp[c(m, n, 0)]);
      std::copy_n(&o.Ez[c(m, n, 0)], N, &Ez[c(m, n, 0)]);
    }
  }
};

struct Currents {
  int M = 0, Nz = 0, N = 0;
  std::vector<cd> jr;       // faces
  std::vector<cd> jp, jz;   // centres
  std::vector<double> rho;  // m = 0 charge density of electrons and bunches (diagnostic)
  std::vector<double> rhoi; // m = 0 charge density of mobile ions (diagnostic)

  void init(int M_, int Nz_, int N_) {
    M = M_; Nz = Nz_; N = N_;
    jr.assign(size_t(M + 1) * Nz * (N + 1), 0.0);
    jp.assign(size_t(M + 1) * Nz * N, 0.0);
    jz.assign(size_t(M + 1) * Nz * N, 0.0);
    rho.assign(size_t(Nz) * N, 0.0);
    rhoi.assign(size_t(Nz) * N, 0.0);
  }
  size_t f(int m, int k, int i) const { return (size_t(m) * Nz + k) * (N + 1) + i; }
  size_t c(int m, int k, int i) const { return (size_t(m) * Nz + k) * N + i; }

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
    std::fill_n(&rho[size_t(k) * N], N, 0.0);
    std::fill_n(&rhoi[size_t(k) * N], N, 0.0);
  }
};

}  // namespace gem
