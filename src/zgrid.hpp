// ============================================================================
//  zgrid.hpp -- longitudinal grid of the window (uniform or free-hand)
//
//    slices      zeta_n,  n = 0..Nz           (zeta_0 = 0 rear end; slice Nz
//                                              is the empty ghost ahead)
//    half-cells  [zeta_k, zeta_k+1], k = 0..Nz-1, length h_k, centre zc_k
//
//  The field equations of the scheme that couple two slices (L, B_z, E_z)
//  are two-point (trapezoidal) rules over one half-cell and hold for any h_k.
//  The equation of R lives on the slice: its current is the linear
//  interpolation of the two neighbouring half-cells to zeta_n,
//      j(zeta_n) = (h_n j_{n-1/2} + h_{n-1} j_{n+1/2}) / (h_{n-1} + h_n)
//  (the plain average on a uniform grid); the half-cell behind the window is a
//  copy of half-cell 0 (h_{-1} = h_0).
//
//  ZetaView is a plain view of the arrays, used by the serial code and by the
//  kernels of the parallel engine. On a uniform grid all functions use the
//  arithmetic of the uniform code (n * h, zeta / h), so a uniform run is
//  unchanged bit for bit.
// ============================================================================
#pragma once
#include <math.h>

#include <vector>

#include "hd.hpp"

namespace gem {

struct ZetaView {
  int Nz = 0;
  int uniform = 1;
  double h = 0.0;                // uniform step
  const double* z = nullptr;     // zeta_0 .. zeta_Nz
  const double* hs = nullptr;    // h_0 .. h_{Nz-1}
  const double* zc = nullptr;    // zc[j] = centre of half-cell j-1, j = 0..Nz+1 (virtual -1 and Nz)
  const double* wlo = nullptr;   // slice n: 2 h_n / (h_{n-1} + h_n)   (weight of half-cell n-1, x 2)
  const double* whi = nullptr;   // slice n: 2 h_{n-1} / (h_{n-1} + h_n) (weight of half-cell n, x 2)

  // number of entries of x[0..n) that are <= v (std::upper_bound)
  GEM_HD static int upper(const double* x, int n, double v) {
    int lo = 0, hi = n;
    while (lo < hi) {
      const int mid = (lo + hi) >> 1;
      if (x[mid] > v) hi = mid; else lo = mid + 1;
    }
    return lo;
  }
  // number of entries of x[0..n) that are < v (std::lower_bound)
  GEM_HD static int lower(const double* x, int n, double v) {
    int lo = 0, hi = n;
    while (lo < hi) {
      const int mid = (lo + hi) >> 1;
      if (x[mid] >= v) hi = mid; else lo = mid + 1;
    }
    return lo;
  }

  GEM_HD double node(int n) const { return uniform ? n * h : z[n]; }
  GEM_HD double len(int k) const { return uniform ? h : hs[k]; }
  GEM_HD double wLo(int n) const { return uniform ? 1.0 : wlo[n]; }
  GEM_HD double wHi(int n) const { return uniform ? 1.0 : whi[n]; }
  // centre of half-cell k, k = -1 .. Nz
  GEM_HD double centre(int k) const { return uniform ? (k + 0.5) * h : zc[k + 1]; }

  // slice interval of zeta: n with zeta_n <= zeta < zeta_n+1 and the linear
  // weight w of slice n+1 (n < 0 or n >= Nz: outside the window)
  GEM_HD int cell(double zeta, double& w) const {
    if (uniform) {
      const double zf = zeta / h;
      const int na = int(floor(zf));
      w = zf - na;
      return na;
    }
    w = 0.0;
    if (zeta < z[0]) return -1;
    if (zeta >= z[Nz]) return Nz;
    const int na = upper(z, Nz + 1, zeta) - 1;
    w = (zeta - z[na]) / hs[na];
    return na;
  }

  // beam deposit, linear between half-cell centres: half-cells j-1 (weight
  // 1-f) and j (weight f), where j is the slice nearest to zeta
  GEM_HD int beamCell(double zeta, double& f) const {
    if (uniform) {
      const double kf = zeta / h - 0.5;
      const int k = int(floor(kf));
      f = kf - k;
      return k + 1;
    }
    f = 0.0;
    const int j = upper(zc, Nz + 2, zeta) - 1;  // zc[j] <= zeta < zc[j+1]
    if (j < 0 || j > Nz) return j;
    f = (zeta - zc[j]) / (zc[j + 1] - zc[j]);
    return j;
  }

  // slice nearest to zeta (= beamCell; the owner of a beam particle under MPI)
  GEM_HD int nearest(double zeta) const {
    if (uniform) return int(floor(zeta / h + 0.5));
    return upper(zc, Nz + 2, zeta) - 1;
  }

  // half-cell k with zeta_k < zeta <= zeta_k+1 (where a re-streamed particle joins the march)
  GEM_HD int joinCell(double zeta) const {
    if (uniform) return int(ceil(zeta / h)) - 1;
    return lower(z, Nz + 1, zeta) - 1;
  }
};

// owner of the arrays (host)
struct ZetaGrid {
  int Nz = 0;
  bool uniform = true;
  double h = 0.0;
  std::vector<double> z, hs, zc, wlo, whi;

  // nodes zeta_0..zeta_Nz (uniform: n * h exactly)
  void build(const std::vector<double>& nodes, bool isUniform, double hUniform) {
    Nz = int(nodes.size()) - 1;
    uniform = isUniform;
    h = hUniform;
    z = nodes;
    hs.resize(Nz);
    for (int k = 0; k < Nz; ++k) hs[k] = uniform ? h : z[k + 1] - z[k];
    zc.resize(Nz + 2);
    zc[0] = z[0] - 0.5 * hs[0];
    for (int k = 0; k < Nz; ++k) zc[k + 1] = uniform ? (k + 0.5) * h : 0.5 * (z[k] + z[k + 1]);
    zc[Nz + 1] = z[Nz] + 0.5 * hs[Nz - 1];
    wlo.resize(Nz);
    whi.resize(Nz);
    for (int n = 0; n < Nz; ++n) {
      const double hm = (n > 0) ? hs[n - 1] : hs[0], hp = hs[n];
      wlo[n] = uniform ? 1.0 : 2.0 * hp / (hm + hp);
      whi[n] = uniform ? 1.0 : 2.0 * hm / (hm + hp);
    }
  }
  ZetaView view() const {
    return ZetaView{Nz, uniform ? 1 : 0, h, z.data(), hs.data(), zc.data(), wlo.data(), whi.data()};
  }
};

}  // namespace gem
