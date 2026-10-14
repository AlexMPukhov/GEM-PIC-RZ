// ============================================================================
//  external.hpp -- external (applied) fields acting on the particles
//
//  Undulators (deck blocks 'undulator'), static magnetic fields in the lab
//  frame, functions of the lab position (x, y, z) of a particle. They act on
//  the particles only (added to the gathered self-fields before the push);
//  the field solver never sees them. With k_u = 2 pi / period, B0 = K k_u
//  (code units: B in m omega0 / e, k in k0), phase psi = k_u (z - z_start)
//  and the envelope g(z) (linear ramps of ramp length at both ends, 0 outside):
//
//    planar:   B_y = B0 g cosh(k_u y) cos(psi),  B_z = -B0 g sinh(k_u y) sin(psi)
//              (curl-free, div-free for g = 1: natural focusing in y); without
//              natural focusing only B_y = B0 g cos(psi)
//    helical:  B_x = B0 g cos(psi),  B_y = B0 g sin(psi)   (on-axis field)
//
//  An electron (q/m = -1) entering with u_x = 0 follows
//  u_x = K g sin(psi) (planar): no drift if the ramp is slow.
//  A taper K(z) (deck: K_function) replaces B0 by K(z) k_u from a table.
// ============================================================================
#pragma once
#include <math.h>

#include "hd.hpp"
#include "particles.hpp"

namespace gem {

struct UndulatorV {
  int helical;    // 0 planar, 1 helical
  int focusing;   // planar: natural focusing (cosh / sinh terms)
  double ku, B0;  // wavenumber, peak field (code units)
  double z0, z1;  // start and end (lab z)
  double ramp;    // length of the linear ramps at both ends
  int nk, koff;   // taper: nk > 0 values B0(z0 + j dzk) at ExternalV::k[koff + j]
  double dzk;
};

struct ExternalV {
  int n = 0;
  const UndulatorV* u = nullptr;
  const double* k = nullptr;  // taper tables
};

GEM_HD inline double undulatorEnvelope(const UndulatorV& U, double z) {
  if (z <= U.z0 || z >= U.z1) return 0.0;
  if (U.ramp <= 0.0) return 1.0;
  return fmin(1.0, fmin((z - U.z0) / U.ramp, (U.z1 - z) / U.ramp));
}

// add the external fields at the lab position (x, y, z) to f
GEM_HD inline void addExternal(const ExternalV& X, double x, double y, double z, EB& f) {
  (void)x;
  for (int k = 0; k < X.n; ++k) {
    const UndulatorV& U = X.u[k];
    const double g = undulatorEnvelope(U, z);
    if (g == 0.0) continue;
    const double psi = U.ku * (z - U.z0);
    double b = U.B0 * g;
    if (U.nk > 0) {  // K(z): linear interpolation in the table
      const double t = (z - U.z0) / U.dzk;
      const int j = t < double(U.nk - 2) ? int(t) : U.nk - 2;
      const double w = t - j;
      b = g * ((1.0 - w) * X.k[U.koff + j] + w * X.k[U.koff + j + 1]);
    }
    if (U.helical) {
      f.Bx += b * cos(psi);
      f.By += b * sin(psi);
    } else if (U.focusing) {
      const double ky = U.ku * y;
      f.By += b * cosh(ky) * cos(psi);
      f.Bz -= b * sinh(ky) * sin(psi);
    } else {
      f.By += b * cos(psi);
    }
  }
}

}  // namespace gem
