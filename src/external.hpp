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

// transport elements (deck blocks quadrupole, plasma_lens, solenoid), on [z0, z1]:
//   quadrupole   B_x = -g y, B_y = -g x       (g > 0: focuses electrons in x, defocuses in y)
//   plasma lens  B_phi = -g r for r < R, -g R^2 / r beyond  (g > 0: focuses electrons)
//   solenoid     B_z = b f(z), B_r = -(r/2) b f'(z)  (f: linear fringes of length `ramp`)
struct MagnetV {
  int kind;       // 0 quadrupole, 1 plasma lens, 2 solenoid
  double z0, z1;
  double g;       // gradient (quadrupole, lens) or field (solenoid), code units
  double R;       // plasma lens radius
  double ramp;    // solenoid fringe
};

struct ExternalV {
  int n = 0;
  const UndulatorV* u = nullptr;
  const double* k = nullptr;  // taper tables
  int nm = 0;
  const MagnetV* m = nullptr;
  GEM_HD bool any() const { return n > 0 || nm > 0; }
};

GEM_HD inline double undulatorEnvelope(const UndulatorV& U, double z) {
  if (z <= U.z0 || z >= U.z1) return 0.0;
  if (U.ramp <= 0.0) return 1.0;
  return fmin(1.0, fmin((z - U.z0) / U.ramp, (U.z1 - z) / U.ramp));
}

GEM_HD inline void addMagnets(const ExternalV& X, double x, double y, double z, EB& f) {
  for (int k = 0; k < X.nm; ++k) {
    const MagnetV& G = X.m[k];
    if (z < G.z0 || z >= G.z1) continue;
    if (G.kind == 0) {
      f.Bx -= G.g * y;
      f.By -= G.g * x;
    } else if (G.kind == 1) {
      const double r2 = x * x + y * y, R2 = G.R * G.R;
      const double s = (r2 <= R2) ? G.g : G.g * R2 / r2;  // B_phi / r
      f.Bx += s * y;   // B_x = -B_phi sin(phi) = g y
      f.By -= s * x;   // B_y =  B_phi cos(phi) = -g x
    } else {
      double fz = 1.0, dfz = 0.0;
      if (G.ramp > 0.0) {
        if (z - G.z0 < G.ramp) { fz = (z - G.z0) / G.ramp; dfz = 1.0 / G.ramp; }
        else if (G.z1 - z < G.ramp) { fz = (G.z1 - z) / G.ramp; dfz = -1.0 / G.ramp; }
      }
      f.Bz += G.g * fz;
      f.Bx -= 0.5 * x * G.g * dfz;
      f.By -= 0.5 * y * G.g * dfz;
    }
  }
}

// add the external fields at the lab position (x, y, z) to f
GEM_HD inline void addExternal(const ExternalV& X, double x, double y, double z, EB& f) {
  if (X.nm) addMagnets(X, x, y, z, f);
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
