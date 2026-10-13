// ============================================================================
//  particles.hpp -- macroparticle record and relativistic Boris pusher
// ============================================================================
#pragma once
#include <math.h>

#include "hd.hpp"

namespace gem {

// species id:
//   >= 0                             index of the bunch in the deck
//   SPECIES_TRAPPED            (-1)  plasma electrons trapped and moved to the beam list
//   SPECIES_PLASMA             (-2)  streaming plasma electrons
//   SPECIES_TRAPPED_ION        (-3)  ions that reached the next s-level (not expected)
//   SPECIES_ION_BASE - s       (-10 ... -19)  streaming ions of ion block s
//   SPECIES_IONIZED_BASE - s   (-20 ... -29)  electrons freed by ionization of ion block s
//   SPECIES_TRAPPED_IONIZED_BASE - s  (-30 ... -39)  such electrons after trapping
// particle status: 0 alive, 1 trapped (reached the next s-level in this slice),
// 2 lost, 3 beam particle returning to streaming, STATUS_WAITING streaming
// particle that joins the march when it reaches its zeta
constexpr int STATUS_WAITING = 4;

constexpr int SPECIES_TRAPPED = -1;
constexpr int SPECIES_PLASMA = -2;
constexpr int SPECIES_TRAPPED_ION = -3;
constexpr int SPECIES_ION_BASE = -10;
constexpr int SPECIES_IONIZED_BASE = -20;
constexpr int SPECIES_TRAPPED_IONIZED_BASE = -30;
GEM_HD inline bool isIon(int sp) { return (sp <= SPECIES_ION_BASE && sp > SPECIES_IONIZED_BASE) || sp == SPECIES_TRAPPED_ION; }
GEM_HD inline bool isIonizationElectron(int sp) { return sp <= SPECIES_IONIZED_BASE && sp > SPECIES_TRAPPED_IONIZED_BASE; }
GEM_HD inline int ionIndex(int sp) { return SPECIES_ION_BASE - sp; }            // for streaming ions
GEM_HD inline int ionizedIndex(int sp) { return SPECIES_IONIZED_BASE - sp; }    // for ionization electrons

struct Particle {
  double x, y, z;     // Cartesian transverse position, z = lab coordinate (= s for GT2)
  double ux, uy, uz;  // gamma v
  double q;           // macroparticle charge [code units], signed
  double qm;          // charge-to-mass ratio in units of e/m_e (-1 electron, +1 positron)
  double zeta;        // fast coordinate
  int status;         // see above
  int sp;             // species id (see above)
  double w;           // ions: number of atoms (in units of the code charge), q = Zi * w
  int Zi;             // ions: charge state
  int sub;            // MPI pipeline: next sub-step of a beam push continued on the next rank
  unsigned long long id;  // parallel engine: key of the random numbers (ionization)
};

struct EB {
  double Ex, Ey, Ez, Bx, By, Bz;
};

// relativistic Boris push, charge-to-mass ratio qm, time step dt
GEM_HD inline void boris(double& ux, double& uy, double& uz, const EB& f, double qm, double dt) {
  const double qh = 0.5 * qm * dt;
  double umx = ux + qh * f.Ex, umy = uy + qh * f.Ey, umz = uz + qh * f.Ez;
  double gm = sqrt(1.0 + umx * umx + umy * umy + umz * umz);
  double tx = qh * f.Bx / gm, ty = qh * f.By / gm, tz = qh * f.Bz / gm;
  double t2 = tx * tx + ty * ty + tz * tz;
  double upx = umx + (umy * tz - umz * ty);
  double upy = umy + (umz * tx - umx * tz);
  double upz = umz + (umx * ty - umy * tx);
  double s = 2.0 / (1.0 + t2);
  double vx = umx + s * (upy * tz - upz * ty);
  double vy = umy + s * (upz * tx - upx * tz);
  double vz = umz + s * (upx * ty - upy * tx);
  ux = vx + qh * f.Ex; uy = vy + qh * f.Ey; uz = vz + qh * f.Ez;
}

GEM_HD inline double gammaOf(const Particle& p) { return sqrt(1.0 + p.ux * p.ux + p.uy * p.uy + p.uz * p.uz); }

}  // namespace gem
