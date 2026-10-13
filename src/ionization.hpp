// ============================================================================
//  ionization.hpp -- ADK tunnel-ionization rates for the levels of one species
//
//  Rate for the transition k -> k+1 (Ammosov, Delone, Krainov, JETP 64, 1191
//  (1986); in the form used e.g. by Chen et al., JCP 236, 220 (2013) and FBPIC,
//  for s-states with l* = n* - 1, m = 0):
//
//    w = omega_a C^2 (U/2U_H) [2 (U/U_H)^{3/2} E_a/E]^{2n*-1} exp[-(2/3)(U/U_H)^{3/2} E_a/E]
//    n* = (k+1) sqrt(U_H/U),   C^2 = 2^{2n*} / (n* Gamma(n*+1) Gamma(n*))
//
//  U = ionization energy of the level (NIST), E = instantaneous |E| at the
//  atom (the laser carrier is resolved, so the static-field rate applies),
//  U_H = 13.6057 eV, E_a = 5.1422e11 V/m, omega_a = 4.1341e16 1/s.
//  For hydrogen the formula reduces to the Landau rate 4 omega_a (E_a/E)
//  exp(-2 E_a / 3E).
//  Beyond E_peak = (2/3)(U/U_H)^{3/2} E_a / (2n*-1) the ADK expression
//  decreases again (barrier-suppression regime); the rate is kept at its
//  maximum there.
// ============================================================================
#pragma once
#include <vector>

#include "units.hpp"

namespace gem {

constexpr double U_H_EV = 13.605693122994;           // Rydberg energy [eV]
constexpr double E_ATOMIC_VM = 5.14220674763e11;     // atomic field [V/m]
constexpr double OMEGA_ATOMIC = 4.1341373335e16;     // atomic frequency [1/s]

// reference implementation in SI units: rate [1/s] for a level with energy U [eV]
// of an ion that has charge Zafter after the ionization, in the field E [V/m]
double adkRateSI(double U_eV, int Zafter, double E_Vm);

class AdkRates {
 public:
  // levels k = Z0 ... Zmax-1 with energies ip_eV[k]; code units from `u`;
  // T = longest time a particle spends in the window (for the skip threshold)
  void setup(const std::vector<double>& ip_eV, int Z0, int Zmax, const Units& u, double T);

  int zmax() const { return Zmax; }
  // rate of k -> k+1 [1/omega0] in the field E [m c omega0 / e]; 0 below the threshold
  double rate(int k, double E) const;
  // below this field no level can be ionized within the window time (probability < 1e-10)
  double minField() const { return Emin_all; }
  // field [m c omega0 / e] at which the rate of level k reaches 1/tau
  double fieldForRate(int k, double rate_code) const;

 private:
  friend class Engine;  // uploads the tables
  struct Level {
    double logpre = 0, power = 0, expo = 0;  // w = exp(logpre - power ln E + expo / E)
    double Epeak = 1e300;                    // maximum of the ADK expression
    double Emin = 0;                         // w(Emin) * T = 1e-10
  };
  int Z0 = 0, Zmax = 0;
  std::vector<Level> lev;  // index k
  double Emin_all = 0;
  double raw(const Level& L, double E) const;
};

}  // namespace gem
