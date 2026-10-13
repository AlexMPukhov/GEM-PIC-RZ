// ============================================================================
//  ionization.cpp -- ADK rates (see ionization.hpp)
// ============================================================================
#include "ionization.hpp"

#include <algorithm>
#include <cmath>

namespace gem {

namespace {
// ln of omega_a C^2 (U/2U_H) (2 (U/U_H)^{3/2} E_a)^{2n*-1}  (SI, E in V/m),
// and the power 2n*-1 and exponent factor -(2/3)(U/U_H)^{3/2} E_a
void adkCoefficients(double U_eV, int Zafter, double& logpre, double& power, double& expo) {
  const double x = U_eV / U_H_EV;
  const double nstar = Zafter / std::sqrt(x);
  const double logC2 = 2.0 * nstar * std::log(2.0) - std::log(nstar) - std::lgamma(nstar + 1.0) - std::lgamma(nstar);
  power = 2.0 * nstar - 1.0;
  logpre = std::log(OMEGA_ATOMIC) + logC2 + std::log(0.5 * x) + power * std::log(2.0 * std::pow(x, 1.5) * E_ATOMIC_VM);
  expo = -(2.0 / 3.0) * std::pow(x, 1.5) * E_ATOMIC_VM;
}
}  // namespace

double adkRateSI(double U_eV, int Zafter, double E_Vm) {
  if (!(E_Vm > 0)) return 0.0;
  double logpre, power, expo;
  adkCoefficients(U_eV, Zafter, logpre, power, expo);
  return std::exp(logpre - power * std::log(E_Vm) + expo / E_Vm);
}

void AdkRates::setup(const std::vector<double>& ip_eV, int Z0_, int Zmax_, const Units& u, double T) {
  Z0 = Z0_;
  Zmax = Zmax_;
  lev.assign(std::max(Zmax, 0), Level());
  // code units: E_SI = E * E0,  w_code = w_SI / omega0
  const double E0 = 2.0 * PI * MC2_MEV * 1e6 / (u.lambda0_um * 1e-6);    // V/m
  const double omega0 = 2.0 * PI * C_LIGHT / (u.lambda0_um * 1e-6);      // 1/s
  Emin_all = 1e300;
  for (int k = Z0; k < Zmax; ++k) {
    Level& L = lev[k];
    double logpre, power, expo;
    adkCoefficients(ip_eV[k], k + 1, logpre, power, expo);
    // w_code(E) = w_SI(E E0) / omega0
    L.power = power;
    L.logpre = logpre - std::log(omega0) - power * std::log(E0);
    L.expo = expo / E0;
    L.Epeak = (power > 0) ? -L.expo / power : 1e300;
    // skip threshold: w(Emin) T = 1e-10 (w is increasing below Epeak)
    double lo = 1e-12, hi = std::min(L.Epeak, 1e6);
    if (raw(L, hi) * T < 1e-10) {
      L.Emin = hi;
    } else {
      for (int it = 0; it < 200; ++it) {
        const double mid = std::sqrt(lo * hi);
        (raw(L, mid) * T < 1e-10 ? lo : hi) = mid;
      }
      L.Emin = lo;
    }
    Emin_all = std::min(Emin_all, L.Emin);
  }
}

double AdkRates::raw(const Level& L, double E) const {
  if (!(E > 0)) return 0.0;
  E = std::min(E, L.Epeak);
  return std::exp(L.logpre - L.power * std::log(E) + L.expo / E);
}

double AdkRates::rate(int k, double E) const {
  if (k < Z0 || k >= Zmax) return 0.0;
  const Level& L = lev[k];
  if (E < L.Emin) return 0.0;
  return raw(L, E);
}

double AdkRates::fieldForRate(int k, double w) const {
  const Level& L = lev[k];
  double lo = 1e-12, hi = std::min(L.Epeak, 1e6);
  if (raw(L, hi) < w) return hi;
  for (int it = 0; it < 200; ++it) {
    const double mid = std::sqrt(lo * hi);
    (raw(L, mid) < w ? lo : hi) = mid;
  }
  return hi;
}

}  // namespace gem
