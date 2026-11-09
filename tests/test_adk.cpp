// ctest adk: the ADK rates of AdkRates against the analytic standard ADK (Ammosov, Delone, Krainov 1986: l* = n* - 1,
// C^2 = 2^{2n*} / (n* Gamma(2n*)), factor 2l + 1 for m = 0), written out independently here, for the bound electrons of Rb
// (5s: l = 0, 5p_1/2 and 5p_3/2: l = 1) at several fields; hydrogen against the Landau rate;
// the nitrogen species with the element table's l; the ratios standard / l* = 0 (Rb: 0.418, 0.719, 0.708).
#include <cmath>
#include <cstdio>
#include <vector>

#include "atomic_data.hpp"
#include "ionization.hpp"
#include "units.hpp"

namespace {
const double UH = 13.605693122994, EA = 5.14220674763e11, WA = 4.1341373335e16;
// standard ADK [1/s] (Z after = 1, m = 0), E [V/m]; kept constant beyond the maximum of the expression
double adk_standard(double U, int l, double E) {
  const double ns = std::sqrt(UH / U), x = std::pow(U / UH, 1.5) * EA;
  const double C2 = std::pow(2.0, 2 * ns) / (ns * std::tgamma(2 * ns));
  const double Ep = (2.0 / 3.0) * x / (2 * ns - 1);
  E = std::min(E, Ep);
  return WA * C2 * (2 * l + 1) * (U / (2 * UH)) * std::pow(2 * x / E, 2 * ns - 1) * std::exp(-(2.0 / 3.0) * x / E);
}
double adk_lstar0(double U, double E) {
  const double ns = std::sqrt(UH / U), x = std::pow(U / UH, 1.5) * EA;
  const double C2 = std::pow(2.0, 2 * ns) / (ns * std::tgamma(ns + 1) * std::tgamma(ns));
  return WA * C2 * (U / (2 * UH)) * std::pow(2 * x / E, 2 * ns - 1) * std::exp(-(2.0 / 3.0) * x / E);
}
}  // namespace

int main() {
  int fails = 0;
  const double lam_um = 0.8, E0 = 2.0 * gem::PI * gem::MC2_MEV * 1e6 / (lam_um * 1e-6), w0 = 2.0 * gem::PI * gem::C_LIGHT / (lam_um * 1e-6);
  gem::Units u;
  u.lambda0_um = lam_um;
  const double U[3] = {4.177128, 4.177128 - 1.559591, 4.177128 - 1.589049};
  const int l[3] = {0, 1, 1};
  const char* name[3] = {"5s", "5p1/2", "5p3/2"};
  double worst = 0, worst_si = 0;
  for (int k = 0; k < 3; ++k) {
    gem::AdkRates R;
    R.setup({U[k]}, 0, 1, u, 1e3, {l[k]});
    for (double f : {0.003, 0.005, 0.01, 0.02}) {
      const double E = f * EA, ref = adk_standard(U[k], l[k], E);
      const double code = R.rate(0, E / E0) * w0;
      if (code > 0) worst = std::max(worst, std::abs(code / ref - 1));
      const double ns = std::sqrt(UH / U[k]), Ep = (2.0 / 3.0) * std::pow(U[k] / UH, 1.5) * EA / (2 * ns - 1);
      if (E < Ep) worst_si = std::max(worst_si, std::abs(gem::adkRateSI(U[k], 1, E, l[k]) / ref - 1));   // adkRateSI: no cap
    }
    std::printf("  %s: standard / (l* = 0) = %.3f\n", name[k], adk_standard(U[k], l[k], 0.01 * EA) / adk_lstar0(U[k], 0.01 * EA));
  }
  std::printf("adk: AdkRates with l against the analytic standard ADK (Rb 5s, 5p1/2, 5p3/2, E = 0.003 .. 0.02 E_a): %.1e, adkRateSI %.1e\n",
              worst, worst_si);
  if (!(worst < 1e-12 && worst_si < 1e-12)) ++fails;
  // the ion species: nitrogen with the l of the element table (removedOrbitalL: 2p, 2p, 2p, 2s, 2s, 1s, 1s) against the
  // analytic standard ADK per level, at a field well below each level's ADK maximum; the ratios to the former l* = 0 form
  {
    const gem::Element* el = gem::findElement("N");
    const std::vector<int> lN = gem::removedOrbitalL("N");
    gem::AdkRates R;
    R.setup(el->ip_eV, 0, 7, u, 1e3, lN);
    double d = 0;
    std::printf("  N, standard / (l* = 0) per level:");
    for (int k = 0; k < 7; ++k) {
      const double x = el->ip_eV[size_t(k)] / UH, ns = (k + 1) / std::sqrt(x), xx = std::pow(x, 1.5) * EA;
      const double E = 0.05 * xx / (2 * ns - 1);
      const double pw = std::pow(2 * xx / E, 2 * ns - 1) * std::exp(-(2.0 / 3.0) * xx / E) * WA * (x / 2);
      const double st = pw * std::pow(2.0, 2 * ns) / (ns * std::tgamma(2 * ns)) * (2 * lN[size_t(k)] + 1);
      const double old = pw * std::pow(2.0, 2 * ns) / (ns * std::tgamma(ns + 1) * std::tgamma(ns));
      const double code = R.rate(k, E / E0) * w0;
      if (code > 0) d = std::max(d, std::abs(code / st - 1));
      std::printf(" %d+ %.3f", k, st / old);
    }
    std::printf("\nadk: nitrogen species with the element's l against the analytic standard ADK: %.1e\n", d);
    if (!(lN.size() == 7 && lN[0] == 1 && lN[3] == 0 && lN[6] == 0 && d < 1e-12)) ++fails;
  }
  // no l given: s electrons (l = 0)
  {
    gem::AdkRates R;
    R.setup({U[0]}, 0, 1, u, 1e3);
    const double E = 0.01 * EA, d = std::abs(R.rate(0, E / E0) * w0 / adk_standard(U[0], 0, E) - 1);
    std::printf("adk: no l given = l = 0: %.1e\n", d);
    if (!(d < 1e-12)) ++fails;
  }
  // hydrogen, standard with l = 0: the Landau rate 4 omega_a (E_a / E) exp(-2 E_a / 3E)
  {
    gem::AdkRates R;
    R.setup({UH}, 0, 1, u, 1e3, {0});
    double w = 0;
    for (double f : {0.02, 0.05, 0.1}) w = std::max(w, std::abs(R.rate(0, f * EA / E0) * w0 / (4 * WA / f * std::exp(-2.0 / (3 * f))) - 1));
    std::printf("adk: hydrogen against the Landau rate: %.1e\n", w);
    if (!(w < 1e-12)) ++fails;
  }
  std::printf(fails ? "adk: FAIL\n" : "adk: ok\n");
  return fails ? 1 : 0;
}
