// ============================================================================
//  units.cpp -- table of accepted units and conversion to code units
// ============================================================================
#include "units.hpp"

#include <stdexcept>

namespace gem {

namespace {
// Base units:  Length um, Time fs, Density cm^-3, Charge C, Energy MeV,
//              Angle rad, Mass m_e, None 1.   Normalized units convert directly.
struct Row {
  const char* name;
  Dim dim;
  double factor;
  bool normalized;
};
const Row TABLE[] = {
    // length
    {"nm", Dim::Length, 1e-3, false},
    {"um", Dim::Length, 1.0, false},
    {"\xC2\xB5m", Dim::Length, 1.0, false},  // "µm"
    {"micron", Dim::Length, 1.0, false},
    {"mm", Dim::Length, 1e3, false},
    {"cm", Dim::Length, 1e4, false},
    {"m", Dim::Length, 1e6, false},
    {"k0^-1", Dim::Length, 1.0, true},
    {"1/k0", Dim::Length, 1.0, true},
    {"c/w0", Dim::Length, 1.0, true},
    {"lambda0", Dim::Length, 2.0 * PI, true},
    // normalized emittance is a length (mm mrad = um)
    {"mm*mrad", Dim::Length, 1.0, false},
    {"mm.mrad", Dim::Length, 1.0, false},
    {"mm-mrad", Dim::Length, 1.0, false},
    {"mm_mrad", Dim::Length, 1.0, false},
    // time
    {"as", Dim::Time, 1e-3, false},
    {"fs", Dim::Time, 1.0, false},
    {"ps", Dim::Time, 1e3, false},
    {"ns", Dim::Time, 1e6, false},
    {"1/w0", Dim::Time, 1.0, true},
    {"T0", Dim::Time, 2.0 * PI, true},
    // density
    {"cm^-3", Dim::Density, 1.0, false},
    {"cm-3", Dim::Density, 1.0, false},
    {"1/cm3", Dim::Density, 1.0, false},
    {"/cm3", Dim::Density, 1.0, false},
    {"cm^{-3}", Dim::Density, 1.0, false},
    {"m^-3", Dim::Density, 1e-6, false},
    {"m-3", Dim::Density, 1e-6, false},
    {"1/m3", Dim::Density, 1e-6, false},
    {"n_c", Dim::Density, 1.0, true},
    {"nc", Dim::Density, 1.0, true},
    // charge
    {"fC", Dim::Charge, 1e-15, false},
    {"pC", Dim::Charge, 1e-12, false},
    {"nC", Dim::Charge, 1e-9, false},
    {"uC", Dim::Charge, 1e-6, false},
    {"C", Dim::Charge, 1.0, false},
    {"e", Dim::Charge, E_CHARGE_C, false},
    // energy (kinetic)
    {"eV", Dim::Energy, 1e-6, false},
    {"keV", Dim::Energy, 1e-3, false},
    {"MeV", Dim::Energy, 1.0, false},
    {"GeV", Dim::Energy, 1e3, false},
    {"TeV", Dim::Energy, 1e6, false},
    // angle
    {"rad", Dim::Angle, 1.0, false},
    {"mrad", Dim::Angle, 1e-3, false},
    {"urad", Dim::Angle, 1e-6, false},
    {"deg", Dim::Angle, PI / 180.0, false},
    // mass (base: electron mass, which is also the code unit)
    {"u", Dim::Mass, U_OVER_ME, true},
    {"amu", Dim::Mass, U_OVER_ME, true},
    {"m_p", Dim::Mass, MP_OVER_ME, true},
    {"m_e", Dim::Mass, 1.0, true},
    {"kg", Dim::Mass, 1.0 / 9.1093837015e-31, true},
    // dimensionless
    {"%", Dim::None, 0.01, false},
};
}  // namespace

const char* dimName(Dim d) {
  switch (d) {
    case Dim::None: return "dimensionless number";
    case Dim::Length: return "length";
    case Dim::Time: return "time";
    case Dim::Density: return "density";
    case Dim::Charge: return "charge";
    case Dim::Energy: return "energy";
    case Dim::Angle: return "angle";
    case Dim::Mass: return "mass";
  }
  return "?";
}

const char* dimUnitsHelp(Dim d) {
  switch (d) {
    case Dim::None: return "no unit, or %";
    case Dim::Length: return "nm, um, mm, cm, m; normalized: k0^-1 (= 1/k0 = c/w0), lambda0";
    case Dim::Time: return "as, fs, ps, ns; normalized: 1/w0, T0 (laser periods)";
    case Dim::Density: return "cm^-3, m^-3; normalized: n_c";
    case Dim::Charge: return "fC, pC, nC, uC, C, e (elementary charges)";
    case Dim::Energy: return "eV, keV, MeV, GeV, TeV (kinetic energy)";
    case Dim::Angle: return "rad, mrad, urad, deg";
    case Dim::Mass: return "u (= amu), m_p, m_e, kg";
  }
  return "?";
}

bool findUnit(const std::string& name, UnitInfo& info) {
  for (const auto& r : TABLE) {
    if (name == r.name) {
      info.dim = r.dim;
      info.factor = r.factor;
      info.normalized = r.normalized;
      return true;
    }
  }
  return false;
}

double toCode(double x, const std::string& unit, Dim d, const Units& u) {
  if (unit.empty()) {
    if (d == Dim::None) return x;
    throw std::invalid_argument(std::string("a ") + dimName(d) + " needs a unit (" + dimUnitsHelp(d) + ")");
  }
  UnitInfo inf;
  if (!findUnit(unit, inf))
    throw std::invalid_argument("unknown unit '" + unit + "' (for a " + dimName(d) + " use: " + dimUnitsHelp(d) + ")");
  if (inf.dim != d)
    throw std::invalid_argument("unit '" + unit + "' is a " + dimName(inf.dim) + ", but a " + dimName(d) +
                                " is expected (" + dimUnitsHelp(d) + ")");
  const double v = x * inf.factor;
  if (inf.normalized) return v;
  switch (d) {
    case Dim::Length: return v / u.k0inv_um();
    case Dim::Time: return v * u.omega0_per_fs();
    case Dim::Density: return v / u.nc_cm3();
    case Dim::Charge: return u.chargeToCode(v);
    case Dim::Energy:
    case Dim::Angle:
    case Dim::Mass:
    case Dim::None: return v;
  }
  return v;
}

}  // namespace gem
