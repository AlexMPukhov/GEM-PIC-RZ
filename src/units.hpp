// ============================================================================
//  units.hpp -- physical constants, code (normalized) units, unit conversion
//
//  Code units (all derived from the reference wavelength lambda0):
//    length  1/k0 = lambda0/(2 pi)      time    1/omega0
//    field   m c omega0 / e             momentum m c  (u = gamma v)
//    density n_c (critical density at lambda0)
//    mass    m_e
//    charge  e n_c / k0^3
// ============================================================================
#pragma once
#include <string>

namespace gem {

constexpr double PI = 3.14159265358979323846;
constexpr double E_CHARGE_C = 1.602176634e-19;  // elementary charge [C]
constexpr double MC2_MEV = 0.51099895;           // electron rest energy [MeV]
constexpr double C_LIGHT = 299792458.0;          // speed of light [m/s]
constexpr double MP_OVER_ME = 1836.15267343;     // proton / electron mass
constexpr double U_OVER_ME = 1822.888486;        // atomic mass unit / electron mass

// Physical dimension of a deck quantity.
enum class Dim { None, Length, Time, Density, Charge, Energy, Angle, Mass, Field, Gradient };

const char* dimName(Dim d);       // "length", "time", ...
const char* dimUnitsHelp(Dim d);  // list of accepted units, for error messages

struct Units {
  double lambda0_um = 0.8;  // reference wavelength [um]

  double k0inv_um() const { return lambda0_um / (2.0 * PI); }
  double k0inv_cm() const { return lambda0_um * 1e-4 / (2.0 * PI); }
  double nc_cm3() const { return 1.1148e21 / (lambda0_um * lambda0_um); }
  double omega0_per_fs() const { return 2.0 * PI * C_LIGHT * 1e-15 / (lambda0_um * 1e-6); }
  // unit of the magnetic field m_e omega0 / e [T]
  double B0_T() const { return 9.1093837015e-31 / 1.602176634e-19 * omega0_per_fs() * 1e15; }

  // charge [C] -> code units, and back to pC
  double chargeToCode(double coulomb) const {
    const double k = k0inv_cm();
    return (coulomb / E_CHARGE_C) / (nc_cm3() * k * k * k);
  }
  double codeToPC(double q) const {
    const double k = k0inv_cm();
    return q * nc_cm3() * k * k * k * E_CHARGE_C * 1e12;
  }
};

// Information about one unit name.
struct UnitInfo {
  Dim dim = Dim::None;
  double factor = 1.0;      // to the base unit of the dimension (see units.cpp)
  bool normalized = false;  // factor converts directly to code units
};

// Look up a unit name (case-sensitive, e.g. "um", "fs", "cm^-3", "pC", "MeV").
// Returns false if the name is unknown.
bool findUnit(const std::string& name, UnitInfo& info);

// Convert a number given in `unit` (empty string: no unit) to code units:
//   Length, Time, Density, Charge -> normalized code units
//   Mass -> electron masses,  Energy -> MeV,  Angle -> rad,
//   None -> plain number ("%" allowed)
// Throws std::invalid_argument with a readable message on a missing,
// unknown or wrong-dimension unit.
double toCode(double x, const std::string& unit, Dim d, const Units& u);

}  // namespace gem
