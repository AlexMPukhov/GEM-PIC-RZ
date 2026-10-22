// ============================================================================
//  config_read_plasma.cpp -- block readers: plasma profile, plasma, ions
// ============================================================================
#include "config_internal.hpp"

namespace gem {
namespace cfgdetail {

void readProfile(const BlockReader& r, const GridConfig& g, PlasmaConfig& P) {
  (void)g;
  // longitudinal profile
  if (r.has("z_points") || r.has("z_values")) {
    P.z_points = r.vec("z_points", Dim::Length);
    P.z_values = r.vec("z_values", Dim::None);
    if (P.z_points.size() != P.z_values.size() || P.z_points.size() < 1)
      r.failKey("z_values", "must have as many entries as z_points");
    for (size_t i = 1; i < P.z_points.size(); ++i)
      if (!(P.z_points[i] > P.z_points[i - 1])) r.failKey("z_points", "must be strictly increasing");
    for (double v : P.z_values)
      if (v < 0) r.failKey("z_values", "must not be negative");
  }

  // any profile: a further factor density_function(x, y, z, r)
  if (r.has("density_function")) P.fn = r.expr("density_function", {"x", "y", "z", "r"});

  // radial profile
  const std::string rad = r.word("radial_profile", "uniform", {"uniform", "parabolic", "table"});
  P.r_max = r.num("r_max", Dim::Length, 1e300);
  P.smooth_length = r.num("smooth_length", Dim::Length, 0.0);
  if (P.smooth_length < 0.0) r.failKey("smooth_length", "must be >= 0");
  P.smooth_length_m0 = r.num("smooth_length_m0", Dim::Length, -1.0);
  if (r.has("smooth_length_m0") && P.smooth_length_m0 < 0.0) r.failKey("smooth_length_m0", "must be >= 0");
  if (rad == "uniform") {
    P.radial = PlasmaConfig::Radial::Uniform;
  } else if (rad == "parabolic") {
    P.radial = PlasmaConfig::Radial::Parabolic;
    P.channel_depth = r.num("channel_depth", Dim::None);
    P.r_channel = r.num("r_channel", Dim::Length);
    if (!(P.r_channel > 0)) r.failKey("r_channel", "must be positive");
  } else {
    P.radial = PlasmaConfig::Radial::Table;
    P.r_points = r.vec("r_points", Dim::Length);
    P.r_values = r.vec("r_values", Dim::None);
    if (P.r_points.size() != P.r_values.size() || P.r_points.size() < 2)
      r.failKey("r_values", "must have as many entries as r_points (at least 2)");
    for (size_t i = 1; i < P.r_points.size(); ++i)
      if (!(P.r_points[i] > P.r_points[i - 1])) r.failKey("r_points", "must be strictly increasing");
  }

}

void readPlasma(const Block& b, const Deck& d, const Units& u, const GridConfig& g, PlasmaConfig& P,
                std::vector<std::string>& warn) {
  BlockReader r(b, d, u);
  P.present = true;
  P.n0 = r.num("density", Dim::Density);
  if (!(P.n0 > 0)) r.failKey("density", "must be positive");
  P.free_electrons = r.flag("free_electrons", true);
  if (P.n0 > 0.1) warn.push_back("plasma.density exceeds 0.1 n_c: the scheme assumes underdense plasma");

  readProfile(r, g, P);

  // macroparticles
  P.ppc_r = int(r.integer("ppc_r", 2));
  P.ppc_phi = int(r.integer("ppc_phi", 8));
  P.ppc_r_axis = int(r.integer("ppc_r_axis", 4));
  P.ppc_phi_axis = int(r.integer("ppc_phi_axis", 16));
  P.r_axis_refine = r.num("r_axis_refine", Dim::Length, 10.0);
  if (P.ppc_r < 1 || P.ppc_phi < 1 || P.ppc_r_axis < 1 || P.ppc_phi_axis < 1)
    r.failKey("ppc_r", "particle numbers per cell must be >= 1");
  if (P.ppc_phi < 4 * (g.M + 1))
    warn.push_back("plasma.ppc_phi < 4 (m_max + 1): azimuthal modes will be noisy");
  if (P.r_max < 1e299 && P.r_max > g.r_absorber)
    warn.push_back("plasma.r_max extends into the radial absorber");
}

void readIon(const Block& b, const Deck& d, const Units& u, const PlasmaConfig& P, IonConfig& I,
             std::vector<std::string>& warn) {
  BlockReader r(b, d, u);
  // ---- identity: chemical element (NIST data) or explicit mass
  if (r.oneOf("element", "mass") == 0) {
    std::vector<std::string> syms;
    for (const auto& e : elements()) {
      std::string s = e.symbol;
      for (auto& ch : s) ch = char(std::tolower(static_cast<unsigned char>(ch)));
      syms.push_back(s);
    }
    const Element* el = findElement(r.word("element", "", syms));
    I.element = el->symbol;
    I.Znuc = el->Znuc;
    I.mass = el->mass_u * U_OVER_ME;
    I.ip_eV = el->ip_eV;
  } else {
    I.mass = r.num("mass", Dim::Mass);
    if (!(I.mass > 0)) r.failKey("mass", "must be positive");
  }
  // ---- ionization
  I.ionization = (r.word("ionization", "none", {"none", "adk"}) == "adk");
  if (r.has("ionization_potentials")) {
    std::vector<double> v = r.vec("ionization_potentials", Dim::Energy);  // MeV
    for (auto& x : v) x *= 1e6;                                            // eV
    for (size_t k = 0; k < v.size(); ++k)
      if (!(v[k] > 0) || (k > 0 && !(v[k] > v[k - 1])))
        r.failKey("ionization_potentials", "must be positive and increasing");
    if (I.Znuc > 0 && int(v.size()) != I.Znuc)
      r.failKey("ionization_potentials",
                "needs " + std::to_string(I.Znuc) + " values for " + I.element + " (charge states 0 ... " +
                    std::to_string(I.Znuc - 1) + ")");
    if (I.Znuc == 0) I.Znuc = int(v.size());
    I.ip_eV = v;
    if (!I.ionization) warn.push_back(b.label() + ": ionization_potentials given, but ionization = none");
  }
  if (I.ionization && I.ip_eV.empty())
    r.failKey("ionization", "needs ionization energies: give 'element' (NIST data) or 'ionization_potentials'");

  // ---- charge state: hydrogen isotopes and ions given by mass default to
  // Z = 1; for heavier elements it is a physics choice and must be given.
  // With ionization the default is the neutral atom.
  if (!r.has("charge_state") && !I.ionization && I.Znuc > 1)
    r.failKey("charge_state", "is required for " + I.element + " (1 ... " + std::to_string(I.Znuc) + ")");
  I.Z = int(r.integer("charge_state", I.ionization ? 0 : 1));
  if (I.Z < 0) r.failKey("charge_state", "must be >= 0");
  if (I.Z == 0 && !I.ionization) r.failKey("charge_state", "0 (neutral atoms) needs ionization = adk");
  if (I.Znuc > 0 && I.Z > I.Znuc)
    r.failKey("charge_state", "exceeds the atomic number " + std::to_string(I.Znuc) +
                                  (I.element.empty() ? std::string("") : " of " + I.element));
  if (I.ionization) {
    I.Zmax = int(r.integer("max_charge_state", I.Znuc));
    if (I.Zmax < I.Z || I.Zmax > I.Znuc)
      r.failKey("max_charge_state", "must lie between charge_state and " + std::to_string(I.Znuc));
    if (I.Zmax == I.Z) warn.push_back(b.label() + ": max_charge_state = charge_state, nothing to ionize");
  } else {
    if (r.has("max_charge_state")) r.failKey("max_charge_state", "is only used with ionization = adk");
    I.Zmax = I.Z;
  }
  I.mobile = r.flag("mobile", true);

  // ---- own density profile (a dopant in part of the plasma, a gas jet, ...)
  I.own = r.word("profile", "plasma", {"plasma", "own"}) == "own";
  if (I.own) {
    if (r.has("fraction")) r.failKey("fraction", "is not used with profile = own (give the atom 'density')");
    I.prof.present = true;
    I.prof.n0 = r.num("density", Dim::Density);
    if (!(I.prof.n0 > 0)) r.failKey("density", "must be positive");
    readProfile(r, GridConfig{}, I.prof);
    I.atoms = I.prof.n0 / P.n0;  // at the reference density, for the summary
  }
  // ---- amount: share of the neutralizing charge, or atom density
  const int which = I.own ? -1 : r.oneOf("fraction", "density", false);
  if (I.own) {
  } else if (which == 1) {
    I.atoms = r.num("density", Dim::Density) / P.n0;
    if (!(I.atoms > 0)) r.failKey("density", "must be positive");
  } else {
    if (I.Z == 0) r.failKey("charge_state", "neutral atoms: give their 'density' (a 'fraction' needs charge_state > 0)");
    const double f = (which == 0) ? r.num("fraction", Dim::None) : 1.0;
    if (!(f > 0) || f > 1.0) r.failKey("fraction", "must be in (0, 1]");
    I.atoms = f / I.Z;
  }

  const std::string def = I.element.empty() ? "ion" + std::to_string(b.index)
                          : (I.Z == 0 ? I.element : I.element + std::to_string(I.Z) + "+");
  I.name = r.str("name", def);
  if (I.mass < 100.0) warn.push_back(b.label() + " (" + I.name + "): mass below 100 m_e - is the unit right?");
}

// shot-noise model (all bunch kinds)

}  // namespace cfgdetail
}  // namespace gem
