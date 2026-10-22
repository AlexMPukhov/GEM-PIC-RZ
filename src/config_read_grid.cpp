// ============================================================================
//  config_read_grid.cpp -- block readers: lambda0, grid, constants, laser, flying focus
// ============================================================================
#include "config_internal.hpp"

namespace gem {
namespace cfgdetail {

double readLambda0(const Deck& deck, const Block* grid) {
  // grid.lambda0, or the wavelength of the first laser; physical units only
  const Block* src = nullptr;
  std::string key;
  if (grid && grid->find("lambda0")) { src = grid; key = "lambda0"; }
  else {
    auto lasers = deck.all("laser");
    if (!lasers.empty() && lasers[0]->find("wavelength")) { src = lasers[0]; key = "wavelength"; }
  }
  if (!src)
    throw DeckError(deck.where(grid ? grid->line : 1) +
                    ": set the reference wavelength grid.lambda0 (it defines the code units)");
  Units dummy;
  BlockReader r(*src, deck, dummy);
  std::vector<double> v;
  std::string unit;
  r.raw(key, v, unit);
  UnitInfo inf;
  if (v.size() != 1) r.failKey(key, "expects a single value");
  if (unit.empty() || !findUnit(unit, inf) || inf.dim != Dim::Length || inf.normalized)
    r.failKey(key, "must be given in physical length units (nm, um, ...)");
  const double um = v[0] * inf.factor;
  if (!(um > 0)) r.failKey(key, "must be positive");
  return um;
}

void readGrid(const Block& b, const Deck& d, const Units& u, GridConfig& g, std::vector<std::string>& warn) {
  BlockReader r(b, d, u);
  if (r.has("lambda0")) r.num("lambda0", Dim::Length);  // marks as used (already read)
  g.M = int(r.integer("m_max", 1));
  if (g.M < 0 || g.M > 16) r.failKey("m_max", "must be between 0 and 16");

  // ---- radial
  const std::string type = r.word("radial_grid", "stretched", {"uniform", "stretched", "table"});
  g.radial = type == "uniform" ? GridConfig::Radial::Uniform
             : type == "table" ? GridConfig::Radial::Table
                               : GridConfig::Radial::Stretched;
  g.stretched = (g.radial == GridConfig::Radial::Stretched);
  if (g.radial == GridConfig::Radial::Table) {
    // free-hand grid: steps dr_points at the radii r_points, linear in between
    const std::vector<double> rp = r.vec("r_points", Dim::Length);
    const std::vector<double> hp = r.vec("dr_points", Dim::Length);
    if (rp.size() != hp.size()) r.failKey("dr_points", "must have as many entries as r_points");
    const std::string err = g.rtable.set(rp, hp);
    if (!err.empty()) r.failKey("r_points", err);
    g.dr_min = g.rtable.minStep();
    g.dr_max = g.rtable.maxStep();
    if (g.rtable.maxRelativeChange() > 0.25)
      warn.push_back("grid: the radial step changes by up to " +
                     std::to_string(int(std::lround(100 * g.rtable.maxRelativeChange()))) +
                     " % from one cell to the next; more than ~20 % reduces the accuracy");
  } else if (g.stretched) {
    g.dr_min = r.num("dr_min", Dim::Length);
    g.dr_max = r.num("dr_max", Dim::Length);
    if (!(g.dr_min > 0) || g.dr_max < g.dr_min) r.failKey("dr_max", "need 0 < dr_min <= dr_max");
    g.wcells = r.num("transition_cells", Dim::None, 8.0);
    if (!(g.wcells > 0)) r.failKey("transition_cells", "must be positive");
    if (r.oneOf("r_transition", "transition_cell") == 1) {
      g.xi0 = r.num("transition_cell", Dim::None);
    } else {
      const double rt = r.num("r_transition", Dim::Length);
      if (!(rt > 0)) r.failKey("r_transition", "must be positive");
      // solve R(xi0) = rt for xi0 (R(xi0) is monotonic in xi0) by bisection
      double lo = 0.0, hi = rt / g.dr_min + 1.0;
      for (int it = 0; it < 200; ++it) {
        const double mid = 0.5 * (lo + hi);
        GridConfig t = g;
        t.xi0 = mid;
        (t.mapping(mid) < rt ? lo : hi) = mid;
      }
      g.xi0 = 0.5 * (lo + hi);
    }
  } else {
    g.dr_min = g.dr_max = r.num("dr", Dim::Length);
    if (!(g.dr_min > 0)) r.failKey("dr", "must be positive");
    g.xi0 = 0.0;
    g.wcells = 8.0;
  }
  if (r.oneOf("nr", "r_max") == 0) {
    g.Nr = int(r.integer("nr"));
  } else {
    const double rmax = r.num("r_max", Dim::Length);
    int N = 1;
    while (g.mapping(N) < rmax * (1.0 - 1e-12) && N < 1000000) ++N;
    g.Nr = N;
  }
  if (g.Nr < 4) r.failKey(r.has("nr") ? "nr" : "r_max", "the radial grid needs at least 4 cells");
  g.r_max = g.mapping(g.Nr);

  // ---- longitudinal window
  // free-hand table (shared with the z grid): steps hk at the points xk, linear in between
  auto readTable = [&](const char* pk, const char* hk, StepMapping& m, const char* what) {
    const std::vector<double> xp = r.vec(pk, Dim::Length);
    const std::vector<double> hp = r.vec(hk, Dim::Length);
    if (xp.size() != hp.size()) r.failKey(hk, std::string("must have as many entries as ") + pk);
    const std::string err = m.set(xp, hp);
    if (!err.empty()) r.failKey(pk, err);
    if (m.maxRelativeChange() > 0.25)
      warn.push_back(std::string("grid: the ") + what + " step changes by up to " +
                     std::to_string(int(std::lround(100 * m.maxRelativeChange()))) +
                     " % from one step to the next; more than ~20 % reduces the accuracy");
  };
  g.zetaTable = r.word("zeta_grid", "uniform", {"uniform", "table"}) == "table";
  if (g.zetaTable) {
    readTable("zeta_points", "dzeta_points", g.ztable, "zeta");
    g.dzeta = g.ztable.minStep();
    g.dzeta_max = g.ztable.maxStep();
    if (r.oneOf("window_length", "nzeta") == 1) {
      g.Nz = int(r.integer("nzeta"));
    } else {
      const double L = r.num("window_length", Dim::Length);
      if (!(L > 0)) r.failKey("window_length", "must be positive");
      int N = 1;  // rounded up to a whole cell, like r_max
      while (g.ztable.coord(N) < L * (1.0 - 1e-12) && N < 10000000) ++N;
      g.Nz = N;
    }
    if (g.Nz < 4) r.failKey(r.has("nzeta") ? "nzeta" : "window_length", "the window needs at least 4 slices");
  } else {
    g.dzeta = g.dzeta_max = r.num("dzeta", Dim::Length);
    if (!(g.dzeta > 0)) r.failKey("dzeta", "must be positive");
    if (r.oneOf("window_length", "nzeta") == 1) {
      g.Nz = int(r.integer("nzeta"));
    } else {
      const double L = r.num("window_length", Dim::Length);
      g.Nz = int(std::lround(L / g.dzeta));
      if (std::fabs(g.Nz * g.dzeta - L) > 1e-6 * L)
        warn.push_back("grid.window_length is not a multiple of dzeta; using " + std::to_string(g.Nz) + " slices");
    }
    if (g.Nz < 4) r.failKey("dzeta", "the window needs at least 4 slices");
  }

  // ---- propagation
  g.zTable = r.word("z_grid", "uniform", {"uniform", "table"}) == "table";
  if (g.zTable) {
    readTable("z_points", "dz_points", g.stable, "propagation");
    g.dz = g.stable.minStep();
    g.dz_max = g.stable.maxStep();
  } else {
    g.dz = g.dz_max = r.num("dz", Dim::Length);
    if (!(g.dz > 0)) r.failKey("dz", "must be positive");
  }
  g.z_end = r.num("z_end", Dim::Length);
  if (!(g.z_end > 0)) r.failKey("z_end", "must be positive");

  // ---- absorber
  g.r_absorber = r.num("r_absorber", Dim::Length, 0.85 * g.r_max);
  g.absorber_strength = r.num("absorber_strength", Dim::None, 0.2);
  if (g.r_absorber >= g.r_max) warn.push_back("grid.r_absorber >= r_max: the radial absorber is off");
}

// block 'constants': named values (numbers with units or expressions of the
// constants above them), in code units, e.g.  w0 = 8 um;  zr = pi * w0^2 / (0.8 um)
void readConstants(const Block& b, const Deck& d, const Units& u, std::map<std::string, double>& C) {
  BlockReader r(b, d, u);
  static const char* reserved[] = {"pi", "x", "y", "z", "r", "zeta", "if"};
  for (const auto& e : b.entries) {
    for (const char* w : reserved)
      if (e.key == w) r.fail(&e, "'" + e.key + "' is reserved (a variable or built-in name)");
    if (e.key.empty() || !(std::isalpha(static_cast<unsigned char>(e.key[0])) || e.key[0] == '_'))
      r.fail(&e, "a constant name must start with a letter");
    C[e.key] = r.expr(e.key, {})(nullptr);  // "8 um", "3e18 cm^-3", "2 * w0": code units
  }
}

void readLaser(const Block& b, const Deck& d, const Units& u, const GridConfig& g, LaserConfig& L,
               std::vector<std::string>& warn) {
  BlockReader r(b, d, u);
  L.a0 = r.num("a0", Dim::None);
  if (r.has("wavelength")) {
    const double lam = r.num("wavelength", Dim::Length);  // code units: 2 pi lambda/lambda0
    if (!(lam > 0)) r.failKey("wavelength", "must be positive");
    L.omega = 2.0 * PI / lam;
  }
  // transverse profile: Gaussian beam (waist, z_focus) or profile_function(r) at the start
  if (r.has("profile_function")) {
    if (r.has("waist")) r.failKey("waist", "give either 'waist' or 'profile_function', not both");
    if (r.has("z_focus"))
      r.failKey("z_focus", "profile_function is the field at the start (z = 0): no focus position can be applied");
    L.profile_fn = r.expr("profile_function", {"r"});
    // 1/e radius of |g| for the checks (relative to the largest value on [0, 1e4])
    double gmax = 0.0;
    for (int k = 0; k <= 4000; ++k) gmax = std::max(gmax, std::fabs(L.profile_fn({0.25 * k})));
    if (!(gmax > 0) || !std::isfinite(gmax)) r.failKey("profile_function", "must be finite and not zero everywhere");
    L.waist = 0.0;
    for (int k = 0; k <= 4000; ++k)
      if (std::fabs(L.profile_fn({0.25 * k})) >= gmax / std::exp(1.0)) L.waist = 0.25 * k;
    L.waist = std::max(L.waist, 0.25);
  } else {
    L.waist = r.num("waist", Dim::Length);
    if (!(L.waist > 0)) r.failKey("waist", "must be positive");
  }
  // longitudinal envelope: Gaussian (tau_fwhm or length_1e) or envelope_function(zeta)
  if (r.has("envelope_function")) {
    if (r.has("tau_fwhm") || r.has("length_1e"))
      r.failKey("envelope_function", "replaces tau_fwhm / length_1e: give only one of them");
    L.envelope_fn = r.expr("envelope_function", {"zeta"});
  } else {
    if (r.oneOf("tau_fwhm", "length_1e") == 0) {
      // intensity FWHM -> 1/e half length of the field: FWHM_I = L sqrt(2 ln 2)
      L.length_1e = r.num("tau_fwhm", Dim::Time) / std::sqrt(2.0 * std::log(2.0));
    } else {
      L.length_1e = r.num("length_1e", Dim::Length);
    }
    if (!(L.length_1e > 0)) r.failKey(r.has("tau_fwhm") ? "tau_fwhm" : "length_1e", "must be positive");
  }
  if (L.envelope_fn.empty()) {
    L.zeta_center = r.num("zeta_center", Dim::Length);
  } else {
    // envelope_function is a function of the window coordinate zeta itself;
    // zeta_center (optional) is then only the reference of the carrier phase
    const double W = g.window();
    double fmax = 0.0, s0 = 0.0, s1 = 0.0, s2 = 0.0;
    for (int k = 0; k <= 4000; ++k) {
      const double z = W * k / 4000.0, f = L.envelope_fn({z});
      if (!std::isfinite(f)) r.failKey("envelope_function", "is not finite in the window");
      fmax = std::max(fmax, std::fabs(f));
      s0 += f * f; s1 += f * f * z; s2 += f * f * z * z;
    }
    if (!(fmax > 0)) r.failKey("envelope_function", "is zero everywhere in the window");
    // centre and 1/e length of the intensity-weighted envelope (for the checks below)
    const double zc = s1 / s0;
    L.length_1e = std::max(2.0 * std::sqrt(std::max(0.0, s2 / s0 - zc * zc)), 1e-3);
    L.zeta_center = r.num("zeta_center", Dim::Length, zc);
    if (std::fabs(L.envelope_fn({W})) > 1e-3 * fmax)
      warn.push_back(b.label() + ": envelope_function is not small at the front of the window");
  }
  L.z_focus = r.num("z_focus", Dim::Length, 0.0);
  L.cep = r.num("cep", Dim::Angle, 0.0);
  const std::string pol = r.word("polarization", "linear", {"linear", "circular"});
  L.pol = (pol == "linear") ? LaserConfig::Pol::Linear : LaserConfig::Pol::Circular;
  L.pol_angle = r.num("pol_angle", Dim::Angle, 0.0);
  L.helicity = int(r.integer("helicity", 1));
  if (L.helicity != 1 && L.helicity != -1) r.failKey("helicity", "must be +1 or -1");
  if (L.pol == LaserConfig::Pol::Circular && r.has("pol_angle"))
    r.failKey("pol_angle", "is only used for linear polarization");

  if (L.zeta_center < 0 || L.zeta_center > g.window())
    r.failKey("zeta_center", "lies outside the window [0, window_length]");
  if (L.envelope_fn.empty() && L.zeta_center + 3.0 * L.length_1e > g.window())
    warn.push_back(b.label() + ": the pulse front is closer than 3 length_1e to the front of the window");
  if (2.0 * L.waist > g.r_max) warn.push_back(b.label() + ": 2 x waist exceeds the radial box r_max");
  double hmax = g.zetaStep(L.zeta_center);  // largest step over the pulse (+- 2 length_1e)
  for (int k = -20; k <= 20; ++k) hmax = std::max(hmax, g.zetaStep(L.zeta_center + 0.1 * k * L.length_1e));
  if (g.zetaTable)
    for (double x : g.ztable.points())
      if (std::fabs(x - L.zeta_center) < 2.0 * L.length_1e) hmax = std::max(hmax, g.zetaStep(x));
  if (hmax * L.omega > 0.5)
    warn.push_back(b.label() + ": fewer than ~12 points per laser wavelength in zeta (dzeta too large)");
}



void readFlyingFocus(const Block& b, const Deck& d, const Units& u, const GridConfig& g, FlyingFocusConfig& F,
                     std::vector<std::string>& warn) {
  BlockReader r(b, d, u);
  F.present = true;
  F.a0 = r.num("a0", Dim::None);
  F.theta = r.num("cone_angle", Dim::Angle);
  if (!(F.theta > 0 && F.theta < 0.5)) r.failKey("cone_angle", "must be between 0 and 0.5 rad");
  if (r.oneOf("tau_fwhm", "length_1e") == 0)
    F.length_1e = r.num("tau_fwhm", Dim::Time) / std::sqrt(2.0 * std::log(2.0));
  else
    F.length_1e = r.num("length_1e", Dim::Length);
  if (!(F.length_1e > 0)) r.failKey(r.has("tau_fwhm") ? "tau_fwhm" : "length_1e", "must be positive");
  F.zeta_center = r.num("zeta_center", Dim::Length);
  if (F.zeta_center < 0 || F.zeta_center > g.window()) r.failKey("zeta_center", "lies outside the window");
  F.cep = r.num("cep", Dim::Angle, 0.0);
  F.background = r.num("background_density", Dim::Density, 0.0);
  if (F.background < 0) r.failKey("background_density", "must be >= 0");
  if (2.0 * F.kappa() - F.kappa() * F.kappa() <= F.background)
    r.failKey("background_density", "too high for this cone angle: no conical wave (k_perp^2 < 0); "
                                    "needs omega_p^2 < omega^2 sin^2(theta)");
  if (g.M < 1) r.error(b.line, "a flying focus needs azimuthal mode 1: set grid.m_max >= 1");
  if (!(g.r_absorber < g.r_max))
    r.error(b.line, "the flying focus is fed through the radial absorber layer: set grid.r_absorber < r_max");
  // the conical wave has the transverse wavelength 2 pi / (k sin theta) at every radius
  const double lperp = 2.0 * PI / std::sin(F.theta);
  if (g.dr_max > lperp / 10.0)
    warn.push_back(b.label() + ": dr_max = " + std::to_string(g.dr_max / lperp).substr(0, 5) +
                   " transverse wavelengths 2 pi / (k sin theta): use dr <= 1/10 of it everywhere");
  // the rings ahead of the axial focus: at radius r the converging wave is ~ r theta / 2 ahead in zeta
  if (F.zeta_center + 3.0 * F.length_1e + 0.5 * g.r_max * F.theta > g.window())
    warn.push_back(b.label() + ": the converging wave at r_max (r theta / 2 ahead) and the pulse need room in front "
                   "of zeta_center: window too short");
}

}  // namespace cfgdetail
}  // namespace gem
