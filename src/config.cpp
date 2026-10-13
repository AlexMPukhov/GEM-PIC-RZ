// ============================================================================
//  config.cpp -- reading the deck blocks into typed configuration structures
// ============================================================================
#include "config.hpp"

#include "atomic_data.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <sstream>

namespace gem {

// ============================================================================
//  BlockReader: typed access to the keys of one block, with units
// ============================================================================
namespace {

class BlockReader {
 public:
  BlockReader(const Block& b, const Deck& d, const Units& u) : b_(b), d_(d), u_(u) {}

  bool has(const std::string& key) const { return b_.find(key) != nullptr; }

  // required number
  double num(const std::string& key, Dim dim) const {
    const Entry* e = get(key, true);
    auto v = numbers(*e, dim);
    if (v.size() != 1) fail(e, "expects a single value, got " + std::to_string(v.size()));
    return v[0];
  }
  // optional number with default (default is already in code units)
  double num(const std::string& key, Dim dim, double def) const {
    return has(key) ? num(key, dim) : def;
  }
  std::vector<double> vec(const std::string& key, Dim dim) const {
    return numbers(*get(key, true), dim);
  }
  long integer(const std::string& key) const {
    const Entry* e = get(key, true);
    if (e->tokens.size() != 1) fail(e, "expects a single integer");
    char* end = nullptr;
    const long v = std::strtol(e->tokens[0].c_str(), &end, 10);
    if (end == e->tokens[0].c_str() || *end != '\0') fail(e, "'" + e->tokens[0] + "' is not an integer");
    return v;
  }
  long integer(const std::string& key, long def) const { return has(key) ? integer(key) : def; }

  std::string word(const std::string& key, const std::string& def, const std::vector<std::string>& allowed) const {
    if (!has(key)) return def;
    const Entry* e = get(key, true);
    if (e->tokens.size() != 1) fail(e, "expects a single word");
    std::string w = e->tokens[0];
    for (auto& ch : w) ch = char(std::tolower(static_cast<unsigned char>(ch)));
    if (std::find(allowed.begin(), allowed.end(), w) == allowed.end()) {
      std::string list;
      for (const auto& a : allowed) list += (list.empty() ? "" : ", ") + a;
      fail(e, "'" + e->tokens[0] + "' is not allowed (choose one of: " + list + ")");
    }
    return w;
  }
  std::string str(const std::string& key, const std::string& def) const {
    if (!has(key)) return def;
    const Entry* e = get(key, true);
    if (e->tokens.size() != 1) fail(e, "expects a single string (use quotes for blanks)");
    return e->tokens[0];
  }
  bool flag(const std::string& key, bool def) const {
    if (!has(key)) return def;
    const Entry* e = get(key, true);
    std::string w = e->tokens.size() == 1 ? e->tokens[0] : "";
    for (auto& ch : w) ch = char(std::tolower(static_cast<unsigned char>(ch)));
    if (w == "true" || w == "yes" || w == "on" || w == "1") return true;
    if (w == "false" || w == "no" || w == "off" || w == "0") return false;
    fail(e, "expects true/false, yes/no or on/off");
  }

  // "exactly one of a, b" -> returns 0 for a, 1 for b
  int oneOf(const std::string& a, const std::string& b, bool required = true) const {
    const bool ha = has(a), hb = has(b);
    if (ha && hb) fail(b_.find(b), "give either '" + a + "' or '" + b + "', not both");
    if (!ha && !hb) {
      if (!required) return -1;
      error(b_.line, "block '" + b_.label() + "' needs '" + a + "' or '" + b + "'");
    }
    return ha ? 0 : 1;
  }

  [[noreturn]] void fail(const Entry* e, const std::string& msg) const {
    throw DeckError(d_.where(e ? e->line : b_.line) + ": " + b_.label() + "." + (e ? e->key : std::string("?")) +
                    ": " + msg);
  }
  [[noreturn]] void error(int line, const std::string& msg) const { throw DeckError(d_.where(line) + ": " + msg); }
  [[noreturn]] void failKey(const std::string& key, const std::string& msg) const {
    const Entry* e = b_.find(key);
    if (e) fail(e, msg);
    error(b_.line, b_.label() + "." + key + ": " + msg);
  }

  const Block& block() const { return b_; }

  // raw number list and unit of an entry (no conversion)
  void raw(const std::string& key, std::vector<double>& v, std::string& unit) const {
    const Entry* e = get(key, true);
    split(*e, v, unit);
  }

 private:
  const Block& b_;
  const Deck& d_;
  const Units& u_;

  const Entry* get(const std::string& key, bool required) const {
    const Entry* e = b_.find(key);
    if (!e) {
      if (required) error(b_.line, "block '" + b_.label() + "' needs key '" + key + "'" + similarKeyHint(key));
      return nullptr;
    }
    e->used = true;
    return e;
  }

  // hint for a missing key: an unread key of the block that differs by <= 2 edits
  std::string similarKeyHint(const std::string& key) const {
    auto dist = [](const std::string& a, const std::string& b) {
      std::vector<size_t> prev(b.size() + 1), cur(b.size() + 1);
      for (size_t j = 0; j <= b.size(); ++j) prev[j] = j;
      for (size_t i = 1; i <= a.size(); ++i) {
        cur[0] = i;
        for (size_t j = 1; j <= b.size(); ++j)
          cur[j] = std::min({prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (a[i - 1] == b[j - 1] ? 0 : 1)});
        std::swap(prev, cur);
      }
      return prev[b.size()];
    };
    for (const auto& e : b_.entries)
      if (!e.used && dist(e.key, key) <= 2)
        return " (found '" + e.key + "' in line " + std::to_string(e.line) + ": misspelled?)";
    return "";
  }

  // numbers followed by an optional unit (last token, or suffix of the last number)
  void split(const Entry& e, std::vector<double>& v, std::string& unit) const {
    v.clear();
    unit.clear();
    UnitInfo inf;
    for (size_t i = 0; i < e.tokens.size(); ++i) {
      const std::string& t = e.tokens[i];
      const bool last = (i + 1 == e.tokens.size());
      if (last && !v.empty() && findUnit(t, inf)) { unit = t; break; }
      char* end = nullptr;
      const double x = std::strtod(t.c_str(), &end);
      if (end == t.c_str()) {
        if (last && !v.empty()) fail(&e, "unknown unit '" + t + "'");
        fail(&e, "'" + t + "' is not a number");
      }
      v.push_back(x);
      if (*end != '\0') {
        const std::string suffix(end);
        if (!last) fail(&e, "'" + t + "': a unit may only follow the last value");
        if (!findUnit(suffix, inf)) fail(&e, "unknown unit '" + suffix + "' in '" + t + "'");
        unit = suffix;
      }
    }
    if (v.empty()) fail(&e, "no number given");
  }

  std::vector<double> numbers(const Entry& e, Dim dim) const {
    std::vector<double> v;
    std::string unit;
    split(e, v, unit);
    try {
      for (auto& x : v) x = toCode(x, unit, dim, u_);
    } catch (const std::invalid_argument& ex) {
      fail(&e, ex.what());
    }
    return v;
  }
};

// log(cosh x), overflow-safe (same expression as the grid mapping)
double logcosh(double x) {
  const double a = std::fabs(x);
  return a + std::log1p(std::exp(-2.0 * a)) - std::log(2.0);
}

double interpTable(const std::vector<double>& xs, const std::vector<double>& ys, double x) {
  if (x <= xs.front()) return ys.front();
  if (x >= xs.back()) return ys.back();
  const size_t i = size_t(std::upper_bound(xs.begin(), xs.end(), x) - xs.begin()) - 1;
  return ys[i] + (ys[i + 1] - ys[i]) * ((x - xs[i]) / (xs[i + 1] - xs[i]));
}

}  // namespace

// ============================================================================
//  helper methods of the config structures
// ============================================================================
double GridConfig::mapping(double x) const {
  if (radial == Radial::Table) return rtable.coord(x);
  return dr_min * x + 0.5 * (dr_max - dr_min) * (x + wcells * (logcosh((x - xi0) / wcells) - logcosh(-xi0 / wcells)));
}

double GridConfig::zetaNode(int n) const { return zetaTable ? ztable.coord(n) : n * dzeta; }

double GridConfig::zetaStep(double zeta) const { return zetaTable ? ztable.step(zeta) : dzeta; }

double GridConfig::stepAt(double s) const {
  if (!zTable) return dz;
  return stable.coord(stable.index(s) + 1.0) - s;
}

int GridConfig::stepsFrom(double s) const {
  int k = 0;
  while (s < z_end - 1e-9) { s = s + stepAt(s); ++k; }  // the arithmetic of the run loop
  return k;
}

double UndulatorConfig::ku() const { return 2.0 * PI / period; }

double UndulatorConfig::resonance(double gamma) const {
  const double kk = (type == Type::Planar) ? 0.5 * K * K : K * K;
  return period * (1.0 + kk) / (2.0 * gamma * gamma);
}

double PlasmaConfig::longitudinal(double z) const {
  if (z_points.empty()) return 1.0;
  return interpTable(z_points, z_values, z);
}

double PlasmaConfig::radialFactor(double r) const {
  switch (radial) {
    case Radial::Uniform: return 1.0;
    case Radial::Parabolic: return 1.0 + channel_depth * r * r / (r_channel * r_channel);
    case Radial::Table: return interpTable(r_points, r_values, r);
  }
  return 1.0;
}

double PlasmaConfig::density(double r, double z) const {
  if (!present || r >= r_max) return 0.0;
  return n0 * longitudinal(z) * radialFactor(r);
}

// ============================================================================
//  block readers
// ============================================================================
namespace {

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

void readLaser(const Block& b, const Deck& d, const Units& u, const GridConfig& g, LaserConfig& L,
               std::vector<std::string>& warn) {
  BlockReader r(b, d, u);
  L.a0 = r.num("a0", Dim::None);
  if (r.has("wavelength")) {
    const double lam = r.num("wavelength", Dim::Length);  // code units: 2 pi lambda/lambda0
    if (!(lam > 0)) r.failKey("wavelength", "must be positive");
    L.omega = 2.0 * PI / lam;
  }
  L.waist = r.num("waist", Dim::Length);
  if (!(L.waist > 0)) r.failKey("waist", "must be positive");
  if (r.oneOf("tau_fwhm", "length_1e") == 0) {
    // intensity FWHM -> 1/e half length of the field: FWHM_I = L sqrt(2 ln 2)
    L.length_1e = r.num("tau_fwhm", Dim::Time) / std::sqrt(2.0 * std::log(2.0));
  } else {
    L.length_1e = r.num("length_1e", Dim::Length);
  }
  if (!(L.length_1e > 0)) r.failKey(r.has("tau_fwhm") ? "tau_fwhm" : "length_1e", "must be positive");
  L.zeta_center = r.num("zeta_center", Dim::Length);
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
  if (L.zeta_center + 3.0 * L.length_1e > g.window())
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

void readPlasma(const Block& b, const Deck& d, const Units& u, const GridConfig& g, PlasmaConfig& P,
                std::vector<std::string>& warn) {
  BlockReader r(b, d, u);
  P.present = true;
  P.n0 = r.num("density", Dim::Density);
  if (!(P.n0 > 0)) r.failKey("density", "must be positive");
  P.free_electrons = r.flag("free_electrons", true);
  if (P.n0 > 0.1) warn.push_back("plasma.density exceeds 0.1 n_c: the scheme assumes underdense plasma");

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

  // radial profile
  const std::string rad = r.word("radial_profile", "uniform", {"uniform", "parabolic", "table"});
  P.r_max = r.num("r_max", Dim::Length, 1e300);
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

  // ---- amount: share of the neutralizing charge, or atom density
  const int which = r.oneOf("fraction", "density", false);
  if (which == 1) {
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

void readBunch(const Block& b, const Deck& d, const Units& u, const GridConfig& g, BunchConfig& B,
               std::vector<std::string>& warn) {
  BlockReader r(b, d, u);
  B.name = r.str("name", "bunch" + std::to_string(b.index));
  const std::string sp = r.word("species", "electron", {"electron", "positron"});
  B.qm = (sp == "electron") ? -1.0 : 1.0;
  // charge: magnitude from the deck, sign from the species
  const double q = std::fabs(r.num("charge", Dim::Charge));
  if (!(q > 0)) r.failKey("charge", "must not be zero");
  B.charge = (B.qm < 0) ? -q : q;
  B.charge_pC = u.codeToPC(q);
  if (r.oneOf("energy", "gamma") == 0) {
    const double E = r.num("energy", Dim::Energy);
    B.gamma = 1.0 + E / MC2_MEV;
  } else {
    B.gamma = r.num("gamma", Dim::None);
  }
  if (!(B.gamma > 1.0)) r.failKey(r.has("gamma") ? "gamma" : "energy", "the bunch must be relativistic (gamma > 1)");
  if (B.gamma < 10.0)
    warn.push_back(b.label() + ": gamma < 10 - the bunch self-field is initialized for v = c and the beam push "
                               "assumes v_z close to c");
  B.energy_spread = r.num("energy_spread", Dim::None, 0.0);
  B.sigma_x = r.num("sigma_x", Dim::Length);
  B.sigma_y = r.num("sigma_y", Dim::Length, B.sigma_x);
  B.emittance_x = r.num("emittance_x", Dim::Length, 0.0);
  B.emittance_y = r.num("emittance_y", Dim::Length, B.emittance_x);
  if (!(B.sigma_x > 0) || !(B.sigma_y > 0)) r.failKey("sigma_x", "must be positive");
  const std::string shape = r.word("distribution_z", "gaussian", {"gaussian", "flattop"});
  if (shape == "gaussian") {
    B.shape = BunchConfig::Shape::Gaussian;
    B.sigma_z = r.num("sigma_z", Dim::Length);
    if (!(B.sigma_z > 0)) r.failKey("sigma_z", "must be positive");
  } else {
    B.shape = BunchConfig::Shape::Flattop;
    B.length = r.num("length", Dim::Length);
    if (!(B.length > 0)) r.failKey("length", "must be positive");
  }
  B.zeta_center = r.num("zeta_center", Dim::Length);
  B.x0 = r.num("x0", Dim::Length, 0.0);
  B.y0 = r.num("y0", Dim::Length, 0.0);
  B.z_waist = r.num("z_waist", Dim::Length, 0.0);
  B.macroparticles = int(r.integer("macroparticles", 10000));
  if (B.macroparticles < 1) r.failKey("macroparticles", "must be >= 1");

  const double half = (B.shape == BunchConfig::Shape::Gaussian) ? 3.0 * B.sigma_z : 0.5 * B.length;
  if (B.zeta_center - half < 0 || B.zeta_center + half > g.window())
    warn.push_back(b.label() + " (" + B.name + ") extends beyond the window (particles outside are removed)");
  // The longitudinal space-charge field of a bunch with v < c reaches the
  // distance d ahead of it with E_z(d) / E_z(head) ~ L r / (gamma d^2). The
  // scheme sets E_z = 0 at the front of the window; a field cut off there acts
  // as a uniform E_z on the whole bunch (it changes the bunch energy).
  if (B.gamma < 1e5) {
    const double L = 2.0 * half, rb = std::max(B.sigma_x, B.sigma_y);
    // (the same reach behind the bunch: the initial field is computed with psi = 0 at the rear end)
    const double front = g.window() - (B.zeta_center + half), rear = B.zeta_center - half;
    const double need = std::sqrt(100.0 * L * rb / B.gamma);  // E_z(gap) / E_z(edge) = 1 %
    if (std::min(front, rear) < need) {
      char buf[256];
      std::snprintf(buf, sizeof buf, "%.3g um from the front and %.3g um from the rear; give it >= %.3g um on both sides",
                    front * u.k0inv_um(), rear * u.k0inv_um(), need * u.k0inv_um());
      warn.push_back(b.label() + " (" + B.name + "): its space-charge field (gamma = " +
                     std::to_string(int(B.gamma)) + ") reaches the ends of the window: the bunch is " + buf);
    }
  }
  if (3.0 * std::max(B.sigma_x, B.sigma_y) > g.r_max)
    warn.push_back(b.label() + " (" + B.name + ") is wider than the radial box");
  if (std::max(B.sigma_x, B.sigma_y) < 2.0 * g.dr_min)
    warn.push_back(b.label() + " (" + B.name + "): sigma < 2 dr_min, the bunch is not resolved radially");
}

void readOutput(const Block& b, const Deck& d, const Units& u, OutputConfig& O) {
  BlockReader r(b, d, u);
  O.dir = r.str("dir", O.dir);
  O.every = int(r.integer("every", O.every));
  O.map_stride = int(r.integer("map_stride", O.map_stride));
  O.axis = r.flag("axis", O.axis);
  O.maps = r.flag("maps", O.maps);
  O.particles = r.flag("particles", O.particles);
  const std::string fmt = r.word("format", "text", {"text", "hdf5", "both"});
  O.text = (fmt != "hdf5");
  O.hdf5 = (fmt != "text");
#if !defined(GEMPIC_HAVE_HDF5)
  if (O.hdf5)
    r.failKey("format", "this build has no HDF5 support (install the HDF5 library, e.g. libhdf5-dev, and rebuild)");
#endif
  O.hdf5_single = r.word("hdf5_precision", "double", {"double", "single"}) == "single";
  O.hdf5_compression = int(r.integer("hdf5_compression", O.hdf5_compression));
  if (O.hdf5_compression < 0 || O.hdf5_compression > 9) r.failKey("hdf5_compression", "must be 0 ... 9");
  O.checkpoint_every = int(r.integer("checkpoint_every", O.checkpoint_every));
  O.checkpoint_keep = int(r.integer("checkpoint_keep", O.checkpoint_keep));
  O.walltime_hours = r.num("walltime_hours", Dim::None, O.walltime_hours);
  if (O.checkpoint_every < 0) r.failKey("checkpoint_every", "must be >= 0");
  if (O.checkpoint_keep < 0) r.failKey("checkpoint_keep", "must be >= 0");
  if (O.walltime_hours < 0) r.failKey("walltime_hours", "must be >= 0");
#if !defined(GEMPIC_HAVE_HDF5)
  if (O.checkpoint_every > 0) r.failKey("checkpoint_every", "checkpoints need HDF5 support (install libhdf5-dev and rebuild)");
  if (O.walltime_hours > 0) r.failKey("walltime_hours", "the checkpoint at the wall-time limit needs HDF5 support");
#endif
  if (O.every < 0) r.failKey("every", "must be >= 0");
  if (O.map_stride < 1) r.failKey("map_stride", "must be >= 1");
}

void readUndulator(const Block& b, const Deck& d, const Units& u, UndulatorConfig& U) {
  BlockReader r(b, d, u);
  const std::string t = r.word("type", "planar", {"planar", "helical"});
  U.type = (t == "planar") ? UndulatorConfig::Type::Planar : UndulatorConfig::Type::Helical;
  U.period = r.num("period", Dim::Length);
  if (!(U.period > 0)) r.failKey("period", "must be positive");
  U.K = r.num("k", Dim::None);  // keys are case-insensitive: K
  if (!(U.K >= 0)) r.failKey("k", "must not be negative");
  U.z_start = r.num("z_start", Dim::Length, 0.0);
  if (r.oneOf("periods", "length") == 0) {
    const long n = r.integer("periods");
    if (n < 1) r.failKey("periods", "must be >= 1");
    U.length = double(n) * U.period;
  } else {
    U.length = r.num("length", Dim::Length);
    if (!(U.length > 0)) r.failKey("length", "must be positive");
  }
  // whole periods: then the linear ramps leave no net transverse momentum (no drift of the orbit)
  const long rp = r.integer("ramp_periods", 1);
  if (rp < 0 || 2.0 * double(rp) * U.period > U.length * (1.0 + 1e-12))
    r.failKey("ramp_periods", "must be >= 0 and fit twice into the length");
  U.ramp = double(rp) * U.period;
  U.focusing = r.flag("natural_focusing", true);
  if (U.type == UndulatorConfig::Type::Helical && rp < 1)
    r.failKey("ramp_periods", "a helical undulator needs ramp_periods >= 1 (a hard edge deflects the orbit by K/gamma)");
  if (U.type == UndulatorConfig::Type::Helical && r.has("natural_focusing"))
    r.failKey("natural_focusing", "applies to planar undulators only (the helical field is the on-axis field)");
}

void readNumerics(const Block& b, const Deck& d, const Units& u, NumericsConfig& N) {
  BlockReader r(b, d, u);
#if defined(GEMPIC_HAVE_CUDA)
  const std::string bdef = "cuda";
#else
  const std::string bdef = "serial";
#endif
  const std::string be = r.word("backend", bdef, {"serial", "openmp", "cuda"});
  N.backend = be == "cuda" ? NumericsConfig::Backend::Cuda
              : be == "openmp" ? NumericsConfig::Backend::OpenMP
                               : NumericsConfig::Backend::Serial;
#if defined(GEMPIC_HAVE_CUDA)
  if (N.backend == NumericsConfig::Backend::OpenMP)
    r.failKey("backend", "this is a CUDA build: use 'cuda' or 'serial' (build without -DGEMPIC_CUDA=ON for 'openmp')");
#else
  if (N.backend == NumericsConfig::Backend::Cuda)
    r.failKey("backend", "this build has no GPU support: configure with cmake -DGEMPIC_CUDA=ON (or use 'openmp')");
#endif
  N.iterations = int(r.integer("iterations", N.iterations));
  N.beam_filter = int(r.integer("beam_filter", N.beam_filter));
  if (N.beam_filter < 0 || N.beam_filter > 8) r.failKey("beam_filter", "must be 0 ... 8");
  N.seed = static_cast<unsigned long>(r.integer("seed", long(N.seed)));
  N.beam_dz_max = r.num("beam_dz_max", Dim::Length, N.beam_dz_max);
  N.vz_min = r.num("vz_min", Dim::None, N.vz_min);
  N.vz_beam = r.num("vz_beam", Dim::None, N.vz_beam);
  if (!(N.vz_min > 0) || !(N.vz_beam < 1) || !(N.vz_min < N.vz_beam))
    r.failKey("vz_beam", "need 0 < vz_min < vz_beam < 1");
  if (N.iterations < 1) r.failKey("iterations", "must be >= 1");
  if (!(N.beam_dz_max > 0)) r.failKey("beam_dz_max", "must be positive");
}

}  // namespace

// ============================================================================
Config Config::fromDeck(const Deck& deck, std::vector<std::string>* warnings) {
  std::vector<std::string> warn;
  deck.checkBlockNames({"grid", "laser", "plasma", "ion", "bunch", "undulator", "output", "numerics"});
  Config c;
  const Block* grid = deck.single("grid");
  if (!grid) throw DeckError(deck.source + ": the deck needs a 'grid' block");
  c.units.lambda0_um = readLambda0(deck, grid);

  readGrid(*grid, deck, c.units, c.grid, warn);
  for (const Block* b : deck.all("laser")) {
    LaserConfig L;
    readLaser(*b, deck, c.units, c.grid, L, warn);
    c.lasers.push_back(L);
  }
  if (!c.lasers.empty() && c.grid.M < 1)
    throw DeckError(deck.where(grid->line) + ": a laser needs azimuthal mode 1: set grid.m_max >= 1");
  if (const Block* b = deck.single("plasma")) readPlasma(*b, deck, c.units, c.grid, c.plasma, warn);
  double fsum = 0.0;
  for (const Block* b : deck.all("ion")) {
    if (!c.plasma.present)
      throw DeckError(deck.where(b->line) + ": an 'ion' block needs a 'plasma' block (ions follow its profile)");
    if (c.ions.size() >= 10) throw DeckError(deck.where(b->line) + ": at most 10 ion blocks are supported");
    IonConfig I;
    readIon(*b, deck, c.units, c.plasma, I, warn);
    fsum += I.fraction();
    c.ions.push_back(I);
  }
  if (!c.plasma.free_electrons && fsum > 0.0) {
    throw DeckError(deck.where(deck.single("plasma")->line) +
                    ": plasma.free_electrons = false, but ion blocks with charge_state > 0 carry positive charge "
                    "(it must be neutralized by free electrons)");
  }
  if (fsum > 1.0 + 1e-9) {
    const auto v = deck.all("ion");
    throw DeckError(deck.where(v.back()->line) + ": the ion fractions add up to " + std::to_string(fsum) +
                    " > 1 (fraction = charge_state x atom density / free-electron density)");
  }
  if (c.plasma.present && !c.plasma.free_electrons && c.ions.empty())
    warn.push_back("plasma.free_electrons = false and no ion blocks: the plasma block has no effect");
  for (const Block* b : deck.all("bunch")) {
    BunchConfig B;
    readBunch(*b, deck, c.units, c.grid, B, warn);
    c.bunches.push_back(B);
  }
  for (const Block* b : deck.all("undulator")) {
    UndulatorConfig U;
    readUndulator(*b, deck, c.units, U);
    c.undulators.push_back(U);
  }
  if (const Block* b = deck.single("output")) readOutput(*b, deck, c.units, c.output);
  if (const Block* b = deck.single("numerics")) readNumerics(*b, deck, c.units, c.numerics);
  deck.checkAllUsed();
  // resolution of an undulator: the steps along z and the sub-steps of the
  // beam push must follow the wiggling (period / 16), the zeta grid the
  // radiation (resonant wavelength / 12)
  for (size_t k = 0; k < c.undulators.size(); ++k) {
    const auto& U = c.undulators[k];
    const std::string lab = "undulator[" + std::to_string(k) + "]";
    if (c.grid.dz_max > U.period / 16.0)
      warn.push_back(lab + ": grid.dz (up to " + std::to_string(c.grid.dz_max / U.period) +
                     " periods) should be <= period / 16 (the transverse current oscillates with the period)");
    if (c.numerics.beam_dz_max > U.period / 16.0)
      warn.push_back(lab + ": numerics.beam_dz_max should be <= period / 16 (sub-steps of the beam push)");
    for (const auto& B : c.bunches) {
      const double lr = U.resonance(B.gamma);
      if (c.grid.dzeta_max > lr / 12.0)
        warn.push_back(lab + ", bunch '" + B.name + "': the resonant wavelength is only " +
                       std::to_string(lr / c.grid.dzeta_max) + " zeta steps (use >= 12, better ~20)");
    }
  }

  if (c.lasers.empty() && c.bunches.empty()) warn.push_back("the deck defines neither a laser nor a bunch");
  if (!c.plasma.present) warn.push_back("no plasma block: vacuum propagation");
  if (warnings) *warnings = warn;
  return c;
}

// ============================================================================
std::string Config::summary() const {
  const double um = units.k0inv_um();
  const double fs = 1.0 / units.omega0_per_fs();
  char buf[512];
  std::ostringstream os;
  auto line = [&](const char* fmt, auto... args) {
    std::snprintf(buf, sizeof buf, fmt, args...);
    os << buf << '\n';
  };
  line("units: lambda0 = %.4g um, 1/k0 = %.5g um, 1/omega0 = %.5g fs, n_c = %.4g cm^-3", units.lambda0_um, um, fs,
       units.nc_cm3());
  const auto& g = grid;
  line("grid : modes m = 0..%d", g.M);
  if (g.radial == GridConfig::Radial::Table) {
    std::string pts;
    for (size_t k = 0; k < g.rtable.points().size(); ++k) {
      std::snprintf(buf, sizeof buf, "%s(%.4g, %.4g)", k ? " " : "", g.rtable.points()[k] * um,
                    g.rtable.steps()[k] * um);
      pts += buf;
    }
    line("       radial table: Nr = %d, dr = %.4g .. %.4g um, r_max = %.4g um, max cell-to-cell change %.1f %%",
         g.Nr, g.dr_min * um, g.dr_max * um, g.r_max * um, 100 * g.rtable.maxRelativeChange());
    line("       reference points (r, dr) [um]: %s", pts.c_str());
  } else if (g.stretched)
    line("       radial stretched: Nr = %d, dr = %.4g .. %.4g um, transition at r = %.4g um (cell %.2f, width %.3g "
         "cells), r_max = %.4g um",
         g.Nr, g.dr_min * um, g.dr_max * um, g.mapping(g.xi0) * um, g.xi0, g.wcells, g.r_max * um);
  else
    line("       radial uniform: Nr = %d, dr = %.4g um, r_max = %.4g um", g.Nr, g.dr_min * um, g.r_max * um);
  auto points = [&](const StepMapping& m, double f) {
    std::string pts;
    for (size_t k = 0; k < m.points().size(); ++k) {
      std::snprintf(buf, sizeof buf, "%s(%.4g, %.4g)", k ? " " : "", m.points()[k] * f, m.steps()[k] * um);
      pts += buf;
    }
    return pts;
  };
  if (g.zetaTable) {
    line("       window table: Nz = %d, dzeta = %.4g .. %.4g um, length = %.4g um, max step-to-step change %.1f %%",
         g.Nz, g.dzeta * um, g.dzeta_max * um, g.window() * um, 100 * g.ztable.maxRelativeChange());
    line("       reference points (zeta, dzeta) [um]: %s", points(g.ztable, um).c_str());
  } else {
    line("       window: Nz = %d, dzeta = %.4g um, length = %.4g um", g.Nz, g.dzeta * um, g.window() * um);
  }
  if (g.zTable) {
    line("       z table: dz = %.4g .. %.4g um, z_end = %.4g mm (%d steps), max step-to-step change %.1f %%",
         g.dz * um, g.dz_max * um, g.z_end * um * 1e-3, g.stepsFrom(0.0), 100 * g.stable.maxRelativeChange());
    line("       reference points (z [mm], dz [um]): %s", points(g.stable, um * 1e-3).c_str());
    line("       absorber r > %.4g um, strength %.3g", g.r_absorber * um, g.absorber_strength);
  } else {
    line("       dz = %.4g um, z_end = %.4g mm (%d steps); absorber r > %.4g um, strength %.3g", g.dz * um,
         g.z_end * um * 1e-3, g.stepsFrom(0.0), g.r_absorber * um, g.absorber_strength);
  }
  for (size_t i = 0; i < lasers.size(); ++i) {
    const auto& L = lasers[i];
    // peak power, linear polarization: P = 21.5 GW (a0 w0 / lambda)^2; circular: x2
    const double lam_um = units.lambda0_um / L.omega;
    double P_TW = 21.5e-3 * std::pow(L.a0 * L.waist * um / lam_um, 2.0);
    if (L.pol == LaserConfig::Pol::Circular) P_TW *= 2.0;
    line("laser[%zu]: a0 = %.4g, lambda = %.4g um, w0 = %.4g um (zR = %.4g mm), tau_FWHM = %.4g fs, "
         "zeta_c = %.4g um, focus at z = %.4g mm, %s, P = %.3g TW",
         i, L.a0, lam_um, L.waist * um, L.zR() * um * 1e-3, L.length_1e * std::sqrt(2.0 * std::log(2.0)) * fs,
         L.zeta_center * um, L.z_focus * um * 1e-3, L.pol == LaserConfig::Pol::Linear ? "linear" : "circular",
         P_TW);
  }
  for (size_t k = 0; k < undulators.size(); ++k) {
    const auto& U = undulators[k];
    std::string res;
    for (const auto& B : bunches) {
      std::snprintf(buf, sizeof buf, "%s%s %.4g um", res.empty() ? "" : ", ", B.name.c_str(),
                    U.resonance(B.gamma) * um);
      res += buf;
    }
    line("undulator[%zu]: %s, period = %.4g mm, K = %.4g (B0 = %.4g T), z = %.4g ... %.4g mm (%.4g periods, "
         "ramps %.3g periods)%s",
         k, U.type == UndulatorConfig::Type::Planar ? "planar" : "helical", U.period * um * 1e-3, U.K,
         U.B0() * units.B0_T(), U.z_start * um * 1e-3, (U.z_start + U.length) * um * 1e-3, U.length / U.period,
         U.ramp / U.period,
         U.type == UndulatorConfig::Type::Planar ? (U.focusing ? ", natural focusing" : ", no focusing") : "");
    if (!res.empty()) line("        resonant wavelength (on axis): %s", res.c_str());
  }
  if (plasma.present) {
    const double kp = std::sqrt(plasma.n0);
    char rbuf[64] = "open (up to the grid edge)";
    if (plasma.r_max < 1e299) std::snprintf(rbuf, sizeof rbuf, "%.4g um", plasma.r_max * um);
    line("plasma: n0 = %.4g cm^-3 (%.4g n_c), lambda_p = %.4g um, radial profile: %s, r_max = %s",
         plasma.n0 * units.nc_cm3(), plasma.n0, 2.0 * PI / kp * um,
         plasma.radial == PlasmaConfig::Radial::Uniform     ? "uniform"
         : plasma.radial == PlasmaConfig::Radial::Parabolic ? "parabolic channel"
                                                            : "table",
         rbuf);
    if (!plasma.z_points.empty()) {
      std::string s;
      for (size_t i = 0; i < plasma.z_points.size(); ++i) {
        std::snprintf(buf, sizeof buf, "%s(%.4g mm, %.3g)", i ? " " : "", plasma.z_points[i] * um * 1e-3,
                      plasma.z_values[i]);
        s += buf;
      }
      line("        longitudinal profile: %s", s.c_str());
    }
    line("        macroparticles per cell: %d x %d, and %d x %d for r < %.4g um", plasma.ppc_r, plasma.ppc_phi,
         plasma.ppc_r_axis, plasma.ppc_phi_axis, plasma.r_axis_refine * um);
  }
  if (plasma.present) {
    if (!plasma.free_electrons) line("%s", "        no free electrons at the start (neutral gas)");
    double fsum = 0.0;
    bool anyIonization = false;
    for (size_t i = 0; i < ions.size(); ++i) {
      const auto& I = ions[i];
      fsum += I.fraction();
      line("ion[%zu] '%s': Z = %d, mass = %.5g u (%.6g m_e), atom density %.4g cm^-3 (%.4g n0), %s", i,
           I.name.c_str(), I.Z, I.mass / U_OVER_ME, I.mass, I.atoms * plasma.n0 * units.nc_cm3(), I.atoms,
           I.mobile ? "mobile" : "fixed");
      if (I.ionization) {
        anyIonization = true;
        std::string lv;
        for (int k = I.Z; k < I.Zmax; ++k) {
          std::snprintf(buf, sizeof buf, "%s%d->%d: %.6g", k > I.Z ? ", " : "", k, k + 1, I.ip_eV[k]);
          lv += buf;
        }
        line("        ADK field ionization up to Z = %d, ionization energies [eV] %s", I.Zmax, lv.c_str());
      }
    }
    if (plasma.free_electrons && fsum < 1.0 - 1e-9)
      line("        immobile ion background: %.3g of the positive charge", 1.0 - fsum);
    if (anyIonization) line("        ionization energies: %s", atomicDataSource());
  }
  for (size_t i = 0; i < bunches.size(); ++i) {
    const auto& B = bunches[i];
    line("bunch[%zu] '%s': %s, Q = %.4g pC, E = %.5g MeV (gamma %.5g), dE/E = %.3g, sigma_x,y = %.4g, %.4g um, "
         "eps_n = %.4g, %.4g um, %s %.4g um at zeta = %.4g um, %d macroparticles",
         i, B.name.c_str(), B.qm < 0 ? "electrons" : "positrons", B.charge_pC, (B.gamma - 1.0) * MC2_MEV, B.gamma,
         B.energy_spread, B.sigma_x * um, B.sigma_y * um, B.emittance_x * um, B.emittance_y * um,
         B.shape == BunchConfig::Shape::Gaussian ? "sigma_z" : "length",
         (B.shape == BunchConfig::Shape::Gaussian ? B.sigma_z : B.length) * um, B.zeta_center * um,
         B.macroparticles);
  }
  if (output.hdf5)
    line("output format: %s (HDF5 %s precision, compression %d)", output.text ? "text + hdf5" : "hdf5",
         output.hdf5_single ? "single" : "double", output.hdf5_compression);
  if (output.checkpoint_every > 0 || output.walltime_hours > 0)
    line("checkpoints: every %d steps (keep %d)%s", output.checkpoint_every, output.checkpoint_keep,
         output.walltime_hours > 0 ? (", stop at " + std::to_string(output.walltime_hours) + " h wall time").c_str()
                                   : "");
  line("output: dir '%s', every %d steps; numerics: %d iterations/slice, beam filter %d, seed %lu, backend %s",
       output.dir.c_str(), output.every, numerics.iterations, numerics.beam_filter, numerics.seed,
       numerics.backend == NumericsConfig::Backend::Cuda     ? "cuda"
       : numerics.backend == NumericsConfig::Backend::OpenMP ? "openmp"
                                                             : "serial");
  return os.str();
}

}  // namespace gem
