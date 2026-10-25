// ============================================================================
//  config.cpp -- reading the deck blocks into typed configuration structures
// ============================================================================
#include "config.hpp"

#include "atomic_data.hpp"
#include "beam_import.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <sstream>

namespace gem {

// ============================================================================
//  BlockReader: typed access to the keys of one block, with units
// ============================================================================
namespace {

// constants of the deck (block 'constants'), visible to every expression
const std::map<std::string, double>* g_consts = nullptr;
const std::map<std::string, double> NO_CONSTS;

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
  // a function of the variables `vars` (a quoted expression, see expr.hpp)
  Expr expr(const std::string& key, const std::vector<std::string>& vars) const {
    const Entry* e = get(key, true);
    std::string text;
    for (const auto& t : e->tokens) text += (text.empty() ? "" : " ") + t;
    try {
      Expr x(text, vars, g_consts ? *g_consts : NO_CONSTS, u_);
      std::vector<double> probe(vars.size(), 0.0);
      (void)x(probe.data());
      return x;
    } catch (const std::invalid_argument& ex) {
      fail(e, std::string(ex.what()) +
                  (e->tokens.size() > 1 ? " (put expressions with commas in quotes: \"...\")" : ""));
    }
  }
  bool hasExpr(const std::string& key) const { return has(key); }

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

  // a number, optionally with a unit suffix ("5um") - anything else is an expression
  static bool isPlainNumber(const std::string& t) {
    char* end = nullptr;
    std::strtod(t.c_str(), &end);
    if (end == t.c_str()) return false;
    UnitInfo inf;
    return *end == '\0' || findUnit(end, inf);
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
    // a single value may be an expression with units and constants, e.g. "2 * w0"
    // (its result is taken in code units, like the deck's own conversions)
    if (e.tokens.size() == 1 && !isPlainNumber(e.tokens[0])) {
      try {
        Expr x(e.tokens[0], {}, g_consts ? *g_consts : NO_CONSTS, u_);
        return {x(nullptr)};
      } catch (const std::invalid_argument& ex) {
        fail(&e, ex.what());
      }
    }
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
  while (s < z_end - endTol()) { s = s + stepAt(s); ++k; }  // the arithmetic of the run loop
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

double PlasmaConfig::density(double x, double y, double z) const {
  const double r = std::sqrt(x * x + y * y);
  const double n = density(r, z);
  if (n <= 0.0 || fn.empty()) return n;
  const double v[4] = {x, y, z, r};
  return std::max(0.0, n * fn(v));
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

// laser from a LASY / openPMD envelope file (laser.file)
void readLaserFromFile(const Block& b, const Deck& d, const Units& u, const GridConfig& g, LaserConfig& L,
                       std::vector<std::string>& warn) {
  BlockReader r(b, d, u);
  for (const char* k : {"waist", "tau_fwhm", "length_1e", "z_focus", "cep", "polarization", "pol_angle", "helicity",
                        "profile_function", "envelope_function", "wavelength"})
    if (r.has(k)) r.failKey(k, "is given by the laser file (remove it, or remove 'file')");
  std::string path = r.str("file", "");
  if (!std::filesystem::exists(path)) {  // relative to the directory of the deck
    const std::filesystem::path alt = std::filesystem::path(d.source).parent_path() / path;
    if (std::filesystem::exists(alt)) path = alt.string();
  }
  auto D = std::make_shared<LaserFileData>();
  try {
    *D = readLaserFile(path, r.str("file_mesh", ""), r.integer("file_iteration", -1), u);
  } catch (const std::exception& e) {
    r.failKey("file", e.what());
  }
  L.omega = D->omegaL;
  // centroid and rms duration of the intensity (all radii, weight r), for t_center and the checks
  double s0 = 0, s1 = 0, s2 = 0;
  const size_t inner = D->cartesian ? size_t(D->ny) * D->nx : size_t(D->nr);
  const int outer = D->cartesian ? 1 : D->ncomp;
  for (int c = 0; c < outer; ++c)
    for (int l = 0; l < D->nt; ++l)
      for (size_t q = 0; q < inner; ++q) {
        double w = std::norm(D->env[(size_t(c) * D->nt + l) * inner + q]);
        if (!D->cartesian) w *= std::max(D->r0 + q * D->dr, 0.25 * D->dr);
        const double t = D->t0 + l * D->dt;
        s0 += w; s1 += w * t; s2 += w * t * t;
      }
  const double tc = s1 / s0;
  L.length_1e = std::max(2.0 * std::sqrt(std::max(0.0, s2 / s0 - tc * tc)), 1e-3);
  L.t_center = r.has("t_center") ? r.num("t_center", Dim::Time) : tc;
  L.zeta_center = r.num("zeta_center", Dim::Length);
  if (L.zeta_center < 0 || L.zeta_center > g.window()) r.failKey("zeta_center", "lies outside the window");
  L.file_z = r.num("file_z", Dim::Length, 0.0);
  L.boundary = r.word("inject", "box", {"box", "boundary"}) == "boundary";
  L.background = r.num("background_density", Dim::Density, 0.0);
  if (L.background < 0) r.failKey("background_density", "must be >= 0");
  if (L.background >= L.omega * L.omega) r.failKey("background_density", "is overcritical for the laser");
  L.n_theta = int(r.integer("n_theta", 0));
  if (L.n_theta < 0) r.failKey("n_theta", "must be >= 0");
  if (L.n_theta > 0 && !D->cartesian) r.failKey("n_theta", "is only used for cartesian files");
  L.time_padding = int(r.integer("time_padding", 2));
  if (L.time_padding < 1 || L.time_padding > 16) r.failKey("time_padding", "must be 1 ... 16");
  // amplitude: as in the file, or rescaled so that the peak |E| is a0 omega_L
  if (r.has("a0")) {
    L.a0 = r.num("a0", Dim::None);
    if (!(L.a0 > 0)) r.failKey("a0", "must be positive");
    L.file_scale = L.a0 * L.omega / D->peak;
  } else {
    L.a0 = D->peak / L.omega;
  }
  L.waist = 0.0;  // not used for file lasers
  if (L.boundary) {
    if (!(g.r_absorber < g.r_max))
      r.failKey("inject", "boundary feeds the field through the radial absorber layer: set grid.r_absorber < r_max");
    if (std::fabs(g.absorber_strength - 1.0) > 1e-12)
      warn.push_back(b.label() + ": inject = boundary works best with grid.absorber_strength = 1 (R = incident "
                                 "field at r_max)");
  }
  // window: the file's time span at zeta_center
  const double front = L.zeta_center + (L.t_center - D->t0), rear = L.zeta_center - (D->t0 + (D->nt - 1) * D->dt - L.t_center);
  if (front > g.window() || rear < 0)
    warn.push_back(b.label() + ": the time span of the file does not fit into the window at zeta_center (it covers "
                   "zeta = " + std::to_string(rear * u.k0inv_um()).substr(0, 7) + " ... " +
                   std::to_string(front * u.k0inv_um()).substr(0, 7) + " um); the part outside is cut");
  double hmax = 0.0;
  for (int k = -20; k <= 20; ++k) hmax = std::max(hmax, g.zetaStep(L.zeta_center + 0.1 * k * L.length_1e));
  if (hmax * L.omega > 0.5)
    warn.push_back(b.label() + ": fewer than ~12 points per laser wavelength in zeta (dzeta too large)");
  const double Rf = D->cartesian ? std::min({-D->x0, D->x0 + (D->nx - 1) * D->dx, -D->y0, D->y0 + (D->ny - 1) * D->dy})
                                 : D->r0 + (D->nr - 1) * D->dr;
  if (L.boundary && Rf < g.r_max)
    warn.push_back(b.label() + ": the file ends at r = " + std::to_string(Rf * u.k0inv_um()).substr(0, 7) +
                   " um, inside the box: the incident field beyond is zero (the file should extend past r_max)");
  L.file = D;
}

void readLaser(const Block& b, const Deck& d, const Units& u, const GridConfig& g, LaserConfig& L,
               std::vector<std::string>& warn) {
  BlockReader r(b, d, u);
  if (r.has("file")) {
    readLaserFromFile(b, d, u, g, L, warn);
    return;
  }
  for (const char* k : {"file_mesh", "file_iteration", "t_center", "file_z", "inject", "background_density", "n_theta",
                        "time_padding"})
    if (r.has(k)) r.failKey(k, "is only used with 'file' (a laser from a LASY / openPMD envelope file)");
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
void readShotNoise(const BlockReader& r, BunchConfig& B) {
  B.fawley = r.word("shot_noise", "random", {"random", "fawley"}) == "fawley";
  if (B.fawley) {
    B.noise_wavelength = r.num("noise_wavelength", Dim::Length);
    if (!(B.noise_wavelength > 0)) r.failKey("noise_wavelength", "must be positive");
    B.beamlet = int(r.integer("beamlet", 4));
    if (B.beamlet < 3) r.failKey("beamlet", "must be >= 3 (macroparticles per beamlet)");
    if (B.macroparticles % B.beamlet != 0)
      r.failKey("macroparticles", "must be a multiple of beamlet = " + std::to_string(B.beamlet) +
                                      " with shot_noise = fawley");
  } else {
    if (r.has("noise_wavelength")) r.failKey("noise_wavelength", "is used with shot_noise = fawley");
    if (r.has("beamlet")) r.failKey("beamlet", "is used with shot_noise = fawley");
  }
  B.axis_refine = int(r.integer("axis_refine", 1));
  if (B.axis_refine < 1 || B.axis_refine > 1000) r.failKey("axis_refine", "must be 1 ... 1000");
  B.axis_radius = r.num("axis_radius", Dim::Length, 0.0);
  if (B.axis_radius < 0) r.failKey("axis_radius", "must be >= 0 (0: the first four radial cells)");
  if (B.axis_refine == 1 && r.has("axis_radius")) r.failKey("axis_radius", "is used with axis_refine > 1");
}

// bunch read from a file (bunch.file): select, place at zeta_center, rescale charge
void readBunchFile(const Block& b, const Deck& d, const Units& u, const GridConfig& g, BunchConfig& B,
                   std::vector<std::string>& warn) {
  BlockReader r(b, d, u);
  B.shape = BunchConfig::Shape::File;
  std::string path = r.str("file", "");
  if (!std::filesystem::exists(path)) {  // relative to the directory of the deck
    const std::filesystem::path alt = std::filesystem::path(d.source).parent_path() / path;
    if (std::filesystem::exists(alt)) path = alt.string();
  }
  ImportedBeam imp;
  try {
    imp = importBeam(path, r.str("file_species", ""), r.integer("file_iteration", -1));
  } catch (const std::exception& e) {
    r.failKey("file", e.what());
  }
  const double L = u.k0inv_um() * 1e-6;  // code length [m]
  const double Qc = 1.602176634e-19 / (u.codeToPC(1.0) * 1e-12);  // one electron in code charge units
  // longitudinal coordinate relative to the charge centroid of the species
  double sw = 0.0, sz = 0.0;
  for (const auto& p : imp.p) { sw += p.w; sz += p.w * p.z; }
  const double zc = sz / sw;
  Expr sel;
  if (r.has("select")) sel = r.expr("select", {"x", "y", "z", "ux", "uy", "uz", "gamma"});
  std::vector<BunchConfig::FileParticle> v;
  v.reserve(imp.p.size());
  for (const auto& p : imp.p) {
    const BunchConfig::FileParticle f{p.x / L, p.y / L, (p.z - zc) / L, p.ux, p.uy, p.uz, p.w * Qc};
    if (!sel.empty()) {
      const double vars[7] = {f.x, f.y, f.zeta, f.ux, f.uy, f.uz, std::sqrt(1.0 + f.ux * f.ux + f.uy * f.uy + f.uz * f.uz)};
      if (sel(vars) == 0.0) continue;
    }
    v.push_back(f);
  }
  if (v.empty()) r.failKey("select", "no particle of the file passes the selection");
  // place the centroid of the selected particles at zeta_center (plus x0, y0)
  B.x0 = r.num("x0", Dim::Length, 0.0);
  B.y0 = r.num("y0", Dim::Length, 0.0);
  const double zcen = r.num("zeta_center", Dim::Length);
  double Q = 0.0, mz = 0.0;
  for (const auto& f : v) { Q += f.q; mz += f.q * f.zeta; }
  mz /= Q;
  double scale = 1.0;
  if (r.has("charge")) {
    const double q = std::fabs(r.num("charge", Dim::Charge));
    if (!(q > 0)) r.failKey("charge", "must not be zero");
    scale = q / Q;
  }
  for (auto& f : v) { f.zeta += zcen - mz; f.x += B.x0; f.y += B.y0; f.q *= scale; }
  Q *= scale;
  // statistics for the summary, the checks and the self-fields
  double g1 = 0, x1 = 0, x2 = 0, y1 = 0, y2 = 0, z2 = 0, zmin = 1e300, zmax = -1e300;
  double xu = 0, ux1 = 0, ux2 = 0, yu = 0, uy1 = 0, uy2 = 0;
  for (const auto& f : v) {
    const double gm = std::sqrt(1.0 + f.ux * f.ux + f.uy * f.uy + f.uz * f.uz), w = f.q / Q;
    g1 += w * gm; x1 += w * f.x; x2 += w * f.x * f.x; y1 += w * f.y; y2 += w * f.y * f.y;
    z2 += w * (f.zeta - zcen) * (f.zeta - zcen); ux1 += w * f.ux; ux2 += w * f.ux * f.ux; xu += w * f.x * f.ux;
    uy1 += w * f.uy; uy2 += w * f.uy * f.uy; yu += w * f.y * f.uy;
    zmin = std::min(zmin, f.zeta); zmax = std::max(zmax, f.zeta);
  }
  double g2 = 0;
  for (const auto& f : v) {
    const double gm = std::sqrt(1.0 + f.ux * f.ux + f.uy * f.uy + f.uz * f.uz);
    g2 += f.q / Q * (gm - g1) * (gm - g1);
  }
  B.qm = imp.charge_sign;
  B.charge = imp.charge_sign * Q;
  B.charge_pC = u.codeToPC(Q);
  B.gamma = g1;
  if (!(B.gamma > 1.0)) r.failKey("file", "the mean gamma of the beam is not > 1");
  B.energy_spread = std::sqrt(g2) / g1;
  B.sigma_x = std::max(1e-12, std::sqrt(std::max(0.0, x2 - x1 * x1)));
  B.sigma_y = std::max(1e-12, std::sqrt(std::max(0.0, y2 - y1 * y1)));
  B.sigma_z = std::max(1e-12, std::sqrt(z2));
  const double cx = x2 - x1 * x1, cu = ux2 - ux1 * ux1, cxu = xu - x1 * ux1;
  B.emittance_x = std::sqrt(std::max(0.0, cx * cu - cxu * cxu));
  const double cy = y2 - y1 * y1, cv = uy2 - uy1 * uy1, cyv = yu - y1 * uy1;
  B.emittance_y = std::sqrt(std::max(0.0, cy * cv - cyv * cyv));
  B.file_center = zcen;
  B.zeta_min = zmin;
  B.zeta_max = zmax;
  B.length = std::max(zmax - zmin, 1e-12);
  B.zeta_center = 0.5 * (zmin + zmax);  // the checks use the extent
  B.macroparticles = int(r.integer("macroparticles", long(v.size())));
  if (B.macroparticles < 1) r.failKey("macroparticles", "must be >= 1");
  B.resample_smoothing = r.num("resample_smoothing", Dim::None, 0.0);
  if (B.resample_smoothing < 0) r.failKey("resample_smoothing", "must be >= 0");
  readShotNoise(r, B);
  if (B.axis_refine > 1 && !B.fawley && B.macroparticles == int(v.size()))
    r.failKey("axis_refine", "needs resampling: set bunch.macroparticles different from the number of particles in "
                             "the file (the file particles as they are cannot be refined)");
  B.file_particles.swap(v);
  char buf[512];
  std::snprintf(buf, sizeof buf, "%s, %s '%s'%s: %zu of %zu particles%s", path.c_str(), imp.format.c_str(),
                imp.species.c_str(), imp.iteration >= 0 ? (" iteration " + std::to_string(imp.iteration)).c_str() : "",
                B.file_particles.size(), imp.p.size(), sel.empty() ? "" : (" (select: " + sel.text() + ")").c_str());
  B.file_info = buf;
  if (B.macroparticles != int(B.file_particles.size()) || B.fawley)
    B.file_info += ", resampled to " + std::to_string(B.macroparticles) + " macroparticles";
  if (B.zeta_min < 0 || B.zeta_max > g.window())
    warn.push_back(b.label() + " (" + B.name + ") extends beyond the window (particles outside are removed)");
  if (std::max(B.sigma_x, B.sigma_y) < 2.0 * g.dr_min)
    warn.push_back(b.label() + " (" + B.name + "): sigma < 2 dr_min, the bunch is not resolved radially");
}

void readBunch(const Block& b, const Deck& d, const Units& u, const GridConfig& g, BunchConfig& B,
               std::vector<std::string>& warn) {
  BlockReader r(b, d, u);
  B.name = r.str("name", "bunch" + std::to_string(b.index));
  if (r.has("file")) {
    readBunchFile(b, d, u, g, B, warn);
    return;
  }
  const std::string sp = r.word("species", "electron", {"electron", "positron"});
  B.qm = (sp == "electron") ? -1.0 : 1.0;
  // charge: magnitude from the deck, sign from the species
  const double q = std::fabs(r.num("charge", Dim::Charge));
  if (!(q > 0)) r.failKey("charge", "must not be zero");
  B.charge = (B.qm < 0) ? -q : q;
  B.charge_pC = u.codeToPC(q);
  const int eg = r.oneOf("energy", "gamma", !r.has("gamma_function"));
  if (eg == 0) {
    const double E = r.num("energy", Dim::Energy);
    B.gamma = 1.0 + E / MC2_MEV;
  } else if (eg == 1) {
    B.gamma = r.num("gamma", Dim::None);
  } else {
    B.gamma = 2.0;  // gamma_function: replaced below by the mean over the profile
  }
  if (!(B.gamma > 1.0)) r.failKey(r.has("gamma") ? "gamma" : "energy", "the bunch must be relativistic (gamma > 1)");
  if (eg >= 0 && B.gamma < 10.0)
    warn.push_back(b.label() + ": gamma < 10 - the beam push along s assumes v_z close to c");
  B.energy_spread = r.num("energy_spread", Dim::None, 0.0);
  const std::string shape_r = r.word("distribution_r", "gaussian", {"gaussian", "function"});
  if (shape_r == "gaussian") {
    B.sigma_x = r.num("sigma_x", Dim::Length);
    B.sigma_y = r.num("sigma_y", Dim::Length, B.sigma_x);
    if (!(B.sigma_x > 0) || !(B.sigma_y > 0)) r.failKey("sigma_x", "must be positive");
  } else {
    if (r.has("sigma_x") || r.has("sigma_y"))
      r.failKey("sigma_x", "is not used with distribution_r = function (the profile_r gives the size)");
    B.profile_r = r.expr("profile_r", {"x", "y", "r", "zeta"});
    B.r_extent = r.num("r_extent", Dim::Length);
    if (!(B.r_extent > 0)) r.failKey("r_extent", "must be positive");
  }
  B.emittance_x = r.num("emittance_x", Dim::Length, 0.0);
  B.emittance_y = r.num("emittance_y", Dim::Length, B.emittance_x);
  if (r.has("ux_function")) B.ux_fn = r.expr("ux_function", {"x", "y", "r", "zeta"});
  if (r.has("uy_function")) B.uy_fn = r.expr("uy_function", {"x", "y", "r", "zeta"});
  const std::string shape = r.word("distribution_z", "gaussian", {"gaussian", "flattop", "function"});
  if (shape == "gaussian") {
    B.shape = BunchConfig::Shape::Gaussian;
    B.sigma_z = r.num("sigma_z", Dim::Length);
    if (!(B.sigma_z > 0)) r.failKey("sigma_z", "must be positive");
  } else if (shape == "flattop") {
    B.shape = BunchConfig::Shape::Flattop;
    B.length = r.num("length", Dim::Length);
    if (!(B.length > 0)) r.failKey("length", "must be positive");
  } else {
    // current profile f(zeta) (any normalization) on [zeta_min, zeta_max]
    B.shape = BunchConfig::Shape::Function;
    B.profile_z = r.expr("profile_z", {"zeta"});
    B.zeta_min = r.num("zeta_min", Dim::Length);
    B.zeta_max = r.num("zeta_max", Dim::Length);
    if (!(B.zeta_max > B.zeta_min)) r.failKey("zeta_max", "must be larger than zeta_min");
    double sum = 0.0;
    for (int k = 0; k <= 1000; ++k) {
      const double v = B.profile_z({B.zeta_min + (B.zeta_max - B.zeta_min) * k / 1000.0});
      if (v < 0 || !std::isfinite(v)) r.failKey("profile_z", "must be finite and >= 0 on [zeta_min, zeta_max]");
      sum += v;
    }
    if (!(sum > 0)) r.failKey("profile_z", "is zero everywhere on [zeta_min, zeta_max]");
    B.zeta_center = 0.5 * (B.zeta_min + B.zeta_max);  // for messages and the box advice
    B.length = B.zeta_max - B.zeta_min;
  }
  if (r.has("gamma_function")) {
    if (shape != "function") r.failKey("gamma_function", "needs distribution_z = function");
    B.gamma_fn = r.expr("gamma_function", {"zeta"});
    for (int k = 0; k <= 100; ++k)
      if (!(B.gamma_fn({B.zeta_min + (B.zeta_max - B.zeta_min) * k / 100.0}) > 1.0))
        r.failKey("gamma_function", "must be > 1 on [zeta_min, zeta_max]");
    if (r.has("energy") || r.has("gamma"))
      r.failKey("gamma_function", "replaces 'energy' / 'gamma': give only one of them");
    // the bunch's gamma (self-fields, box advice, messages): mean over the current profile
    double sw = 0.0, sg = 0.0;
    for (int k = 0; k <= 1000; ++k) {
      const double z = B.zeta_min + (B.zeta_max - B.zeta_min) * k / 1000.0, w = B.profile_z({z});
      sw += w;
      sg += w * B.gamma_fn({z});
    }
    B.gamma = sg / sw;
    if (B.gamma < 10.0)
      warn.push_back(b.label() + ": gamma < 10 - the beam push along s assumes v_z close to c");
  }
  if (shape != "function") B.zeta_center = r.num("zeta_center", Dim::Length);
  if (!B.profile_r.empty()) {
    // bound of f for the rejection sampling and the rms sizes (sigma_x, sigma_y: messages, checks, box advice),
    // from a polar table over the disk at 21 positions along the bunch
    const double zc = B.zeta_center;
    const double hz = (B.shape == BunchConfig::Shape::Gaussian) ? 3.0 * B.sigma_z : 0.5 * B.length;
    const int nrr = 2000, nph = 64, nzz = 21;
    double fmax = 0.0, w = 0.0, wx = 0.0, wy = 0.0;
    for (int kz = 0; kz < nzz; ++kz) {
      const double z = zc - hz + 2.0 * hz * kz / (nzz - 1);
      for (int i = 0; i < nrr; ++i) {
        const double rr = B.r_extent * (i + 0.5) / nrr;
        for (int j = 0; j < nph; ++j) {
          const double ph = 2.0 * PI * j / nph, x = rr * std::cos(ph), y = rr * std::sin(ph);
          const double f = B.profile_r({x, y, rr, z});
          if (!(f >= 0) || !std::isfinite(f)) r.failKey("profile_r", "must be finite and >= 0 in the disk r < r_extent");
          fmax = std::max(fmax, f);
          w += f * rr;
          wx += f * rr * x * x;
          wy += f * rr * y * y;
        }
      }
    }
    if (!(fmax > 0)) r.failKey("profile_r", "is zero everywhere in the disk r < r_extent");
    B.profile_r_max = 1.1 * fmax;  // margin for peaks between the table points; checked while sampling
    B.sigma_x = std::max(1e-12, std::sqrt(wx / w));
    B.sigma_y = std::max(1e-12, std::sqrt(wy / w));
  }
  B.x0 = r.num("x0", Dim::Length, 0.0);
  B.y0 = r.num("y0", Dim::Length, 0.0);
  B.z_waist = r.num("z_waist", Dim::Length, 0.0);
  B.macroparticles = int(r.integer("macroparticles", 10000));
  if (B.macroparticles < 1) r.failKey("macroparticles", "must be >= 1");
  readShotNoise(r, B);

  const double half = (B.shape == BunchConfig::Shape::Gaussian) ? 3.0 * B.sigma_z : 0.5 * B.length;
  if (B.zeta_center - half < 0 || B.zeta_center + half > g.window())
    warn.push_back(b.label() + " (" + B.name + ") extends beyond the window (particles outside are removed)");
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
  // format: one or more of text, hdf5, openpmd (both = text, hdf5)
  if (r.has("format")) {
    O.text = O.hdf5 = O.openpmd = false;
    const Entry* e = b.find("format");
    for (std::string w : e->tokens) {
      for (auto& ch : w) ch = char(std::tolower(static_cast<unsigned char>(ch)));
      if (w == "text") O.text = true;
      else if (w == "hdf5") O.hdf5 = true;
      else if (w == "openpmd") O.openpmd = true;
      else if (w == "both") O.text = O.hdf5 = true;
      else r.fail(e, "'" + w + "' is not allowed (choose one or more of: text, hdf5, openpmd, both)");
    }
    e->used = true;
  }
#if !defined(GEMPIC_HAVE_HDF5)
  if (O.hdf5 || O.openpmd)
    r.failKey("format", "this build has no HDF5 support (install the HDF5 library, e.g. libhdf5-dev, and rebuild)");
#endif
  // openPMD: fields interpolated to a uniform r-z grid (default: the smallest steps)
  O.openpmd_dr = r.num("openpmd_dr", Dim::Length, 0.0);
  O.openpmd_dzeta = r.num("openpmd_dzeta", Dim::Length, 0.0);
  if (O.openpmd_dr < 0) r.failKey("openpmd_dr", "must be positive");
  if (O.openpmd_dzeta < 0) r.failKey("openpmd_dzeta", "must be positive");
  const std::string og = r.word("openpmd_grid", "both", {"uniform", "native", "both"});
  O.openpmd_uniform = (og != "native");
  O.openpmd_native = (og != "uniform");
  if (!O.openpmd && r.has("openpmd_grid")) r.failKey("openpmd_grid", "needs format = openpmd");
  if (!O.openpmd && (r.has("openpmd_dr") || r.has("openpmd_dzeta")))
    r.failKey(r.has("openpmd_dr") ? "openpmd_dr" : "openpmd_dzeta", "needs format = openpmd");
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
  const std::string fel = r.word("fel", "auto", {"auto", "true", "false", "yes", "no", "on", "off"});
  O.fel_mode = (fel == "auto") ? -1 : (fel == "true" || fel == "yes" || fel == "on") ? 1 : 0;
  O.fel_wavelength = r.num("fel_wavelength", Dim::Length, 0.0);
  if (O.fel_wavelength < 0) r.failKey("fel_wavelength", "must be positive");
  O.fel_slice = r.num("fel_slice", Dim::Length, 0.0);
  if (O.fel_slice < 0) r.failKey("fel_slice", "must be positive");
  O.fel_harmonics = int(r.integer("fel_harmonics", 1));
  if (O.fel_harmonics < 1 || O.fel_harmonics > 9) r.failKey("fel_harmonics", "must be between 1 and 9");
  if (O.map_stride < 1) r.failKey("map_stride", "must be >= 1");
  // LASY envelope files of the laser
  O.lasy = r.flag("lasy", false);
  O.lasy_every = int(r.integer("lasy_every", 0));
  O.lasy_dr = r.num("lasy_dr", Dim::Length, 0.0);
  O.lasy_dt = r.num("lasy_dt", Dim::Time, 0.0);
  O.lasy_r_max = r.num("lasy_r_max", Dim::Length, 0.0);
  if (r.has("lasy_wavelength")) {
    const double lam = r.num("lasy_wavelength", Dim::Length);
    if (!(lam > 0)) r.failKey("lasy_wavelength", "must be positive");
    O.lasy_omega = 2.0 * PI / lam;
  }
  O.lasy_band = r.num("lasy_bandwidth", Dim::None, 0.5);
  for (const char* k : {"lasy_every", "lasy_dr", "lasy_dt", "lasy_r_max", "lasy_wavelength", "lasy_bandwidth"})
    if (!O.lasy && r.has(k)) r.failKey(k, "needs lasy = true");
  if (O.lasy_every < 0) r.failKey("lasy_every", "must be >= 0");
  if (O.lasy_dr < 0 || O.lasy_dt < 0 || O.lasy_r_max < 0) r.failKey("lasy", "grid steps must be positive");
  if (!(O.lasy_band > 0 && O.lasy_band < 1)) r.failKey("lasy_bandwidth", "must be between 0 and 1");
#if !defined(GEMPIC_HAVE_HDF5)
  if (O.lasy) r.failKey("lasy", "LASY files need HDF5 support (install libhdf5-dev and rebuild)");
#endif
}

void readUndulator(const Block& b, const Deck& d, const Units& u, UndulatorConfig& U) {
  BlockReader r(b, d, u);
  const std::string t = r.word("type", "planar", {"planar", "helical"});
  U.type = (t == "planar") ? UndulatorConfig::Type::Planar : UndulatorConfig::Type::Helical;
  U.period = r.num("period", Dim::Length);
  if (!(U.period > 0)) r.failKey("period", "must be positive");
  U.z_start = r.num("z_start", Dim::Length, 0.0);
  if (r.oneOf("periods", "length") == 0) {
    const long n = r.integer("periods");
    if (n < 1) r.failKey("periods", "must be >= 1");
    U.length = double(n) * U.period;
  } else {
    U.length = r.num("length", Dim::Length);
    if (!(U.length > 0)) r.failKey("length", "must be positive");
  }
  // keys are case-insensitive: K, K_function
  if (r.oneOf("k", "k_function") == 0) {
    U.K = r.num("k", Dim::None);
    if (!(U.K >= 0)) r.failKey("k", "must not be negative");
  } else {
    // taper: K(z) of the lab position z (the ramps multiply it as usual)
    U.K_fn = r.expr("k_function", {"z"});
    for (int k = 0; k <= 1000; ++k) {
      const double v = U.K_fn({U.z_start + U.length * k / 1000.0});
      if (!(v >= 0) || !std::isfinite(v)) r.failKey("k_function", "must be finite and >= 0 inside the undulator");
    }
    U.K = U.K_fn({U.z_start});  // entrance value: resonance and messages
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

void readMagnet(const Block& b, const Deck& d, const Units& u, MagnetConfig& G, std::vector<std::string>& warn) {
  BlockReader r(b, d, u);
  G.label = b.label();
  G.kind = b.name == "quadrupole" ? MagnetConfig::Kind::Quadrupole
           : b.name == "plasma_lens" ? MagnetConfig::Kind::PlasmaLens : MagnetConfig::Kind::Solenoid;
  G.z_start = r.num("z_start", Dim::Length);
  G.length = r.num("length", Dim::Length);
  if (!(G.length > 0)) r.failKey("length", "must be positive");
  if (G.kind == MagnetConfig::Kind::Solenoid) {
    G.field = r.num("field", Dim::Field);
    G.ramp = r.num("ramp", Dim::Length, 0.0);
    if (G.ramp < 0 || 2.0 * G.ramp > G.length) r.failKey("ramp", "must be >= 0 and fit twice into the length");
  } else {
    G.gradient = r.num("gradient", Dim::Gradient);
    if (G.kind == MagnetConfig::Kind::PlasmaLens) {
      G.radius = r.num("radius", Dim::Length, 1e300);
      if (!(G.radius > 0)) r.failKey("radius", "must be positive");
    }
  }
  (void)warn;
}

void readChicane(const Block& b, const Deck& d, const Units& u, ChicaneConfig& C) {
  BlockReader r(b, d, u);
  C.z = r.num("z", Dim::Length);
  if (!(C.z > 0)) r.failKey("z", "must be > 0 (the chicane acts when the beam passes z)");
  C.r56 = r.num("r56", Dim::Length);
  C.delay = r.num("delay", Dim::Length, 0.5 * std::fabs(C.r56));
  if (r.oneOf("energy", "gamma") == 0) C.gamma_ref = 1.0 + r.num("energy", Dim::Energy) / MC2_MEV;
  else C.gamma_ref = r.num("gamma", Dim::None);
  if (!(C.gamma_ref > 1.0)) r.failKey(r.has("gamma") ? "gamma" : "energy", "the reference must be relativistic");
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
  const std::string td = r.word("tridiag", "auto", {"auto", "thomas", "partition"});
  N.tridiag = td == "thomas" ? NumericsConfig::Tridiag::Thomas
              : td == "partition" ? NumericsConfig::Tridiag::Partition : NumericsConfig::Tridiag::Auto;
  N.iterations = int(r.integer("iterations", N.iterations));
  N.beam_filter = int(r.integer("beam_filter", N.beam_filter));
  if (N.beam_filter < 0 || N.beam_filter > 8) r.failKey("beam_filter", "must be 0 ... 8");
  N.beam_smooth_length = r.num("beam_smooth_length", Dim::Length, 0.0);
  if (N.beam_smooth_length < 0) r.failKey("beam_smooth_length", "must be >= 0");
  {
    const std::string bf = r.word("beam_field", "auto", {"auto", "level", "mean"});
    N.beam_field = bf == "level" ? BEAM_FIELD_LEVEL : bf == "mean" ? BEAM_FIELD_MEAN : -1;
  }
  N.beam_vay = r.word("beam_pusher", "vay", {"vay", "boris"}) == "vay";
  N.seed = static_cast<unsigned long>(r.integer("seed", long(N.seed)));
  N.beam_dz_max = r.num("beam_dz_max", Dim::Length, N.beam_dz_max);
  N.vz_min = r.num("vz_min", Dim::None, N.vz_min);
  N.vz_beam = r.num("vz_beam", Dim::None, N.vz_beam);
  N.trap_advance = r.num("trap_advance", Dim::None, N.trap_advance);
  if (!(N.trap_advance >= 1.0)) r.failKey("trap_advance", "must be >= 1 (in steps dz)");
  N.trap_gamma = r.num("trap_gamma", Dim::None, N.trap_gamma);
  if (!(N.trap_gamma >= 0.0)) r.failKey("trap_gamma", "must be >= 0");
  N.theta = r.num("theta", Dim::None, N.theta);
  if (!(N.theta >= 0.5 && N.theta <= 1.0)) r.failKey("theta", "must be in [0.5, 1]");
  N.theta_stiff = r.num("theta_stiff", Dim::None, N.theta_stiff);
  if (!(N.theta_stiff >= 0.0)) r.failKey("theta_stiff", "must be >= 0");
  N.beam_consistent = r.flag("beam_consistent", N.beam_consistent);
  if (!(N.vz_min > 0) || !(N.vz_beam < 1) || !(N.vz_min < N.vz_beam))
    r.failKey("vz_beam", "need 0 < vz_min < vz_beam < 1");
  if (N.iterations < 1) r.failKey("iterations", "must be >= 1");
  if (!(N.beam_dz_max > 0)) r.failKey("beam_dz_max", "must be positive");
}

}  // namespace

// ============================================================================
//  box size for bunches with v < c
//
//  Radial: the field of a bunch of length L reaches the radius ~ gamma L (rest
//  frame: rod of length gamma L). A wall at r_max < gamma L shields the
//  longitudinal space charge, g = 1 + 2 ln(min(r_max, gamma L)/a) (a beam
//  pipe); r_max >= gamma L gives the free-space field. In an undulator the
//  radiation needs r_max >= 5 w with w^2 = a^2 + (lambda_r Z / (4 pi a))^2 at
//  the end of the undulator (Z: distance from its entrance).
//  (Checked with a gamma = 100 bunch: r_max = gamma L shields the energy
//  spread from the space charge by 1.4 %, 2 gamma L by 0.3 %.)
//  Window: the longitudinal field reaches the distance d beyond the ends of
//  the bunch. In free space E_z(d) = Q/(4 pi eps0 gamma^2 d^2), against
//  E_z(edge) ~ Q/(4 pi eps0 a max(a, gamma L)) (rod or disk in the rest
//  frame); inside the wall it decays as exp(-2.405 gamma d / r_max). 1 %:
//      d = min( 10 sqrt(a max(a, gamma L)) / gamma, ln(100) r_max / (2.405 gamma) ).
//  Inside a plasma the field is screened over c/omega_p: no advice then.
//  The bunch slips back by int (1 + K_eff^2)/(2 gamma^2) dz (K_eff^2 = K^2/2
//  planar, K^2 helical), so the gap behind it is checked at the end of the run.
// ============================================================================
namespace {
std::string fmtUm(double x, const Units& u) {
  char buf[64];
  const double um = x * u.k0inv_um();
  if (um >= 1e4) std::snprintf(buf, sizeof buf, "%.3g mm", um * 1e-3);
  else std::snprintf(buf, sizeof buf, "%.3g um", um);
  return buf;
}
}  // namespace

std::vector<Config::BoxAdvice> Config::boxAdvice() const {
  std::vector<BoxAdvice> out;
  if (plasma.present) return out;  // fields screened over c / omega_p
  for (const auto& B : bunches) {
    if (!(B.gamma < 1e5)) continue;
    BoxAdvice A;
    A.name = B.name;
    const double half = (B.shape == BunchConfig::Shape::Gaussian) ? 3.0 * B.sigma_z : 0.5 * B.length;
    const double L = 2.0 * half, a = std::max(B.sigma_x, B.sigma_y), g2 = 2.0 * B.gamma * B.gamma;
    A.rmax_need = B.gamma * L;
    A.rmax_why = "gamma L";
    for (const auto& U : undulators) {
      const double Z = std::max(0.0, std::min(grid.z_end, U.z_start + U.length) - U.z_start);
      const double lr = U.resonance(B.gamma), d = lr * Z / (4.0 * PI * a);
      const double w5 = 5.0 * std::sqrt(a * a + d * d);
      if (w5 > A.rmax_need) { A.rmax_need = w5; A.rmax_why = "5 x the size of the undulator radiation"; }
    }
    A.rmax_now = grid.r_max;
    A.gap_need = std::min(10.0 * std::sqrt(a * std::max(a, B.gamma * L)) / B.gamma,
                          std::log(100.0) * grid.r_max / (2.405 * B.gamma));
    // slippage over the run
    double slip = grid.z_end / g2;
    for (const auto& U : undulators) {
      const double Z = std::max(0.0, std::min(grid.z_end, U.z_start + U.length) - U.z_start);
      double k2 = 0.0;  // mean K^2 over the part inside the run (taper: K(z))
      for (int j = 0; j < 200; ++j) {
        const double Kz = U.Kat(U.z_start + (j + 0.5) / 200.0 * Z);
        k2 += Kz * Kz / 200.0;
      }
      slip += Z * ((U.type == UndulatorConfig::Type::Planar) ? 0.5 * k2 : k2) / g2;
    }
    A.slip = slip;
    A.front_now = grid.window() - (B.zeta_center + half);
    A.rear_end = B.zeta_center - half - slip;
    // smallest window that keeps both gaps
    A.window_need = A.gap_need + slip + L + A.gap_need;
    out.push_back(A);
  }
  return out;
}

// ============================================================================
Config Config::fromDeck(const Deck& deck, std::vector<std::string>* warnings) {
  std::vector<std::string> warn;
  deck.checkBlockNames({"constants", "grid", "laser", "flying_focus", "plasma", "ion", "bunch", "undulator", "quadrupole",
                        "plasma_lens", "solenoid", "chicane", "output", "numerics"});
  Config c;
  g_consts = &c.constants;  // the expressions of all blocks see the constants
  struct Reset { ~Reset() { g_consts = nullptr; } } reset;
  const Block* grid = deck.single("grid");
  if (!grid) throw DeckError(deck.source + ": the deck needs a 'grid' block");
  c.units.lambda0_um = readLambda0(deck, grid);

  if (const Block* b = deck.single("constants")) readConstants(*b, deck, c.units, c.constants);
  readGrid(*grid, deck, c.units, c.grid, warn);
  for (const Block* b : deck.all("laser")) {
    LaserConfig L;
    readLaser(*b, deck, c.units, c.grid, L, warn);
    c.lasers.push_back(L);
  }
  if (!c.lasers.empty() && c.grid.M < 1)
    throw DeckError(deck.where(grid->line) + ": a laser needs azimuthal mode 1: set grid.m_max >= 1");
  if (const Block* b = deck.single("flying_focus")) readFlyingFocus(*b, deck, c.units, c.grid, c.flying_focus, warn);
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
  std::vector<std::string> planarStepWarnings;
  for (const Block* b : deck.all("undulator")) {
    UndulatorConfig U;
    readUndulator(*b, deck, c.units, U);
    // planar undulator, high current (benchmarks/harmonics: 4 kA, gamma 1000, K = 1.5), beam pushed in the
    // fields of s^0 (numerics.beam_field = level): with the field step lambda_u / 16 the beam blew up
    // transversely (eps_n 0 -> 2.8 um in 25 periods; lambda_u / 24: 0.13 um, / 32: 0.04 um, / 64: 0.02 um);
    // the helical case is far less sensitive. With beam_field = mean (default with undulators): 0.04 um at / 16.
    const double dzMax = c.grid.zTable ? c.grid.dz_max : c.grid.dz;
    double qmax = 0.0;
    for (const auto& B : c.bunches) qmax = std::max(qmax, B.charge_pC);
    if (U.type == UndulatorConfig::Type::Planar && qmax >= 1.0 && dzMax > U.period / 32.0 * 1.0001)
      planarStepWarnings.push_back(b->label() + ": planar undulator with dz = period / " +
                                   std::to_string(int(std::lround(U.period / dzMax))) +
                                   " and numerics.beam_field = level: use dz <= period / 32 with intense beams (at "
                                   "period / 16 a 4 kA beam blew up transversely), or beam_field = mean");
    c.undulators.push_back(U);
  }
  for (const char* kind : {"quadrupole", "plasma_lens", "solenoid"})
    for (const Block* b : deck.all(kind)) {
      MagnetConfig G;
      readMagnet(*b, deck, c.units, G, warn);
      c.magnets.push_back(G);
    }
  for (const Block* b : deck.all("chicane")) {
    ChicaneConfig C;
    readChicane(*b, deck, c.units, C);
    c.chicanes.push_back(C);
  }
  std::sort(c.chicanes.begin(), c.chicanes.end(), [](const ChicaneConfig& a, const ChicaneConfig& b) { return a.z < b.z; });
  if (const Block* b = deck.single("output")) readOutput(*b, deck, c.units, c.output);
  // FEL diagnostics: on with undulators unless switched off; wavelength of the bunching
  {
    auto& O = c.output;
    O.fel = (O.fel_mode == 1) || (O.fel_mode == -1 && !c.undulators.empty() && !c.bunches.empty());
    if (O.fel && !(O.fel_wavelength > 0)) {
      if (c.undulators.empty() || c.bunches.empty())
        throw DeckError(deck.source + ": output.fel needs output.fel_wavelength (no undulator and bunch to take the "
                        "resonant wavelength from)");
      O.fel_wavelength = c.undulators[0].resonance(c.bunches[0].gamma);
    }
    if (!(O.fel_slice > 0)) O.fel_slice = 20.0 * O.fel_wavelength;
  }
  if (c.output.lasy && !(c.output.lasy_omega > 0)) c.output.lasy_omega = c.lasers.empty() ? 1.0 : c.lasers[0].omega;
  if (c.output.lasy && c.grid.M < 1)
    throw DeckError(deck.source + ": output.lasy needs the laser modes: set grid.m_max >= 1");
  if (const Block* b = deck.single("numerics")) readNumerics(*b, deck, c.units, c.numerics);
  deck.checkAllUsed();
  // box size for bunches with v < c (Config::boxAdvice): warn if too small
  for (const auto& A : c.boxAdvice()) {
    if (A.rmax_now < A.rmax_need)
      warn.push_back("bunch '" + A.name + "': the radial boundary is closer than the reach of its fields: r_max = " +
                     fmtUm(A.rmax_now, c.units) + " < " + fmtUm(A.rmax_need, c.units) + " (" + A.rmax_why +
                     "); a wall closer than gamma L shields the space charge like a beam pipe");
    if (A.front_now < A.gap_need || A.rear_end < A.gap_need)
      warn.push_back("bunch '" + A.name + "': its space-charge field reaches the ends of the window (gap ahead " +
                     fmtUm(A.front_now, c.units) + ", behind at the end " + fmtUm(A.rear_end, c.units) +
                     "; needed " + fmtUm(A.gap_need, c.units) + "): E_z = 0 at the front acts as a uniform E_z");
  }
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
  // transport elements: the sub-steps of the beam push must resolve the element and the betatron phase
  for (const auto& G : c.magnets) {
    if (G.length < 4.0 * c.numerics.beam_dz_max)
      warn.push_back(G.label + ": shorter than 4 sub-steps of the beam push (numerics.beam_dz_max)");
    for (const auto& B : c.bunches) {
      const double k = (G.kind == MagnetConfig::Kind::Solenoid) ? 0.5 * std::fabs(G.field) / B.gamma
                                                                  : std::sqrt(std::fabs(G.gradient) / B.gamma);
      if (k * c.numerics.beam_dz_max > 0.1)
        warn.push_back(G.label + ", bunch '" + B.name + "': numerics.beam_dz_max exceeds 0.1 / k_beta (" +
                       fmtUm(0.1 / k, c.units) + ")");
    }
  }

  // beam push in the mean field of s^0 and s^1: default with undulators (the beam's self-field turns with
  // the wiggle; pushed in the fields of s^0 only, the beam blows up for dz > lambda_u/32)
  c.numerics.beam_field_auto = c.numerics.beam_field < 0;
  if (c.numerics.beam_field < 0) c.numerics.beam_field = c.undulators.empty() ? BEAM_FIELD_LEVEL : BEAM_FIELD_MEAN;
  if (c.numerics.beam_field == BEAM_FIELD_LEVEL) warn.insert(warn.end(), planarStepWarnings.begin(), planarStepWarnings.end());
  if (c.lasers.empty() && c.bunches.empty() && !c.flying_focus.present)
    warn.push_back("the deck defines neither a laser nor a bunch");
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
    if (L.file) {
      const auto& D = *L.file;
      const double tfs = 1.0 / units.omega0_per_fs();
      line("laser[%zu]: from file %s (iteration %ld, mesh %s, %s%s)", i, D.path.c_str(), D.iteration, D.mesh.c_str(),
           D.cartesian ? "cartesian" : "thetaMode", D.software.empty() ? "" : (", written by " + D.software).c_str());
      line("        lambda_L = %.4g um, peak |E| = %.4g m c omega0/e (a0 = %.4g%s), polarization (%.3g%+.3gi, "
           "%.3g%+.3gi), %s", units.lambda0_um / L.omega, L.a0 * L.omega, L.a0, L.file_scale != 1.0 ? ", rescaled" : "",
           D.px.real(), D.px.imag(), D.py.real(), D.py.imag(), D.field.c_str());
      line("        t = %.4g ... %.4g fs (%d points), %s; t_center = %.4g fs at zeta_c = %.4g um", D.t0 * tfs,
           (D.t0 + (D.nt - 1) * D.dt) * tfs, D.nt,
           D.cartesian ? (std::to_string(D.nx) + " x " + std::to_string(D.ny) + " points in x, y").c_str()
                       : ("r <= " + std::to_string((D.r0 + (D.nr - 1) * D.dr) * um).substr(0, 7) + " um (" +
                          std::to_string(D.nr) + " points), " + std::to_string((D.ncomp + 1) / 2) + " modes").c_str(),
           L.t_center * tfs, L.zeta_center * um);
      line("        plane of the file z = %.4g mm; %s%s", L.file_z * um * 1e-3,
           L.boundary ? "fed through the absorber layer at every step (inject = boundary)" : "initial field (inject = box)",
           L.background > 0 ? (", propagation in n = " + std::to_string(L.background * units.nc_cm3()).substr(0, 9) +
                               " cm^-3").c_str() : "");
      continue;
    }
    // peak power, linear polarization: P = 21.5 GW (a0 w0 / lambda)^2; circular: x2
    const double lam_um = units.lambda0_um / L.omega;
    double P_TW = 21.5e-3 * std::pow(L.a0 * L.waist * um / lam_um, 2.0);
    if (L.pol == LaserConfig::Pol::Circular) P_TW *= 2.0;
    line("laser[%zu]: a0 = %.4g, lambda = %.4g um, w0 = %.4g um (zR = %.4g mm), tau_FWHM = %.4g fs, "
         "zeta_c = %.4g um, focus at z = %.4g mm, %s, P = %.3g TW",
         i, L.a0, lam_um, L.waist * um, L.zR() * um * 1e-3, L.length_1e * std::sqrt(2.0 * std::log(2.0)) * fs,
         L.zeta_center * um, L.z_focus * um * 1e-3, L.pol == LaserConfig::Pol::Linear ? "linear" : "circular",
         P_TW);
    if (!L.profile_fn.empty())
      line("        profile_function(r) = %s (field at z = 0; w0, zR, P: estimates from its 1/e radius)",
           L.profile_fn.text().c_str());
    if (!L.envelope_fn.empty())
      line("        envelope_function(zeta) = %s (tau_FWHM: estimate)", L.envelope_fn.text().c_str());
  }
  for (const auto& A : boxAdvice())
    line("box for bunch '%s' (v < c): r_max >= %s (%s; now %s), gaps >= %s at both ends (now %s ahead, %s behind "
         "after %s slippage): window >= %s (now %s)",
         A.name.c_str(), fmtUm(A.rmax_need, units).c_str(), A.rmax_why.c_str(), fmtUm(A.rmax_now, units).c_str(),
         fmtUm(A.gap_need, units).c_str(), fmtUm(A.front_now, units).c_str(), fmtUm(A.rear_end, units).c_str(),
         fmtUm(A.slip, units).c_str(), fmtUm(A.window_need, units).c_str(), fmtUm(grid.window(), units).c_str());
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
    if (!U.K_fn.empty()) line("        taper K_function(z) = %s (K above: at the entrance)", U.K_fn.text().c_str());
    if (!res.empty()) line("        resonant wavelength (on axis): %s", res.c_str());
  }
  for (const auto& G : magnets) {
    if (G.kind == MagnetConfig::Kind::Solenoid)
      line("%s: B_z = %.4g T, z = %.4g ... %.4g mm, fringes %.4g mm", G.label.c_str(), G.field * units.B0_T(),
           G.z_start * um * 1e-3, (G.z_start + G.length) * um * 1e-3, G.ramp * um * 1e-3);
    else
      line("%s: gradient %.4g T/m (%s electrons), z = %.4g ... %.4g mm%s", G.label.c_str(),
           G.gradient * units.B0_T() / (um * 1e-6),
           G.kind == MagnetConfig::Kind::PlasmaLens ? (G.gradient > 0 ? "focuses" : "defocuses")
                                                    : (G.gradient > 0 ? "focuses in x, defocuses in y" : "focuses in y, defocuses in x"),
           G.z_start * um * 1e-3, (G.z_start + G.length) * um * 1e-3,
           G.radius < 1e299 ? (", radius " + std::to_string(G.radius * um) + " um").c_str() : "");
  }
  for (const auto& C : chicanes)
    line("chicane at z = %.4g mm: R56 = %.4g um, delay %.4g um, reference energy %.5g MeV", C.z * um * 1e-3, C.r56 * um,
         C.delay * um, (C.gamma_ref - 1.0) * MC2_MEV);
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
    if (!plasma.fn.empty()) line("        times density_function(x, y, z, r) = %s", plasma.fn.text().c_str());
    line("        macroparticles per cell: %d x %d, and %d x %d for r < %.4g um", plasma.ppc_r, plasma.ppc_phi,
         plasma.ppc_r_axis, plasma.ppc_phi_axis, plasma.r_axis_refine * um);
  }
  if (plasma.present) {
    if (!plasma.free_electrons) line("%s", "        no free electrons at the start (neutral gas)");
    if (plasma.smooth_length > 0.0 || plasma.smooth_length_m0 > 0.0) {
      const double am0 = plasma.smooth_length_m0 >= 0.0 ? plasma.smooth_length_m0 : plasma.smooth_length;
      line("        radial smoothing of the plasma sources over a = %.4g um (%.4g k0^-1, %.3g c/omega_p) for m >= 1, "
           "%.4g um (%.4g k0^-1) for m = 0", plasma.smooth_length * units.k0inv_um(), plasma.smooth_length,
           plasma.smooth_length * std::sqrt(plasma.n0), am0 * units.k0inv_um(), am0);
    }
    double fsum = 0.0;
    bool anyIonization = false;
    for (size_t i = 0; i < ions.size(); ++i) {
      const auto& I = ions[i];
      fsum += I.fraction();
      line("ion[%zu] '%s': Z = %d, mass = %.5g u (%.6g m_e), atom density %.4g cm^-3 (%.4g n0), %s", i,
           I.name.c_str(), I.Z, I.mass / U_OVER_ME, I.mass, I.atoms * plasma.n0 * units.nc_cm3(), I.atoms,
           I.mobile ? "mobile" : "fixed");
      if (I.own)
        line("        own profile%s%s%s (atom density above: reference; with %d electrons per atom)",
             I.prof.z_points.empty() ? "" : ", longitudinal table",
             I.prof.radial == PlasmaConfig::Radial::Uniform ? "" : ", radial profile",
             I.prof.fn.empty() ? "" : (", density_function " + I.prof.fn.text()).c_str(), I.Z);
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
         B.shape == BunchConfig::Shape::Gaussian || B.shape == BunchConfig::Shape::File ? "sigma_z" : "length",
         (B.shape == BunchConfig::Shape::Gaussian || B.shape == BunchConfig::Shape::File ? B.sigma_z : B.length) * um,
         (B.shape == BunchConfig::Shape::File ? B.file_center : B.zeta_center) * um,
         B.macroparticles);
    if (B.shape == BunchConfig::Shape::File)
      line("        from %s (rms length %.4g um, extent %.4g ... %.4g um)", B.file_info.c_str(), B.sigma_z * um,
           B.zeta_min * um, B.zeta_max * um);
    if (B.fawley)
      line("        shot noise: Fawley beamlets of %d macroparticles at %.4g um", B.beamlet, B.noise_wavelength * um);
    if (B.axis_refine > 1)
      line("        axis refinement: %d times the macroparticles at r < %s, with 1/%d of the charge", B.axis_refine,
           B.axis_radius > 0 ? (std::to_string(B.axis_radius * um) + " um").c_str() : "4 radial cells",
           B.axis_refine);
    if (B.shape == BunchConfig::Shape::Function)
      line("        profile_z(zeta) = %s on %.4g ... %.4g um", B.profile_z.text().c_str(), B.zeta_min * um,
           B.zeta_max * um);
    if (!B.gamma_fn.empty())
      line("        gamma_function(zeta) = %s (gamma above: mean over the profile)", B.gamma_fn.text().c_str());
    if (!B.profile_r.empty())
      line("        profile_r(x, y, r, zeta) = %s in r < %.4g um (sigma_x,y above: rms of the profile)",
           B.profile_r.text().c_str(), B.r_extent * um);
    if (!B.ux_fn.empty() || !B.uy_fn.empty())
      line("        mean momenta: u_x = %s, u_y = %s (energy kept)", B.ux_fn.empty() ? "0" : B.ux_fn.text().c_str(),
           B.uy_fn.empty() ? "0" : B.uy_fn.text().c_str());
  }
  if (!constants.empty()) {
    std::string s;
    for (const auto& kv : constants) {
      std::snprintf(buf, sizeof buf, "%s%s = %.6g", s.empty() ? "" : ", ", kv.first.c_str(), kv.second);
      s += buf;
    }
    line("constants (code units): %s", s.c_str());
  }
  if (output.hdf5 || output.openpmd) {
    std::string f = output.text ? "text" : "";
    if (output.hdf5) f += f.empty() ? "hdf5" : " + hdf5";
    if (output.openpmd) f += f.empty() ? "openpmd" : " + openpmd";
    line("output format: %s (HDF5 %s precision, compression %d)", f.c_str(), output.hdf5_single ? "single" : "double",
         output.hdf5_compression);
  }
  if (flying_focus.present) {
    const auto& F = flying_focus;
    line("flying focus: luminal Bessel pulse, a0 = %.4g on the axis, cone angle %.4g mrad (core radius %.4g um, "
         "transverse wavelength %.4g um), length_1e = %.4g um at zeta = %.4g um; fed through the absorber layer r > %.4g um",
         F.a0, F.theta * 1e3, 2.405 / std::sin(F.theta) * um, 2.0 * PI / std::sin(F.theta) * um, F.length_1e * um,
         F.zeta_center * um, grid.r_absorber * um);
    if (F.background > 0) line("        incident field in a uniform plasma of %.4g n_c", F.background);
  }
  if (output.fel)
    line("FEL diagnostics: fel.txt every step, bunching at %.5g um", output.fel_wavelength * um);
  if (output.fel && output.fel_harmonics > 1)
      line("        fel_harmonics.txt: harmonics 1 ... %d, modes m = 0 ... %d, slices of %.4g um", output.fel_harmonics,
           grid.M, output.fel_slice * um);
  if (output.openpmd) {
    const double dr = output.openpmd_dr > 0 ? output.openpmd_dr : grid.dr_min;
    double dzt = output.openpmd_dzeta;
    if (!(dzt > 0)) {
      dzt = 1e300;
      for (int n = 0; n < grid.Nz; ++n) dzt = std::min(dzt, grid.zetaNode(n + 1) - grid.zetaNode(n));
    }
    const long nr = long(std::floor(grid.r_max / dr + 1e-9)) + 1, nz = long(std::floor(grid.window() / dzt + 1e-9)) + 1;
    if (output.openpmd_uniform)
      line("openPMD: %s/openpmd/openpmd_%%T.h5, fields on a uniform grid dr = %.4g um, dzeta = %.4g um (%ld x %ld "
           "points, %d modes)%s",
           output.dir.c_str(), dr * um, dzt * um, nr, nz, 2 * grid.M + 1,
           output.openpmd_native ? ", and on the native grid (records native_*)" : "");
    else
      line("openPMD: %s/openpmd/openpmd_%%T.h5, fields on the native grid (records native_*, %d modes)",
           output.dir.c_str(), 2 * grid.M + 1);
  }
  if (output.checkpoint_every > 0 || output.walltime_hours > 0)
    line("checkpoints: every %d steps (keep %d)%s", output.checkpoint_every, output.checkpoint_keep,
         output.walltime_hours > 0 ? (", stop at " + std::to_string(output.walltime_hours) + " h wall time").c_str()
                                   : "");
  if (numerics.trap_gamma > 0)
    line("trapping: a streaming particle goes to the beam list only with gamma >= %.4g (numerics.trap_gamma)",
         numerics.trap_gamma);
  if (numerics.theta_stiff > 0 || numerics.beam_consistent)
    line("beam self-fields: theta per radial mode %s, consistent beam sources %s",
         numerics.theta_stiff > 0 ? ("on (theta_stiff = " + std::to_string(numerics.theta_stiff) + ")").c_str() : "off",
         numerics.beam_consistent ? "on" : "off");
  line("output: dir '%s', every %d steps; numerics: %d iterations/slice, theta %.3g, trap advance %.3g dz, beam "
       "filter %d, beam pusher %s, beam field %s, seed %lu, backend %s",
       output.dir.c_str(), output.every, numerics.iterations, numerics.theta, numerics.trap_advance, numerics.beam_filter,
       numerics.beam_vay ? "vay" : "boris", numerics.beam_field == BEAM_FIELD_MEAN ? "mean of s^0, s^1" : "s^0",
       numerics.seed,
       numerics.backend == NumericsConfig::Backend::Cuda     ? "cuda"
       : numerics.backend == NumericsConfig::Backend::OpenMP ? "openmp"
                                                             : "serial");
  return os.str();
}

}  // namespace gem
