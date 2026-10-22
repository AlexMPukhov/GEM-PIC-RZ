// ============================================================================
//  config_internal.hpp -- shared part of the deck readers (src/config*.cpp): BlockReader (typed access to the
//  keys of one block, with units), the deck constants, helpers, and the block readers
// ============================================================================
#pragma once
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
namespace cfgdetail {

// ============================================================================
//  BlockReader: typed access to the keys of one block, with units
// ============================================================================

// constants of the deck (block 'constants'), visible to every expression
inline const std::map<std::string, double>* g_consts = nullptr;
inline const std::map<std::string, double> NO_CONSTS;

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
inline double logcosh(double x) {
  const double a = std::fabs(x);
  return a + std::log1p(std::exp(-2.0 * a)) - std::log(2.0);
}

inline double interpTable(const std::vector<double>& xs, const std::vector<double>& ys, double x) {
  if (x <= xs.front()) return ys.front();
  if (x >= xs.back()) return ys.back();
  const size_t i = size_t(std::upper_bound(xs.begin(), xs.end(), x) - xs.begin()) - 1;
  return ys[i] + (ys[i + 1] - ys[i]) * ((x - xs[i]) / (xs[i + 1] - xs[i]));
}

// a length in um or mm, for messages
inline std::string fmtUm(double x, const Units& u) {
  char buf[64];
  const double um = x * u.k0inv_um();
  if (um >= 1e4) std::snprintf(buf, sizeof buf, "%.3g mm", um * 1e-3);
  else std::snprintf(buf, sizeof buf, "%.3g um", um);
  return buf;
}

// block readers (src/config_read_*.cpp)
double readLambda0(const Deck& deck, const Block* grid);
void readGrid(const Block& b, const Deck& d, const Units& u, GridConfig& g, std::vector<std::string>& warn);
void readConstants(const Block& b, const Deck& d, const Units& u, std::map<std::string, double>& C);
void readLaser(const Block& b, const Deck& d, const Units& u, const GridConfig& g, LaserConfig& L, std::vector<std::string>& warn);
void readFlyingFocus(const Block& b, const Deck& d, const Units& u, const GridConfig& g, FlyingFocusConfig& F, std::vector<std::string>& warn);
void readProfile(const BlockReader& r, const GridConfig& g, PlasmaConfig& P);
void readPlasma(const Block& b, const Deck& d, const Units& u, const GridConfig& g, PlasmaConfig& P, std::vector<std::string>& warn);
void readIon(const Block& b, const Deck& d, const Units& u, const PlasmaConfig& P, IonConfig& I, std::vector<std::string>& warn);
void readShotNoise(const BlockReader& r, BunchConfig& B);
void readBunchFile(const Block& b, const Deck& d, const Units& u, const GridConfig& g, BunchConfig& B, std::vector<std::string>& warn);
void readBunch(const Block& b, const Deck& d, const Units& u, const GridConfig& g, BunchConfig& B, std::vector<std::string>& warn);
void readOutput(const Block& b, const Deck& d, const Units& u, OutputConfig& O);
void readUndulator(const Block& b, const Deck& d, const Units& u, UndulatorConfig& U);
void readMagnet(const Block& b, const Deck& d, const Units& u, MagnetConfig& G, std::vector<std::string>& warn);
void readChicane(const Block& b, const Deck& d, const Units& u, ChicaneConfig& C);
void readNumerics(const Block& b, const Deck& d, const Units& u, NumericsConfig& N);

}  // namespace cfgdetail
}  // namespace gem
