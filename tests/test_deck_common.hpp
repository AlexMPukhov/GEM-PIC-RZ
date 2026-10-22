// ============================================================================
//  test_deck_common.hpp -- checks, error expectations and the base deck of the deck tests (tests/test_deck*.cpp)
// ============================================================================
#pragma once
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <functional>
#include <string>

#include "atomic_data.hpp"
#include "config.hpp"
#include "deck.hpp"
#include "grid.hpp"
#include "ionization.hpp"
#include "mapping.hpp"
#include "units.hpp"
#include "zgrid.hpp"
#include "smoothing.hpp"
#include "solver.hpp"
#include "external.hpp"
#include "eigen.hpp"
#include "expr.hpp"

using namespace gem;

extern int failures, checks;

#define CHECK(cond)                                                          \
  do {                                                                       \
    ++checks;                                                                \
    if (!(cond)) {                                                           \
      ++failures;                                                            \
      std::printf("FAILED line %d: %s\n", __LINE__, #cond);                  \
    }                                                                        \
  } while (0)

#define CHECK_NEAR(a, b, rel)                                                                   \
  do {                                                                                          \
    ++checks;                                                                                   \
    const double a_ = (a), b_ = (b);                                                            \
    if (!(std::fabs(a_ - b_) <= (rel) * std::max(std::fabs(a_), std::fabs(b_)) + 1e-300)) {     \
      ++failures;                                                                               \
      std::printf("FAILED line %d: %s = %.12g, expected %.12g\n", __LINE__, #a, a_, b_);         \
    }                                                                                           \
  } while (0)

// the code must throw DeckError whose message contains `part`
inline void expectError(int line, const std::function<void()>& f, const std::string& part) {
  ++checks;
  try {
    f();
  } catch (const DeckError& e) {
    if (std::string(e.what()).find(part) == std::string::npos) {
      ++failures;
      std::printf("FAILED line %d: error message\n   '%s'\n   does not contain '%s'\n", line, e.what(), part.c_str());
    }
    return;
  } catch (const std::exception& e) {
    ++failures;
    std::printf("FAILED line %d: wrong exception type: %s\n", line, e.what());
    return;
  }
  ++failures;
  std::printf("FAILED line %d: no error, expected one containing '%s'\n", line, part.c_str());
}
#define EXPECT_ERROR(code, part) expectError(__LINE__, [&] { code; }, part)

// a small valid deck used by many tests (10 lines)
inline const char* BASE = R"(
grid {
  lambda0 = 0.8 um
  radial_grid = uniform
  dr = 0.5 um
  nr = 40
  dzeta = 0.05 um
  window_length = 40 um
  dz = 10 um
  z_end = 1 mm
}
)";

inline Config cfgOf(const std::string& text) { return Config::fromDeck(Deck::parseString(text, "test.deck")); }

void test_deck_basic(const Units& u, double k0);
void test_deck_grids(const Units& u, double k0);
void test_deck_physics(const Units& u, double k0);
void test_deck_beams(const Units& u, double k0);
