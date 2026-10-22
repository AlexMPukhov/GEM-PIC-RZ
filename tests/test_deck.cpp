// ============================================================================
//  test_deck.cpp -- unit tests of the deck parser, units and configuration
//  (no external framework; exit code = number of failed checks)
// ============================================================================
#include "test_deck_common.hpp"

int failures = 0, checks = 0;

int main() {
  const Units u;  // lambda0 = 0.8 um
  const double k0 = u.k0inv_um();

  test_deck_basic(u, k0);
  test_deck_grids(u, k0);
  test_deck_physics(u, k0);
  test_deck_beams(u, k0);

  std::printf("%d checks, %d failed\n", checks, failures);
  return failures;
}
