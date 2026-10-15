// ============================================================================
//  test_deck.cpp -- unit tests of the deck parser, units and configuration
//  (no external framework; exit code = number of failed checks)
// ============================================================================
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

static int failures = 0, checks = 0;

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
static void expectError(int line, const std::function<void()>& f, const std::string& part) {
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
static const char* BASE = R"(
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

static Config cfgOf(const std::string& text) { return Config::fromDeck(Deck::parseString(text, "test.deck")); }

int main() {
  const Units u;  // lambda0 = 0.8 um
  const double k0 = u.k0inv_um();

  // ------------------------------------------------------------ syntax
  {
    Deck d = Deck::parseString(
        "# comment\n"
        "Grid {   ! other comment\n"
        "  A = 1, 2  3 um   # list with a unit\n"
        "  s = \"two words\"; t = 5nm\n"
        "}\n"
        "laser { a0 = 2 }\n"
        "laser\n{\n a0 = 3\n}\n",
        "syntax.deck");
    CHECK(d.blocks.size() == 3);
    CHECK(d.blocks[0].name == "grid");          // case-insensitive names
    CHECK(d.blocks[0].find("a") != nullptr);    // case-insensitive keys
    CHECK(d.blocks[0].find("a")->tokens.size() == 4);
    CHECK(d.blocks[0].find("s")->tokens[0] == "two words");
    CHECK(d.blocks[0].find("t")->tokens[0] == "5nm");
    CHECK(d.all("laser").size() == 2);
    CHECK(d.blocks[1].label() == "laser[0]" && d.blocks[2].label() == "laser[1]");
    CHECK(d.blocks[0].label() == "grid");
    CHECK(d.blocks[2].find("a0")->line == 9);
  }
  EXPECT_ERROR(Deck::parseString("grid { a = 1\n", "x.deck"), "x.deck:1: block 'grid' is not closed");
  EXPECT_ERROR(Deck::parseString("grid { a = 1\n a = 2 }", "x.deck"), "x.deck:2: key 'a' given twice");
  EXPECT_ERROR(Deck::parseString("grid { a = }", "x.deck"), "missing value for key 'a'");
  EXPECT_ERROR(Deck::parseString("grid { a = 1 b = 2 }", "x.deck"), "unexpected '='");
  EXPECT_ERROR(Deck::parseString("grid a = 1", "x.deck"), "expected '{' after block name 'grid'");
  EXPECT_ERROR(Deck::parseString("grid { a = \"open }", "x.deck"), "unterminated string");
  EXPECT_ERROR(Deck::parseString("grid { a = 1\nlaser { b = 2 }", "x.deck"), "block 'grid' (line 1) is not closed with '}' before block 'laser'");

  // ------------------------------------------------------------ units
  CHECK_NEAR(toCode(1.0, "um", Dim::Length, u), 1.0 / k0, 1e-14);
  CHECK_NEAR(toCode(1.0, "mm", Dim::Length, u), 1e3 / k0, 1e-14);
  CHECK_NEAR(toCode(800.0, "nm", Dim::Length, u), 2.0 * PI, 1e-14);
  CHECK_NEAR(toCode(1.0, "lambda0", Dim::Length, u), 2.0 * PI, 1e-14);
  CHECK(toCode(3.0, "k0^-1", Dim::Length, u) == 3.0);
  CHECK_NEAR(toCode(1.0, "T0", Dim::Time, u), 2.0 * PI, 1e-14);
  CHECK_NEAR(toCode(1.0, "fs", Dim::Time, u), 2.0 * PI * 0.299792458 / 0.8, 1e-12);  // omega0 * 1 fs
  CHECK_NEAR(toCode(1.742e21, "cm^-3", Dim::Density, u), 1.742e21 / u.nc_cm3(), 1e-14);
  CHECK_NEAR(toCode(1e24, "m^-3", Dim::Density, u), 1e18 / u.nc_cm3(), 1e-14);
  CHECK_NEAR(u.codeToPC(toCode(5.0, "pC", Dim::Charge, u)), 5.0, 1e-13);
  CHECK_NEAR(u.codeToPC(toCode(1.0, "nC", Dim::Charge, u)), 1000.0, 1e-13);
  CHECK_NEAR(toCode(2.0, "GeV", Dim::Energy, u), 2000.0, 1e-14);
  CHECK_NEAR(toCode(90.0, "deg", Dim::Angle, u), PI / 2, 1e-14);
  CHECK_NEAR(toCode(0.1, "%", Dim::None, u), 1e-3, 1e-14);
  CHECK_NEAR(toCode(0.3, "mm*mrad", Dim::Length, u), 0.3 / k0, 1e-14);
  bool threw = false;
  try { toCode(1.0, "", Dim::Length, u); } catch (const std::invalid_argument&) { threw = true; }
  CHECK(threw);
  threw = false;
  try { toCode(1.0, "fs", Dim::Length, u); } catch (const std::invalid_argument&) { threw = true; }
  CHECK(threw);

  // ------------------------------------------------------------ minimal config
  {
    std::vector<std::string> warn;
    Config c = Config::fromDeck(Deck::parseString(BASE, "t"), &warn);
    CHECK(c.grid.Nr == 40);
    CHECK_NEAR(c.grid.r_max, 20.0 / k0, 1e-12);
    CHECK(c.grid.Nz == 800);
    CHECK_NEAR(c.grid.dzeta, 0.05 / k0, 1e-12);
    CHECK_NEAR(c.grid.z_end, 1000.0 / k0, 1e-12);
    CHECK(!c.plasma.present && c.lasers.empty() && c.bunches.empty());
    CHECK(warn.size() == 2);  // no laser/bunch, vacuum
  }

  // ------------------------------------------------------------ stretched grid
  {
    Config c = cfgOf(R"(grid {
      lambda0 = 0.8 um;  dr_min = 40 nm;  dr_max = 0.4 um
      r_transition = 2.5 um;  transition_cells = 6;  r_max = 24 um
      dzeta = 30 nm;  nzeta = 100;  dz = 5 um;  z_end = 50 um })");
    const auto& g = c.grid;
    CHECK(g.stretched);
    CHECK_NEAR(g.mapping(g.xi0), 2.5 / k0, 1e-9);  // transition radius hit
    CHECK(g.mapping(g.Nr) >= 24.0 / k0 * (1 - 1e-12));
    CHECK(g.mapping(g.Nr - 1) < 24.0 / k0);
    CHECK_NEAR(g.mapping(1) - g.mapping(0), 0.04 / k0, 0.02);  // step near the axis
    CHECK_NEAR(g.mapping(g.Nr) - g.mapping(g.Nr - 1), 0.4 / k0, 0.02);
  }

  // ------------------------------------------------------------ lasers, plasma, bunches
  {
    Config c = cfgOf(std::string(BASE) + R"(
      laser { a0 = 2; waist = 8 um; tau_fwhm = 30 fs; zeta_center = 30 um; z_focus = 1 mm }
      laser { a0 = 0.5; wavelength = 400 nm; waist = 5 um; length_1e = 2 um; zeta_center = 20 um
              polarization = circular; helicity = -1 }
      plasma { density = 1e18 cm^-3
               z_points = 0, 100, 900 um;  z_values = 0, 1, 0.5
               radial_profile = parabolic; channel_depth = 0.5; r_channel = 10 um; r_max = 15 um }
      bunch { name = driver; charge = 1 nC; energy = 1 GeV; sigma_x = 5 um; sigma_z = 10 um
              zeta_center = 20 um; macroparticles = 100 }
      bunch { species = positron; charge = 10 pC; gamma = 200; sigma_x = 1 um; sigma_y = 2 um
              distribution_z = flattop; length = 3 um; zeta_center = 10 um; energy_spread = 1 % })");
    CHECK(c.lasers.size() == 2 && c.bunches.size() == 2);
    CHECK_NEAR(c.lasers[0].length_1e * std::sqrt(2 * std::log(2.0)), toCode(30, "fs", Dim::Time, c.units), 1e-12);
    CHECK_NEAR(c.lasers[0].z_focus, 1000.0 / k0, 1e-12);
    CHECK_NEAR(c.lasers[1].omega, 2.0, 1e-12);
    CHECK(c.lasers[1].pol == LaserConfig::Pol::Circular && c.lasers[1].helicity == -1);
    const auto& P = c.plasma;
    CHECK(P.present);
    CHECK_NEAR(P.n0, 1e18 / c.units.nc_cm3(), 1e-12);
    CHECK_NEAR(P.longitudinal(50.0 / k0), 0.5, 1e-12);
    CHECK_NEAR(P.longitudinal(500.0 / k0), 0.75, 1e-12);
    CHECK_NEAR(P.longitudinal(-1.0), 0.0, 1e-12);
    CHECK_NEAR(P.longitudinal(1e9), 0.5, 1e-12);
    CHECK_NEAR(P.radialFactor(10.0 / k0), 1.5, 1e-12);
    CHECK(P.density(16.0 / k0, 200.0 / k0) == 0.0);
    CHECK(c.bunches[0].name == "driver" && c.bunches[0].qm < 0 && c.bunches[0].charge < 0);
    CHECK_NEAR(c.bunches[0].charge_pC, 1000.0, 1e-12);
    CHECK_NEAR(c.bunches[0].gamma, 1.0 + 1000.0 / MC2_MEV, 1e-12);
    CHECK(c.bunches[1].name == "bunch1" && c.bunches[1].qm > 0 && c.bunches[1].charge > 0);
    CHECK(c.bunches[1].shape == BunchConfig::Shape::Flattop);
    CHECK_NEAR(c.bunches[1].energy_spread, 0.01, 1e-12);
    CHECK_NEAR(c.bunches[1].sigma_y, 2.0 / k0, 1e-12);
  }

  // ------------------------------------------------------------ errors with line numbers
  EXPECT_ERROR(cfgOf(std::string(BASE) + "laser { a0 = 1\n waist = 8 um\n tau_fwhm = 30 fs\n zeta_centre = 5 um }"),
               "needs key 'zeta_center' (found 'zeta_centre' in line 15: misspelled?)");
  EXPECT_ERROR(cfgOf(std::string(BASE) + "laser { a0 = 1\n waist = 8 um\n tau_fwhm = 30 fs\n zeta_center = 5 um\n"
                                         " z_focs = 1 mm }"),
               "test.deck:16: unknown key 'z_focs' in block 'laser'");
  EXPECT_ERROR(cfgOf(std::string(BASE) + "laser { a0 = 1; waist = 8; tau_fwhm = 30 fs; zeta_center = 5 um }"),
               "laser.waist: a length needs a unit");
  EXPECT_ERROR(cfgOf(std::string(BASE) + "laser { a0 = 1; waist = 8 fs; tau_fwhm = 30 fs; zeta_center = 5 um }"),
               "unit 'fs' is a time, but a length is expected");
  EXPECT_ERROR(cfgOf(std::string(BASE) + "laser { a0 = 1; waist = 8 microns; tau_fwhm = 30 fs; zeta_center = 5 um }"),
               "unknown unit 'microns'");
  EXPECT_ERROR(cfgOf(std::string(BASE) + "laser { a0 = 1; waist = 8 um; zeta_center = 5 um }"),
               "needs 'tau_fwhm' or 'length_1e'");
  EXPECT_ERROR(cfgOf(std::string(BASE) + "laser { a0 = 1; waist = 8 um; tau_fwhm = 3 fs; length_1e = 1 um; "
                                         "zeta_center = 5 um }"),
               "give either 'tau_fwhm' or 'length_1e'");
  EXPECT_ERROR(cfgOf(std::string(BASE) + "laser { a0 = 1; waist = 8 um; tau_fwhm = 3 fs; zeta_center = 50 um }"),
               "outside the window");
  EXPECT_ERROR(cfgOf(std::string(BASE) + "\nplasma { density = 1e18 cm^-3 }\nplasma { density = 1e18 cm^-3 }"),
               "block 'plasma' may appear only once");
  EXPECT_ERROR(cfgOf(std::string(BASE) + "\nlasers { a0 = 1 }"), "unknown block 'lasers'");
  EXPECT_ERROR(cfgOf("laser { a0 = 1 }"), "needs a 'grid' block");
  EXPECT_ERROR(cfgOf(std::string(BASE) + "plasma { density = 1e18 cm^-3; radial_profile = gauss }"),
               "'gauss' is not allowed (choose one of: uniform, parabolic, table)");
  EXPECT_ERROR(cfgOf(std::string(BASE) + "plasma { density = 1e18 cm^-3; z_points = 0 1 mm; z_values = 1 }"),
               "must have as many entries as z_points");
  EXPECT_ERROR(cfgOf(std::string(BASE) + "bunch { charge = 1 pC; gamma = 100; sigma_x = 1 um; zeta_center = 5 um }"),
               "needs key 'sigma_z'");
  EXPECT_ERROR(cfgOf(std::string(BASE) + "output { every = 2.5 }"), "'2.5' is not an integer");
  {
    // m_max = 0 is not allowed with a laser
    std::string deck = BASE;
    deck.insert(deck.find("  dr ="), "  m_max = 0\n");
    EXPECT_ERROR(cfgOf(deck + "laser { a0 = 1; waist = 8 um; tau_fwhm = 30 fs; zeta_center = 30 um }"),
                 "set grid.m_max >= 1");
  }
  // line numbers: BASE has 11 lines, the misspelled key is in line 13
  EXPECT_ERROR(cfgOf(std::string(BASE) + "output {\n  evry = 3\n}"), "test.deck:13:");

  // ------------------------------------------------------------ overrides
  {
    Deck d = Deck::parseString(std::string(BASE) +
                                   "laser { a0 = 1; waist = 8 um; tau_fwhm = 30 fs; zeta_center = 30 um }\n"
                                   "laser { a0 = 1; waist = 8 um; tau_fwhm = 30 fs; zeta_center = 20 um }",
                               "t");
    d.applyOverride("grid.z_end=2mm");
    d.applyOverride("laser[1].a0=3.5");
    d.applyOverride("laser[0].z_focus=0.5 mm");  // adds a key
    Config c = Config::fromDeck(d);
    CHECK_NEAR(c.grid.z_end, 2000.0 / k0, 1e-12);
    CHECK_NEAR(c.lasers[1].a0, 3.5, 1e-15);
    CHECK_NEAR(c.lasers[0].a0, 1.0, 1e-15);
    CHECK_NEAR(c.lasers[0].z_focus, 500.0 / k0, 1e-12);
    EXPECT_ERROR(d.applyOverride("laser[2].a0=1"), "has no block laser[2]");
    EXPECT_ERROR(d.applyOverride("a0=1"), "must have the form block.key=value");
    // dump -> parse gives the same configuration
    Config c2 = Config::fromDeck(Deck::parseString(d.dump(), "dump"));
    CHECK(c2.grid.Nz == c.grid.Nz && c2.lasers.size() == 2);
    CHECK(c2.lasers[1].a0 == c.lasers[1].a0 && c2.grid.z_end == c.grid.z_end);
  }

  // ------------------------------------------------------------ lambda0 from the first laser
  {
    Config c = cfgOf(R"(grid { radial_grid = uniform; dr = 1 um; nr = 10; dzeta = 50 nm; nzeta = 400
                              dz = 5 um; z_end = 1 mm }
                       laser { a0 = 1; wavelength = 1.0 um; waist = 5 um; tau_fwhm = 20 fs; zeta_center = 10 um })");
    CHECK_NEAR(c.units.lambda0_um, 1.0, 1e-15);
    CHECK_NEAR(c.lasers[0].omega, 1.0, 1e-14);
    EXPECT_ERROR(cfgOf("grid { lambda0 = 6.28 k0^-1; radial_grid = uniform; dr = 1 um; nr = 10; dzeta = 50 nm\n"
                       "nzeta = 400; dz = 5 um; z_end = 1 mm }"),
                 "must be given in physical length units");
  }

  // ------------------------------------------------------------ ion species
  {
    const std::string P = std::string(BASE) + "plasma { density = 1e18 cm^-3 }\n";
    Config c = cfgOf(P + "ion { element = N; charge_state = 5; fraction = 0.3 }\n"
                         "ion { element = H; fraction = 0.5 }\n"
                         "ion { mass = 4 u; charge_state = 2; name = alpha; fraction = 0.2 }\n");
    CHECK(c.ions.size() == 3);
    CHECK(c.ions[0].name == "N5+" && c.ions[0].Z == 5);
    CHECK_NEAR(c.ions[0].mass, 14.007 * U_OVER_ME, 1e-12);
    CHECK_NEAR(c.ions[0].qm(), 5.0 / (14.007 * U_OVER_ME), 1e-12);
    CHECK(c.ions[1].name == "H1+" && c.ions[1].Z == 1);
    CHECK(c.ions[2].name == "alpha");
    CHECK_NEAR(c.ions[2].mass, 4.0 * U_OVER_ME, 1e-12);
    CHECK_NEAR(toCode(1.0, "m_p", Dim::Mass, u), MP_OVER_ME, 1e-14);
    EXPECT_ERROR(cfgOf(P + "ion { element = N }"), "ion.charge_state: is required for N (1 ... 7)");
    EXPECT_ERROR(cfgOf(P + "ion { element = N; charge_state = 8 }"), "exceeds the atomic number 7 of N");
    EXPECT_ERROR(cfgOf(P + "ion { element = H; mass = 1 u }"), "give either 'element' or 'mass'");
    EXPECT_ERROR(cfgOf(P + "ion { element = Xy }"), "'Xy' is not allowed");
    EXPECT_ERROR(cfgOf(P + "ion { mass = 1 um }"), "unit 'um' is a length, but a mass is expected");
    EXPECT_ERROR(cfgOf(P + "ion { element = H; fraction = 0.7 }\nion { element = D; fraction = 0.6 }"),
                 "ion fractions add up to 1.3");
    EXPECT_ERROR(cfgOf(std::string(BASE) + "ion { element = H }"), "needs a 'plasma' block");
  }

  // ------------------------------------------------------------ deposition shape on a stretched grid
  {
    Config c = cfgOf(R"(grid {
      lambda0 = 0.8 um;  dr_min = 0.5 k0^-1;  dr_max = 10 k0^-1
      r_transition = 30 k0^-1;  r_max = 600 k0^-1
      dzeta = 5 k0^-1;  nzeta = 10;  dz = 10 k0^-1;  z_end = 10 k0^-1 })");
    RadialGrid g;
    g.build(c.grid);
    // volumes equal the solver measures (discrete Gauss law)
    CHECK(g.maxVolumeMismatch < 0.02);  // only a clamped axis interval may deviate
    int exact = 0;
    for (int i = 1; i < g.N - 1; ++i) exact += std::fabs(g.Vc[i] / g.Ac[i] - 1.0) < 1e-12;
    CHECK(exact == g.N - 2);
    // a uniform density (sampled uniformly in r^2) is deposited uniformly
    std::vector<double> Qc(g.N, 0.0);
    const double R = g.rho[g.N] * 0.999;
    const int S = 2000000;
    for (int k = 0; k < S; ++k) {
      const double r = R * std::sqrt((k + 0.5) / S);
      auto w = g.locC(r);
      if (w.mode == 0) { Qc[w.a] += w.wa; Qc[w.a + 1] += w.wb; }
      else Qc[w.a] += 1.0;
    }
    const double dens = S / (0.5 * R * R);  // particles per unit "area" (int r dr)
    double worst = 0.0;
    for (int i = 0; i < g.N - 2; ++i) worst = std::max(worst, std::fabs(Qc[i] / g.Vc[i] / dens - 1.0));
    CHECK(worst < 2e-3);
    if (worst >= 2e-3) std::printf("   uniform density deviation %.2e\n", worst);
    // source smoothing (plasma.smooth_length): the vector filter commutes with the discrete
    // divergence, div F_vec j = F_scalar div j, for every mode; constants are kept away from the wall
    {
      const int M = 2;
      FieldSolver fs;
      fs.setup(g, M);
      SourceFilter F;
      F.setup(g, fs, M, 2.0);
      CHECK(F.on());
      for (int m = 0; m <= M; ++m) {
        const cd im(0.0, double(m));
        std::vector<cd> jr(size_t(g.N) + 1, 0.0), jp(size_t(g.N), 0.0);
        for (int i = 1; i <= g.N; ++i) jr[i] = cd(std::sin(0.37 * i) + 0.2 * std::cos(1.7 * i), 0.3 * std::sin(2.1 * i));
        for (int i = 0; i < g.N; ++i) jp[i] = cd(std::cos(0.53 * i), std::sin(0.11 * i * i));
        auto div = [&](const std::vector<cd>& r, const std::vector<cd>& p) {
          std::vector<cd> d(size_t(g.N));
          for (int i = 0; i < g.N; ++i)
            d[i] = (g.rho[i + 1] * r[i + 1] - g.rho[i] * r[i]) / g.Ac[i] + (im / g.c[i]) * p[i];
          return d;
        };
        std::vector<cd> d0 = div(jr, jp);
        F.scalar(m, d0.data());
        F.vector(m, jr.data(), jp.data());
        const std::vector<cd> d1 = div(jr, jp);
        double err = 0.0, ref = 0.0;
        for (int i = 0; i < g.N; ++i) { err = std::max(err, std::abs(d1[i] - d0[i])); ref = std::max(ref, std::abs(d0[i])); }
        CHECK(err < 1e-11 * ref);
        if (!(err < 1e-11 * ref)) std::printf("   smoothing m = %d: div mismatch %.2e of %.2e\n", m, err, ref);
      }
      std::vector<double> one(size_t(g.N), 1.0);
      F.scalar0(one.data());
      double dev = 0.0;
      for (int i = 0; g.c[i] < 0.5 * g.rho[g.N]; ++i) dev = std::max(dev, std::fabs(one[i] - 1.0));
      CHECK(dev < 1e-10);
      // a ring is spread over ~a: the filtered delta has unit integral (inside the wall)
      std::vector<double> ring(size_t(g.N), 0.0);
      const int ir = 40;
      ring[ir] = 1.0 / g.Ac[ir];
      F.scalar0(ring.data());
      double q = 0.0;
      for (int i = 0; i < g.N; ++i) q += ring[i] * g.Ac[i];
      CHECK_NEAR(q, 1.0, 1e-10);
      SourceFilter off;
      off.setup(g, fs, M, 0.0);
      CHECK(!off.on());
    }
  }

  // ------------------------------------------------------------ backends
  {
    const std::string G = "grid { lambda0 = 0.8 um; radial_grid = uniform; dr = 0.1 um; nr = 10; dzeta = 30 nm; nzeta = 10;"
                          " dz = 5 um; z_end = 5 um }\n";
    using B = NumericsConfig::Backend;
#if defined(GEMPIC_HAVE_CUDA)
    CHECK(cfgOf(G).numerics.backend == B::Cuda);  // default of a CUDA build
    CHECK(cfgOf(G + "numerics { backend = serial }").numerics.backend == B::Serial);
    EXPECT_ERROR(cfgOf(G + "numerics { backend = openmp }"), "this is a CUDA build");
#else
    CHECK(cfgOf(G).numerics.backend == B::Serial);  // default of a CPU build
    CHECK(cfgOf(G + "numerics { backend = openmp }").numerics.backend == B::OpenMP);
    EXPECT_ERROR(cfgOf(G + "numerics { backend = cuda }"), "no GPU support");
#endif
    EXPECT_ERROR(cfgOf(G + "numerics { backend = opencl }"), "backend");
    // beam current filter along zeta
    CHECK(cfgOf(G).numerics.beam_filter == 1);
    CHECK(cfgOf(G + "numerics { beam_filter = 0 }").numerics.beam_filter == 0);
    CHECK(cfgOf(G + "numerics { beam_filter = 3 }").numerics.beam_filter == 3);
    EXPECT_ERROR(cfgOf(G + "numerics { beam_filter = 9 }"), "beam_filter");
    EXPECT_ERROR(cfgOf(G + "numerics { beam_filter = -1 }"), "beam_filter");
  }

  // ------------------------------------------------------------ output format
  {
    const std::string G = "grid { lambda0 = 0.8 um; radial_grid = uniform; dr = 0.1 um; nr = 10; dzeta = 30 nm;"
                          " nzeta = 10; dz = 5 um; z_end = 5 um }\n";
    const OutputConfig o0 = cfgOf(G).output;
    CHECK(o0.text && !o0.hdf5);  // default: text only
    EXPECT_ERROR(cfgOf(G + "output { format = netcdf }"), "format");
#if defined(GEMPIC_HAVE_HDF5)
    const OutputConfig o1 = cfgOf(G + "output { format = hdf5; hdf5_precision = single; hdf5_compression = 4 }").output;
    CHECK(!o1.text && o1.hdf5 && o1.hdf5_single && o1.hdf5_compression == 4);
    const OutputConfig o2 = cfgOf(G + "output { format = both }").output;
    CHECK(o2.text && o2.hdf5 && !o2.hdf5_single && o2.hdf5_compression == 0);
    EXPECT_ERROR(cfgOf(G + "output { format = hdf5; hdf5_compression = 12 }"), "0 ... 9");
    const OutputConfig o3 = cfgOf(G + "output { checkpoint_every = 50; checkpoint_keep = 3; walltime_hours = 23.5 }").output;
    CHECK(o3.checkpoint_every == 50 && o3.checkpoint_keep == 3 && o3.walltime_hours == 23.5);
    CHECK(o0.checkpoint_every == 0 && o0.checkpoint_keep == 2 && o0.walltime_hours == 0.0);  // defaults
    EXPECT_ERROR(cfgOf(G + "output { checkpoint_every = -1 }"), "must be >= 0");
#else
    EXPECT_ERROR(cfgOf(G + "output { format = hdf5 }"), "no HDF5 support");
    EXPECT_ERROR(cfgOf(G + "output { checkpoint_every = 10 }"), "need HDF5 support");
#endif
  }

  // ------------------------------------------------------------ free-hand mapping
  {
    StepMapping M1;
    CHECK(M1.set({0.0}, {0.5}).empty());  // one point: uniform grid
    CHECK_NEAR(M1.coord(10), 5.0, 1e-14);
    CHECK_NEAR(M1.index(7.5), 15.0, 1e-14);
    StepMapping M;
    CHECK(M.set({1.0, 3.0, 10.0}, {0.1, 0.1, 1.0}).empty());
    CHECK_NEAR(M.step(0.5), 0.1, 1e-14);     // constant before the first point
    CHECK_NEAR(M.step(2.0), 0.1, 1e-14);
    CHECK_NEAR(M.step(6.5), 0.55, 1e-14);    // linear in between
    CHECK_NEAR(M.step(50), 1.0, 1e-14);      // constant after the last point
    CHECK_NEAR(M.coord(30), 3.0, 1e-12);     // 30 cells of 0.1 up to r = 3
    // cell size around x equals h(x): inside a segment up to b^2/24 (curvature
    // of the exponential), in a cell containing a kink up to ~b/8
    double worstInv = 0, worstInside = 0, worstKink = 0;
    for (double xi = 0.0; xi < 60.0; xi += 0.37) {
      const double x = M.coord(xi);
      worstInv = std::max(worstInv, std::fabs(M.index(x) - xi));
      const double xl = M.coord(xi - 0.5), xr = M.coord(xi + 0.5);
      const double err = std::fabs((xr - xl) / M.step(x) - 1.0);
      bool kink = false;
      for (double p : {1.0, 3.0, 10.0}) kink = kink || (xl < p && p < xr);
      (kink ? worstKink : worstInside) = std::max(kink ? worstKink : worstInside, err);
    }
    CHECK(worstInv < 1e-10);
    CHECK(worstInside < 1e-3);
    CHECK(worstKink < 0.03);
    CHECK_NEAR(M.maxRelativeChange(), 0.9 / 7.0, 1e-14);
    StepMapping bad;
    CHECK(!bad.set({0.0, 1.0}, {0.1}).empty());
    CHECK(!bad.set({2.0, 1.0}, {0.1, 0.1}).empty());
    CHECK(!bad.set({0.0, 1.0}, {0.1, 0.0}).empty());

    Config c = cfgOf(R"(grid {
      lambda0 = 0.8 um;  radial_grid = table
      r_points  = 0, 1, 3, 20 um
      dr_points = 40, 40, 400, 400 nm
      r_max = 24 um;  dzeta = 30 nm;  nzeta = 100;  dz = 5 um;  z_end = 50 um })");
    const auto& g = c.grid;
    CHECK(g.radial == GridConfig::Radial::Table && !g.stretched);
    CHECK_NEAR(g.dr_min, 0.04 / k0, 1e-12);
    CHECK_NEAR(g.dr_max, 0.4 / k0, 1e-12);
    CHECK_NEAR(g.mapping(25), 1.0 / k0, 1e-9);  // 25 cells of 40 nm
    CHECK(g.mapping(g.Nr) >= 24.0 / k0 * (1 - 1e-12) && g.mapping(g.Nr - 1) < 24.0 / k0);
    RadialGrid rg;
    rg.build(c.grid);
    CHECK_NEAR(rg.rho[1] - rg.rho[0], 0.04 / k0, 1e-9);
    CHECK_NEAR(rg.rho[rg.N] - rg.rho[rg.N - 1], 0.4 / k0, 1e-9);
    CHECK(rg.maxVolumeMismatch < 0.02);
    std::vector<std::string> warn;
    Config::fromDeck(Deck::parseString(R"(grid {
      lambda0 = 0.8 um;  radial_grid = table;  r_points = 0, 0.2 um;  dr_points = 10, 100 nm
      nr = 50;  dzeta = 30 nm;  nzeta = 100;  dz = 5 um;  z_end = 50 um })", "t"), &warn);
    bool warned = false;
    for (const auto& w : warn) warned = warned || w.find("from one cell to the next") != std::string::npos;
    CHECK(warned);  // 45 % per cell
    EXPECT_ERROR(cfgOf("grid { lambda0 = 0.8 um; radial_grid = table; r_points = 0, 1 um; dr_points = 40 nm\n"
                       "nr = 10; dzeta = 30 nm; nzeta = 100; dz = 5 um; z_end = 50 um }"),
                 "must have as many entries as r_points");
    EXPECT_ERROR(cfgOf("grid { lambda0 = 0.8 um; radial_grid = table; r_points = 0, 1 um; dr_points = 40, 0 nm\n"
                       "nr = 10; dzeta = 30 nm; nzeta = 100; dz = 5 um; z_end = 50 um }"),
                 "grid steps must be positive");
    EXPECT_ERROR(cfgOf("grid { lambda0 = 0.8 um; radial_grid = table; r_points = 0, 1 um; dr_points = 40, 40 nm\n"
                       "dr_min = 40 nm; nr = 10; dzeta = 30 nm; nzeta = 100; dz = 5 um; z_end = 50 um }"),
                 "unknown key 'dr_min'");
  }

  // ------------------------------------------------------------ free-hand zeta and z grids
  {
    const double k0 = 2.0 * PI / 0.8;
    const std::string R = "grid { lambda0 = 0.8 um; radial_grid = uniform; dr = 0.1 um; nr = 10\n";
    Config c = cfgOf(R + "zeta_grid = table; zeta_points = 0, 10, 20 um; dzeta_points = 100, 100, 25 nm\n"
                         "window_length = 30 um\n"
                         "z_grid = table; z_points = 0, 1 mm; dz_points = 10, 50 um; z_end = 2 mm }");
    const auto& g = c.grid;
    CHECK(g.zetaTable && g.zTable);
    CHECK_NEAR(g.dzeta, 0.025 * k0, 1e-12);
    CHECK_NEAR(g.dzeta_max, 0.1 * k0, 1e-12);
    CHECK_NEAR(g.zetaNode(100), 10.0 * k0, 1e-9);           // 100 steps of 100 nm
    CHECK(g.window() >= 30.0 * k0 * (1 - 1e-12) && g.zetaNode(g.Nz - 1) < 30.0 * k0);  // rounded up
    CHECK_NEAR(g.zetaNode(g.Nz) - g.zetaNode(g.Nz - 1), 0.025 * k0, 1e-9);
    CHECK_NEAR(g.zetaStep(15.0 * k0), 0.0625 * k0, 1e-12);  // linear between the points
    // the step from s is the cell X(xi(s) + 1) - s: about dz at its midpoint
    CHECK(std::fabs(g.stepAt(0.0) / (10.0 * k0) - 1.0) < 0.03);
    CHECK_NEAR(g.stepAt(1500.0 * k0), 50.0 * k0, 1e-6);       // constant beyond the last point
    {
      // the levels follow the table: s_k = X(k)
      double s = 0.0, worst = 0.0;
      for (int k = 0; k < 40; ++k) {
        const double h = g.stepAt(s);
        worst = std::max(worst, std::fabs(h - g.stable.step(s + 0.5 * h)) / h);  // midpoint step
        s += h;
      }
      CHECK(worst < 0.01);
    }
    CHECK(g.stepsFrom(0.0) > 40 && g.stepsFrom(0.0) < 200);
    // uniform: unchanged arithmetic
    Config u = cfgOf(R + "dzeta = 30 nm; nzeta = 100; dz = 5 um; z_end = 50 um }");
    CHECK(!u.grid.zetaTable && !u.grid.zTable);
    CHECK(u.grid.zetaNode(7) == 7 * u.grid.dzeta && u.grid.stepAt(123.0) == u.grid.dz);
    CHECK(u.grid.stepsFrom(0.0) == 10);
    // the longitudinal grid of the solver
    ZetaGrid zg;
    std::vector<double> nodes(g.Nz + 1);
    for (int n = 0; n <= g.Nz; ++n) nodes[n] = g.zetaNode(n);
    zg.build(nodes, false, g.dzeta);
    const ZetaView zv = zg.view();
    double w = 0;
    const int na = zv.cell(nodes[150] + 0.25 * zg.hs[150], w);
    CHECK(na == 150 && std::fabs(w - 0.25) < 1e-12);
    double f = 0;
    const int j = zv.beamCell(zg.zc[151], f);  // centre of half-cell 150
    CHECK(j == 151 && std::fabs(f) < 1e-12);
    CHECK(zv.nearest(nodes[150] + 0.4 * zg.hs[150]) == 150 && zv.nearest(nodes[150] + 0.6 * zg.hs[150]) == 151);
    CHECK(zv.joinCell(nodes[150]) == 149 && zv.joinCell(nodes[150] + 1e-9) == 150);
    CHECK(std::fabs(zv.wLo(150) + zv.wHi(150) - 2.0) < 1e-14);
    // linear interpolation of half-cell values to slice n is exact for a linear function
    {
      const int n = 250;
      const double lin = 0.5 * (zv.wLo(n) * zv.centre(n - 1) + zv.wHi(n) * zv.centre(n));
      CHECK(std::fabs(lin - zv.node(n)) < 1e-9);
    }
    std::vector<std::string> warn;
    Config::fromDeck(Deck::parseString(R + "zeta_grid = table; zeta_points = 0, 1 um; dzeta_points = 10, 400 nm\n"
                                           "nzeta = 100; dz = 5 um; z_end = 50 um }", "t"), &warn);
    bool warned = false;
    for (const auto& s : warn) warned = warned || s.find("zeta step changes") != std::string::npos;
    CHECK(warned);
    EXPECT_ERROR(cfgOf(R + "zeta_grid = table; zeta_points = 0, 1 um; dzeta_points = 10 nm\n"
                           "nzeta = 100; dz = 5 um; z_end = 50 um }"),
                 "must have as many entries as zeta_points");
    EXPECT_ERROR(cfgOf(R + "zeta_grid = table; zeta_points = 0 um; dzeta_points = 10 nm; dzeta = 10 nm\n"
                           "nzeta = 100; dz = 5 um; z_end = 50 um }"),
                 "unknown key 'dzeta'");
    EXPECT_ERROR(cfgOf(R + "dzeta = 30 nm; nzeta = 100; z_grid = table; z_points = 0, 1 mm; dz_points = 5 um\n"
                           "z_end = 50 um }"),
                 "must have as many entries as z_points");
    EXPECT_ERROR(cfgOf(R + "dzeta = 30 nm; nzeta = 100; z_grid = table; z_points = 0 mm; dz_points = 0 um\n"
                           "z_end = 50 um }"),
                 "grid steps must be positive");
    EXPECT_ERROR(cfgOf(R + "dzeta = 30 nm; nzeta = 100; z_grid = adaptive; dz = 5 um; z_end = 50 um }"), "z_grid");
  }

  // ------------------------------------------------------------ particle pushers
  {
    // gamma = 1000 along z in E_x = beta B_y: E + v x B = 0 (Vay exact, Boris not)
    const double uz0 = std::sqrt(1000.0 * 1000.0 - 1.0), beta = uz0 / 1000.0;
    EB f{beta * 1e-3, 0, 0, 0, 1e-3, 0};
    double ax = 0, ay = 0, az = uz0, bx = 0, by = 0, bz = uz0;
    for (int k = 0; k < 100; ++k) {
      vay(ax, ay, az, f, -1.0, 5e5);
      boris(bx, by, bz, f, -1.0, 5e5);
    }
    CHECK(std::fabs(ax) < 1e-9 && std::fabs(az - uz0) < 1e-9 * uz0);
    CHECK(std::fabs(bx) > 1e-3);  // the Boris push drifts (the reason for the Vay push)
    // pure magnetic field: |u| conserved, gyration as Boris for small steps
    double cx = 1.0, cy = 0.0, cz = 0.0, dx = 1.0, dy = 0.0, dz = 0.0;
    EB m{0, 0, 0, 0, 0, 0.01};
    for (int k = 0; k < 1000; ++k) { vay(cx, cy, cz, m, -1.0, 0.1); boris(dx, dy, dz, m, -1.0, 0.1); }
    CHECK(std::fabs(std::hypot(cx, cy) - 1.0) < 1e-12);
    CHECK(std::fabs(cx - dx) < 1e-4 && std::fabs(cy - dy) < 1e-4);
  }

  // ------------------------------------------------------------ symmetric tridiagonal eigenproblem
  {
    const int n = 40;
    std::vector<double> d(n, 2.0), e(n, -1.0), Z;
    CHECK(symTridiagEigen(d, e, Z));
    std::vector<double> ev = d, ex(n);
    std::sort(ev.begin(), ev.end());
    for (int k = 0; k < n; ++k) ex[k] = 2.0 - 2.0 * std::cos((k + 1) * PI / (n + 1));
    double worst = 0, ortho = 0, resid = 0;
    for (int k = 0; k < n; ++k) worst = std::max(worst, std::fabs(ev[k] - ex[k]));
    for (int a = 0; a < n; ++a)
      for (int b = 0; b < n; ++b) {
        double sdot = 0;
        for (int i = 0; i < n; ++i) sdot += Z[size_t(i) * n + a] * Z[size_t(i) * n + b];
        ortho = std::max(ortho, std::fabs(sdot - (a == b ? 1.0 : 0.0)));
      }
    for (int k = 0; k < n; ++k)  // T z = lambda z
      for (int i = 0; i < n; ++i) {
        double t = 2.0 * Z[size_t(i) * n + k];
        if (i > 0) t -= Z[size_t(i - 1) * n + k];
        if (i + 1 < n) t -= Z[size_t(i + 1) * n + k];
        resid = std::max(resid, std::fabs(t - d[k] * Z[size_t(i) * n + k]));
      }
    CHECK(worst < 1e-12 && ortho < 1e-12 && resid < 1e-12);
  }

  // ------------------------------------------------------------ undulator
  {
    const std::string G = "grid { lambda0 = 0.8 um; radial_grid = uniform; dr = 10 um; nr = 50; dzeta = 25 nm; nzeta = 100;"
                          " dz = 0.5 mm; z_end = 0.1 m }\n"
                          "bunch { name = b; charge = 1 fC; gamma = 100; sigma_x = 50 um; sigma_z = 0.5 um;"
                          " zeta_center = 1.5 um; macroparticles = 100 }\n";
    Config c = cfgOf(G + "undulator { period = 14.2222 mm; K = 0.5; periods = 50 }");
    CHECK(c.undulators.size() == 1);
    const auto& U = c.undulators[0];
    const double k0 = 2.0 * PI / 0.8;
    CHECK(U.type == UndulatorConfig::Type::Planar && U.focusing);
    CHECK_NEAR(U.length, 50 * 14222.2 * k0, 1e-6 * U.length);
    CHECK_NEAR(U.ramp, U.period, 1e-12);
    CHECK_NEAR(U.resonance(100.0) / k0, 0.8, 1e-5);           // lambda_r = 0.8 um
    CHECK_NEAR(U.B0() * c.units.B0_T(), 0.3765, 2e-4);         // K = 0.934 B[T] lambda_u[cm]
    // the field: B_y = B0 cos(k_u z) on axis, zero outside, ramp
    UndulatorV v{0, 1, U.ku(), U.B0(), U.z_start, U.z_start + U.length, U.ramp, 0, 0, 0.0};
    ExternalV X{1, &v};
    EB f{};
    addExternal(X, 0.0, 0.0, 10.0 * U.period, f);
    CHECK_NEAR(f.By, U.B0(), 1e-9 * U.B0());
    EB f2{};
    addExternal(X, 0.0, 0.0, 0.5 * U.period, f2);             // half-way up the ramp, cos = -1
    CHECK_NEAR(f2.By, -0.5 * U.B0(), 1e-9 * U.B0());
    EB f3{};
    addExternal(X, 0.0, 0.0, U.length + 1.0, f3);
    CHECK(f3.By == 0.0 && f3.Bx == 0.0);
    EB f4{};
    addExternal(X, 0.0, 20.0, 10.25 * U.period, f4);         // off axis: B_z = -B0 sinh(k_u y) sin(psi)
    CHECK_NEAR(f4.Bz, -U.B0() * std::sinh(U.ku() * 20.0), 1e-9 * U.B0());
    Config h = cfgOf(G + "undulator { type = helical; period = 1 cm; K = 1; length = 0.2 m; ramp_periods = 2 }");
    CHECK(h.undulators[0].type == UndulatorConfig::Type::Helical);
    CHECK_NEAR(h.undulators[0].resonance(100.0) / k0, 1e4 * 2.0 / 2e4, 1e-9);  // (1 + K^2) lambda_u / (2 gamma^2)
    std::vector<std::string> warn;
    Config::fromDeck(Deck::parseString(G + "undulator { period = 5 mm; K = 0.2; periods = 20 }", "t"), &warn);
    bool wdz = false, wzeta = false;
    for (const auto& w : warn) {
      wdz = wdz || w.find("grid.dz") != std::string::npos;
      wzeta = wzeta || w.find("resonant wavelength") != std::string::npos;
    }
    CHECK(wdz && wzeta);  // dz = 0.5 mm > period / 16; lambda_r = 0.255 um = 10 steps of 25 nm (< 12)
    EXPECT_ERROR(cfgOf(G + "undulator { period = 1 cm; K = 1; periods = 10; ramp_periods = 6 }"), "ramp_periods");
    EXPECT_ERROR(cfgOf(G + "undulator { period = 1 cm; K = 1; periods = 10; length = 10 cm }"), "not both");
    EXPECT_ERROR(cfgOf(G + "undulator { type = helical; period = 1 cm; K = 1; periods = 10; ramp_periods = 0 }"),
                 "helical undulator needs ramp_periods");
    EXPECT_ERROR(cfgOf(G + "undulator { type = helical; period = 1 cm; K = 1; periods = 10; natural_focusing = no }"),
                 "planar undulators only");
    EXPECT_ERROR(cfgOf(G + "undulator { period = 1 cm; K = -1; periods = 10 }"), "must not be negative");
  }

  // ------------------------------------------------------------ NIST data and ADK rates
  {
    // table integrity: Znuc levels, increasing energies
    int bad = 0;
    for (const auto& e : elements()) {
      if (int(e.ip_eV.size()) != e.Znuc) ++bad;
      for (size_t k = 1; k < e.ip_eV.size(); ++k)
        if (!(e.ip_eV[k] > e.ip_eV[k - 1])) ++bad;
    }
    CHECK(bad == 0);
    CHECK(elements().size() == 21);
    CHECK(findElement("n") && findElement("N")->ip_eV[5] == 552.06741);
    CHECK(findElement("H")->ip_eV[0] == 13.598434599702);
    CHECK(findElement("Ar")->ip_eV[0] == 15.7596119 && findElement("Ar")->ip_eV[17] == 4426.22407);
    CHECK(findElement("Au")->ip_eV.size() == 79 && findElement("Au")->ip_eV[78] == 93254.62);
    CHECK(findElement("Xx") == nullptr);
    // ADK = Landau rate for hydrogen: w = 4 omega_a / F exp(-2 / 3F), F = E / E_a
    for (double F : {0.02, 0.05, 0.1}) {
      const double landau = 4.0 * OMEGA_ATOMIC / F * std::exp(-2.0 / (3.0 * F));
      CHECK_NEAR(adkRateSI(U_H_EV, 1, F * E_ATOMIC_VM), landau, 1e-12);
    }
    // hydrogen-like ion Z = 2 (U = 4 U_H): w = 4 omega_a Z^5 / F exp(-2 Z^3 / 3F)
    for (double F : {0.2, 0.5}) {
      const double landau = 4.0 * OMEGA_ATOMIC * 32.0 / F * std::exp(-16.0 / (3.0 * F));
      CHECK_NEAR(adkRateSI(4.0 * U_H_EV, 2, F * E_ATOMIC_VM), landau, 1e-12);
    }
    // code units: rate(E) = w_SI(E E0) / omega0
    AdkRates R;
    R.setup(findElement("N")->ip_eV, 5, 7, u, 1000.0);
    const double E0 = 2 * PI * MC2_MEV * 1e6 / (0.8e-6), om0 = 2 * PI * C_LIGHT / 0.8e-6;
    for (double E : {1.5, 2.0, 3.0}) {
      CHECK_NEAR(R.rate(5, E), adkRateSI(552.06741, 6, E * E0) / om0, 1e-12);
      CHECK_NEAR(R.rate(6, E), adkRateSI(667.0461377, 7, E * E0) / om0, 1e-12);
    }
    CHECK(R.rate(4, 2.0) == 0.0 && R.rate(7, 2.0) == 0.0);  // levels outside Z0 ... Zmax-1
    CHECK(R.rate(5, 0.1) == 0.0);                          // far below the threshold
    CHECK(R.minField() > 0.1 && R.minField() < 1.5);
    // beyond the peak of the ADK expression the rate stays at its maximum
    const double big = R.rate(5, 1e4), bigger = R.rate(5, 1e5);
    CHECK(big > 0 && bigger == big);
  }

  // ------------------------------------------------------------ ionizable ion blocks
  {
    const std::string P = std::string(BASE) + "plasma { density = 1e18 cm^-3 }\n";
    Config c = cfgOf(P + "ion { element = N; charge_state = 5; ionization = adk; density = 1e16 cm^-3; mobile = false }\n"
                         "ion { element = He; ionization = adk; density = 1e17 cm^-3 }\n");
    CHECK(c.ions[0].ionization && c.ions[0].Z == 5 && c.ions[0].Zmax == 7 && !c.ions[0].mobile);
    CHECK_NEAR(c.ions[0].atoms, 0.01, 1e-12);
    CHECK_NEAR(c.ions[0].fraction(), 0.05, 1e-12);
    CHECK(c.ions[1].Z == 0 && c.ions[1].Zmax == 2 && c.ions[1].name == "He");  // neutral by default
    CHECK(c.ions[1].fraction() == 0.0);
    Config c2 = cfgOf(P + "ion { mass = 2 u; ionization = adk; ionization_potentials = 10, 30 eV; density = 1e17 cm^-3 }");
    CHECK(c2.ions[0].Znuc == 2 && c2.ions[0].Zmax == 2);
    CHECK_NEAR(c2.ions[0].ip_eV[1], 30.0, 1e-12);
    Config c3 = cfgOf(std::string(BASE) + "plasma { density = 1e18 cm^-3; free_electrons = false }\n"
                                          "ion { element = Ar; ionization = adk; max_charge_state = 8; density = 1e18 cm^-3 }");
    CHECK(!c3.plasma.free_electrons && c3.ions[0].Zmax == 8);
    EXPECT_ERROR(cfgOf(P + "ion { element = He; ionization = adk }"), "neutral atoms: give their 'density'");
    EXPECT_ERROR(cfgOf(P + "ion { element = He; charge_state = 0 }"), "0 (neutral atoms) needs ionization = adk");
    EXPECT_ERROR(cfgOf(P + "ion { element = N; charge_state = 5; ionization = adk; max_charge_state = 8 }"),
                 "must lie between charge_state and 7");
    EXPECT_ERROR(cfgOf(P + "ion { element = N; charge_state = 5; max_charge_state = 7 }"),
                 "is only used with ionization = adk");
    EXPECT_ERROR(cfgOf(P + "ion { mass = 2 u; ionization = adk; density = 1e17 cm^-3 }"),
                 "needs ionization energies");
    EXPECT_ERROR(cfgOf(P + "ion { element = He; ionization = adk; ionization_potentials = 20 eV; density = 1e17 cm^-3 }"),
                 "needs 2 values for He");
    EXPECT_ERROR(cfgOf(P + "ion { mass = 2 u; ionization = adk; ionization_potentials = 30, 10 eV; density = 1e17 cm^-3 }"),
                 "must be positive and increasing");
    EXPECT_ERROR(cfgOf(std::string(BASE) + "plasma { density = 1e18 cm^-3; free_electrons = false }\n"
                                           "ion { element = N; charge_state = 5; ionization = adk; density = 1e16 cm^-3 }"),
                 "plasma.free_electrons = false, but ion blocks with charge_state > 0");
    EXPECT_ERROR(cfgOf(P + "ion { element = N; charge_state = 5; density = 1e18 cm^-3 }"), "ion fractions add up to 5");
  }

  // ------------------------------------------------------------ overrides may create a missing block
  {
    Deck d = Deck::parseString(BASE, "t");
    d.applyOverride("numerics.seed=7");
    Config c = Config::fromDeck(d);
    CHECK(c.numerics.seed == 7);
  }

  // ------------------------------------------------------------ function parser
  {
    Units u;
    u.lambda0_um = 0.8;
    const double k0 = 2.0 * PI / 0.8;  // code length per um
    const std::map<std::string, double> none;
    auto ev = [&](const std::string& t) { return Expr(t, {}, none, u)(nullptr); };
    CHECK_NEAR(ev("1 + 2*3"), 7.0, 1e-15);
    CHECK_NEAR(ev("2^3^2"), 512.0, 1e-15);         // right-associative
    CHECK_NEAR(ev("-2^2"), -4.0, 1e-15);           // -(2^2)
    CHECK_NEAR(ev("2**-1"), 0.5, 1e-15);
    CHECK_NEAR(ev("(1 < 2) && (3 >= 3)"), 1.0, 1e-15);
    CHECK_NEAR(ev("1 > 2 ? 5 : 6"), 6.0, 1e-15);
    CHECK_NEAR(ev("if(0, 1, 2) + !0"), 3.0, 1e-15);
    CHECK_NEAR(ev("gauss(1, 1)"), std::exp(-0.5), 1e-15);
    CHECK_NEAR(ev("rect(0.4, 1) + rect(0.6, 1)"), 1.0, 1e-15);
    CHECK_NEAR(ev("atan2(1, 1)"), PI / 4, 1e-15);
    CHECK_NEAR(ev("min(3, max(1, 2)) + step(0) + sign(-3)"), 2.0, 1e-15);
    CHECK_NEAR(ev("sqrt(16) + ln(exp(2)) + log10(100) + abs(-1)"), 9.0, 1e-15);
    CHECK_NEAR(ev("10 um"), 10.0 * k0, 1e-14);     // numbers with units -> code units
    CHECK_NEAR(ev("10um / (2 um)"), 5.0, 1e-14);
    CHECK_NEAR(ev("1e-3 mm"), 1.0 * k0, 1e-14);
    CHECK_NEAR(ev("50%"), 0.5, 1e-15);
    CHECK_NEAR(ev("2 * pi"), 2.0 * PI, 1e-15);
    Expr f("exp(-(r/w)^2) * (z > 0)", {"z", "r"}, {{"w", 2.0}}, u);
    CHECK_NEAR(f({1.0, 2.0}), std::exp(-1.0), 1e-15);
    CHECK(f({-1.0, 2.0}) == 0.0);
    CHECK_NEAR(Expr("R * 2", {"r"}, none, u)({3.0}), 6.0, 1e-15);  // names are case-insensitive
    auto bad = [&](const std::string& t, const std::string& part) {
      ++checks;
      try {
        Expr(t, {"x"}, none, u);
        ++failures;
        std::printf("FAILED: '%s' compiled\n", t.c_str());
      } catch (const std::invalid_argument& e) {
        if (std::string(e.what()).find(part) == std::string::npos) {
          ++failures;
          std::printf("FAILED: '%s': message '%s' lacks '%s'\n", t.c_str(), e.what(), part.c_str());
        }
      }
    };
    bad("x + y", "unknown name 'y'");
    bad("foo(x)", "unknown function 'foo'");
    bad("(x + 1", "')' expected");
    bad("min(x)", "',' expected");
    bad("x +", "unexpected end");
    bad("x 3", "unexpected '3'");
    bad("2 parsec", "unexpected 'p'");

    // deck: block constants, expressions as values, profile functions
    const std::string G = "grid { lambda0 = 0.8 um; radial_grid = uniform; dr = 0.5 um; nr = 40; dzeta = 0.05 um;"
                          " window_length = 40 um; dz = 10 um; z_end = 1 mm; m_max = 1 }\n";
    Config c = cfgOf("constants { w0 = 8 um; L = 2 * w0; nb = 3e18 cm^-3 }\n" + G +
                     "laser { a0 = 1; waist = w0; length_1e = \"L / 4\"; zeta_center = 20 um }\n"
                     "plasma { density = nb; density_function = \"(1 + 0.5 * cos(z / w0)) * (x > 0 ? 1 : 0.5)\" }\n"
                     "bunch { charge = 1 pC; distribution_z = function; profile_z = \"rect(zeta - 10 um, 4 um)\";"
                     " zeta_min = 7 um; zeta_max = 13 um; gamma_function = \"100 + (zeta - 10 um) / (1 um)\";"
                     " sigma_x = 2 um; macroparticles = 10 }\n");
    CHECK_NEAR(c.constants.at("w0"), 8.0 * k0, 1e-14);
    CHECK_NEAR(c.constants.at("l"), 16.0 * k0, 1e-14);
    CHECK_NEAR(c.lasers[0].waist, 8.0 * k0, 1e-14);
    CHECK_NEAR(c.lasers[0].length_1e, 4.0 * k0, 1e-14);
    const double nb = c.plasma.n0;
    CHECK_NEAR(nb * c.units.nc_cm3(), 3e18, 1e-12);
    CHECK_NEAR(c.plasma.density(1.0, 0.0, 0.0), 1.5 * nb, 1e-14);
    CHECK_NEAR(c.plasma.density(-1.0, 0.0, 0.0), 0.75 * nb, 1e-14);
    CHECK_NEAR(c.plasma.density(1.0, 0.0, PI * 8.0 * k0), 0.5 * nb, 1e-12);
    const auto& B = c.bunches[0];
    CHECK(B.shape == BunchConfig::Shape::Function);
    CHECK_NEAR(B.zeta_center, 10.0 * k0, 1e-14);
    CHECK_NEAR(B.gamma, 100.0, 1e-3);  // mean of gamma(zeta) over the symmetric profile
    CHECK(B.profile_z({9.0 * k0}) == 1.0 && B.profile_z({12.5 * k0}) == 0.0);
    // the laser envelope and profile as functions
    Config l = cfgOf(G + "laser { a0 = 1; profile_function = \"exp(-(r / (5 um))^4)\";"
                         " envelope_function = \"sin(pi * (zeta - 10 um) / (20 um))^2 * rect(zeta - 20 um, 20 um)\" }");
    CHECK_NEAR(l.lasers[0].zeta_center, 20.0 * k0, 1e-6);
    CHECK_NEAR(l.lasers[0].waist, 5.0 * k0, 0.02);  // 1/e radius of the profile (tabulated)
    CHECK(!l.lasers[0].envelope_fn.empty() && !l.lasers[0].profile_fn.empty());
    EXPECT_ERROR(cfgOf(G + "laser { a0 = 1; profile_function = \"1\"; z_focus = 1 mm; length_1e = 3 um; zeta_center = 20 um }"),
                 "no focus position");
    EXPECT_ERROR(cfgOf(G + "laser { a0 = 1; waist = 3 um; envelope_function = \"1\"; length_1e = 3 um }"),
                 "give only one");
    EXPECT_ERROR(cfgOf(G + "plasma { density = 1e18 cm^-3; density_function = \"exp(-q)\" }"), "unknown name 'q'");
    EXPECT_ERROR(cfgOf(G + "plasma { density = 1e18 cm^-3; density_function = gauss(z, 3) }"), "in quotes");
    EXPECT_ERROR(cfgOf("constants { pi = 3 }\n" + G), "reserved");
    EXPECT_ERROR(cfgOf(G + "bunch { charge = 1 pC; gamma = 10; distribution_z = function; profile_z = \"-1\";"
                           " zeta_min = 7 um; zeta_max = 13 um; sigma_x = 2 um }"),
                 ">= 0");
    // undulator taper K(z)
    const std::string GU = "grid { lambda0 = 0.8 um; radial_grid = uniform; dr = 10 um; nr = 50; dzeta = 25 nm;"
                           " nzeta = 100; dz = 0.5 mm; z_end = 0.1 m }\n";
    Config t = cfgOf(GU + "undulator { period = 1 cm; K_function = \"1 - 0.1 * z / (10 cm)\"; periods = 10;"
                          " ramp_periods = 0; z_start = 0 cm }");
    const auto& U = t.undulators[0];
    CHECK_NEAR(U.K, 1.0, 1e-14);
    CHECK_NEAR(U.Kat(5e4 * k0), 0.95, 1e-12);
    EXPECT_ERROR(cfgOf(GU + "undulator { period = 1 cm; K = 1; K_function = \"1\"; periods = 10 }"), "not both");
    // the field from the taper table (as the simulation builds it)
    std::vector<double> tab;
    const int nk = 321;
    const double dzk = U.length / (nk - 1);
    for (int j = 0; j < nk; ++j) tab.push_back(U.Kat(U.z_start + j * dzk) * U.ku());
    UndulatorV v{0, 0, U.ku(), U.B0(), U.z_start, U.z_start + U.length, 0.0, nk, 0, dzk};
    ExternalV X{1, &v, tab.data()};
    EB fz{};
    addExternal(X, 0.0, 0.0, 5.0 * U.period, fz);  // z = 5 cm: K = 0.95, cos(10 pi) = 1
    CHECK_NEAR(fz.By, 0.95 * U.ku(), 1e-9);
  }

  // ------------------------------------------------------------ shot noise, beams from files
  {
    const std::string Bk = "bunch { charge = 1 pC; gamma = 100; sigma_x = 2 um; sigma_z = 1 um; zeta_center = 20 um;";
    Config c = cfgOf(std::string(BASE) + Bk + " macroparticles = 400; shot_noise = fawley; noise_wavelength = 0.8 um }");
    CHECK(c.bunches[0].fawley && c.bunches[0].beamlet == 4);
    CHECK_NEAR(c.bunches[0].noise_wavelength, 2.0 * PI, 1e-14);
    EXPECT_ERROR(cfgOf(std::string(BASE) + Bk + " macroparticles = 401; shot_noise = fawley; noise_wavelength = 0.8 um }"),
                 "multiple of beamlet");
    EXPECT_ERROR(cfgOf(std::string(BASE) + Bk + " shot_noise = fawley }"), "noise_wavelength");
    EXPECT_ERROR(cfgOf(std::string(BASE) + Bk + " noise_wavelength = 1 um }"), "shot_noise = fawley");
#if defined(GEMPIC_HAVE_HDF5)
    EXPECT_ERROR(cfgOf(std::string(BASE) + "bunch { file = /nonexistent.h5; zeta_center = 10 um }"), "cannot open");
#else
    EXPECT_ERROR(cfgOf(std::string(BASE) + "bunch { file = /nonexistent.h5; zeta_center = 10 um }"), "needs HDF5");
#endif
  }

  // ------------------------------------------------------------ transport elements
  {
    Config c = cfgOf(std::string(BASE) + "quadrupole { z_start = 1 cm; length = 5 cm; gradient = 100 T/m }\n"
                     "solenoid { z_start = 0 m; length = 10 cm; field = 2 T; ramp = 1 cm }\n"
                     "chicane { z = 0.5 mm; r56 = 50 um; energy = 1 GeV }\nchicane { z = 0.2 mm; r56 = -1 um; gamma = 100 }");
    CHECK(c.magnets.size() == 2 && c.chicanes.size() == 2);
    const double L = 0.8e-6 / (2.0 * PI);  // code length [m]
    CHECK_NEAR(c.magnets[0].gradient * c.units.B0_T() / L, 100.0, 1e-12);
    CHECK_NEAR(c.magnets[1].field * c.units.B0_T(), 2.0, 1e-12);
    CHECK_NEAR(c.chicanes[0].z * L, 0.2e-3, 1e-12);          // sorted by z
    CHECK_NEAR(c.chicanes[1].delay * L, 25e-6, 1e-12);       // default |R56| / 2
    // fields: quadrupole focusing for electrons, solenoid div B = 0 in the fringe
    MagnetV q{0, 0.0, 1.0, 2.0, 1e300, 0.0}, so{2, 0.0, 10.0, 3.0, 1e300, 2.0};
    ExternalV X;
    X.nm = 1; X.m = &q;
    EB f{};
    addExternal(X, 0.5, 0.25, 0.5, f);
    CHECK_NEAR(f.By, -1.0, 1e-15);
    CHECK_NEAR(f.Bx, -0.5, 1e-15);
    X.m = &so;
    EB a{}, b{};
    const double h = 1e-4;
    addExternal(X, 0.3, 0.2, 1.0 - h, a);
    addExternal(X, 0.3, 0.2, 1.0 + h, b);
    EB e{};
    addExternal(X, 0.3, 0.2, 1.0, e);
    CHECK(std::fabs(e.Bx / 0.3 + e.By / 0.2 + (b.Bz - a.Bz) / (2 * h)) < 1e-9);  // div B = 0 in the fringe
    CHECK(std::fabs(e.Bz - 1.5) < 1e-12);
    EXPECT_ERROR(cfgOf(std::string(BASE) + "solenoid { z_start = 0 m; length = 1 cm; field = 2 T; ramp = 1 cm }"), "fit twice");
    EXPECT_ERROR(cfgOf(std::string(BASE) + "quadrupole { z_start = 0 m; length = 1 cm; gradient = 2 T }"), "magnetic field");
  }

  // ------------------------------------------------------------ output formats
#if defined(GEMPIC_HAVE_HDF5)
  {
    Config o = cfgOf(std::string(BASE) + "output { format = text, openpmd; openpmd_dr = 1 um }");
    CHECK(o.output.text && o.output.openpmd && !o.output.hdf5);
    CHECK_NEAR(o.output.openpmd_dr, 2.0 * PI / 0.8, 1e-14);
    Config o2 = cfgOf(std::string(BASE) + "output { format = both }");
    CHECK(o2.output.text && o2.output.hdf5 && !o2.output.openpmd);
    EXPECT_ERROR(cfgOf(std::string(BASE) + "output { format = vtk }"), "'vtk' is not allowed");
    EXPECT_ERROR(cfgOf(std::string(BASE) + "output { openpmd_dr = 1 um }"), "needs format = openpmd");
  }
#endif

  // ------------------------------------------------------------ example decks must be valid
#ifdef GEM_SOURCE_DIR
  for (const auto& f : std::filesystem::directory_iterator(std::string(GEM_SOURCE_DIR) + "/examples")) {
    if (f.path().extension() != ".deck") continue;
    ++checks;
    try {
      Config::fromDeck(Deck::parseFile(f.path().string()));
    } catch (const std::exception& e) {
      ++failures;
      std::printf("FAILED: example %s: %s\n", f.path().filename().c_str(), e.what());
    }
  }
  ++checks;
  try {
    Config::fromDeck(Deck::parseFile(std::string(GEM_SOURCE_DIR) + "/tests/regression.deck"));
  } catch (const std::exception& e) {
    ++failures;
    std::printf("FAILED: regression.deck: %s\n", e.what());
  }
#endif

  std::printf("%d checks, %d failed\n", checks, failures);
  return failures;
}
