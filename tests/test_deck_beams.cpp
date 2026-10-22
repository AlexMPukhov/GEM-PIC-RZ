// ============================================================================
//  test_deck_beams.cpp -- deck and unit tests: shot noise and beams, transport elements, output formats, example decks
// ============================================================================
#include "test_deck_common.hpp"

void test_deck_beams([[maybe_unused]] const Units& u, [[maybe_unused]] double k0) {
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
}
