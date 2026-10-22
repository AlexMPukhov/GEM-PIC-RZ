// ============================================================================
//  test_deck_basic.cpp -- deck and unit tests: syntax, units, minimal config, grids, lasers, plasma, errors, overrides, ions, deposition, backends
// ============================================================================
#include "test_deck_common.hpp"

void test_deck_basic([[maybe_unused]] const Units& u, [[maybe_unused]] double k0) {
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
}
