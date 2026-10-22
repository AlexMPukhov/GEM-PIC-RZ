// ============================================================================
//  test_deck_physics.cpp -- deck and unit tests: tridiagonal solver, undulator, NIST data, ionizable ions, overrides, function parser
// ============================================================================
#include "test_deck_common.hpp"

void test_deck_physics([[maybe_unused]] const Units& u, [[maybe_unused]] double k0) {
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
}
