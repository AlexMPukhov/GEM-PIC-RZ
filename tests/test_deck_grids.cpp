// ============================================================================
//  test_deck_grids.cpp -- deck and unit tests: output format, free-hand mappings and zeta grids, particle pushers
// ============================================================================
#include "test_deck_common.hpp"

void test_deck_grids([[maybe_unused]] const Units& u, [[maybe_unused]] double k0) {
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
}
