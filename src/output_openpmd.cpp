// ============================================================================
//  output_openpmd.cpp -- openPMD output (output.format = openpmd | both), see src/output_hdf5.cpp
// ============================================================================
#include "output_hdf5_internal.hpp"

#if defined(GEMPIC_HAVE_HDF5)

namespace gem {
using namespace h5detail;

// ============================================================================
//  openPMD 1.1.0 output (output.format = openpmd): <dir>/openpmd/openpmd_<step>.h5
//
//  fileBased iteration encoding, basePath /data/%T/. Readable with
//  openPMD-viewer, openPMD-api, VisualPIC, ...
//  meshes   E, B (components r, t, z), n_e (and Zn_i with ion blocks):
//           geometry thetaMode, dataset [2M+1][Nr][Nz] (C order, axes r, z):
//           index 0 the mode m = 0, index 2m-1 / 2m the coefficients of
//           cos(m theta) / sin(m theta) ("imag=+"):
//               F = F_0 + sum_m ( C_m cos(m theta) + S_m sin(m theta) ),
//           C_m = 2 Re F_m, S_m = -2 Im F_m for the code's modes F_m e^{i m theta}.
//           The fields are interpolated (as in the particle gather) from the
//           staggered, possibly non-uniform r-zeta grid to the uniform nodes
//           r_j = j dr, zeta_k = k dzeta (output.openpmd_dr, openpmd_dzeta).
//           Lab coordinate z = zeta + s, time t = s / c (gridGlobalOffset[1] = s).
//  particles  one species per bunch (its name), "trapped" (self-injected plasma
//           electrons) and "inj_<ion>" (electrons from ionization): position
//           x, y, z (z: zeta, positionOffset z = s), momentum (per real
//           particle), weighting (real particles per macroparticle), charge and
//           mass (constant records).
// ============================================================================
void Simulation::writeOpenPMD() {
  namespace fs = std::filesystem;
  const auto& O = cfg.output;
  const std::string dir = O.dir + "/openpmd";
  fs::create_directories(dir);
  const std::string name = dir + "/openpmd_" + std::to_string(step_no) + ".h5";
  H5Eset_auto2(H5E_DEFAULT, nullptr, nullptr);
  Writer W(name, O.hdf5_single, O.hdf5_compression, true);

  const Units& u = cfg.units;
  const double L = u.k0inv_um() * 1e-6;                // length unit [m]
  const double E0 = MC2_MEV * 1e6 / L;                 // [V/m]
  const double B0 = E0 / C_LIGHT;                      // [T]
  const double T0 = 1.0 / (u.omega0_per_fs() * 1e15);  // [s]
  const double nc = u.nc_cm3() * 1e6;                  // [m^-3]
  const double Q0 = u.codeToPC(1.0) * 1e-12;           // [C]
  const double ME = 9.1093837015e-31;                  // [kg]
  const int N = g.N;

  // ---- root
  const hid_t root = W.file();
  Writer::attr(root, "openPMD", std::string("1.1.0"));
  Writer::attrU32(root, "openPMDextension", 0u);
  Writer::attr(root, "basePath", std::string("/data/%T/"));
  Writer::attr(root, "meshesPath", std::string("meshes/"));
  Writer::attr(root, "particlesPath", std::string("particles/"));
  Writer::attr(root, "iterationEncoding", std::string("fileBased"));
  Writer::attr(root, "iterationFormat", std::string("openpmd_%T.h5"));
  Writer::attr(root, "software", std::string("GEM-PIC-RZ"));
  Writer::attr(root, "softwareVersion", std::string("1"));
  Writer::attr(root, "author", std::string("GEM-PIC-RZ user"));
  {
    char buf[64];
    const std::time_t t = std::time(nullptr);
    std::strftime(buf, sizeof buf, "%Y-%m-%d %H:%M:%S %z", std::localtime(&t));
    Writer::attr(root, "date", std::string(buf));
  }
  Writer::attr(root, "comment", std::string("GEM-PIC-RZ (Galilean GT2 scheme): fields of the window at z = zeta + s, "
                                            "t = s / c; input deck in /data/<T>/input_deck"));
  const std::string base = "/data/" + std::to_string(step_no);
  {
    Handle it(W.groups(base), H5Gclose, "group " + base);
    Writer::attr(it, "time", s0);
    Writer::attr(it, "dt", Delta);
    Writer::attr(it, "timeUnitSI", T0);
  }
  W.text(base + "/input_deck", deckText);
  auto unitDim = [](double Lx, double Mx, double Tx, double Ix) { return std::vector<double>{Lx, Mx, Tx, Ix, 0, 0, 0}; };

  // ---- output grid
  double dr = O.openpmd_dr > 0 ? O.openpmd_dr : 1e300, dzo = O.openpmd_dzeta > 0 ? O.openpmd_dzeta : 1e300;
  if (!(O.openpmd_dr > 0)) for (int i = 0; i < N; ++i) dr = std::min(dr, g.rho[i + 1] - g.rho[i]);
  if (!(O.openpmd_dzeta > 0)) for (int n = 0; n < Nz; ++n) dzo = std::min(dzo, zv.len(n));
  const int nr = int(std::floor(g.rho[N] / dr + 1e-9)) + 1;
  const int nz = int(std::floor(zv.node(Nz) / dzo + 1e-9)) + 1;
  const int NM = 2 * M + 1;

  // radial stencils (the gather of the code, see tools/h5_extract.cpp)
  struct Rst { int modeC, ac, modeF, af; double wc, wf, xr; };
  std::vector<Rst> rs(nr);
  for (int j = 0; j < nr; ++j) {
    const double r = j * dr;
    Rst& t = rs[j];
    t = Rst{0, 0, 0, 0, 0.0, 0.0, r / g.c[0]};
    if (r < g.c[0]) t.modeC = 1;
    else if (r >= g.c[N - 1]) { t.modeC = 2; t.ac = N - 1; }
    else {
      t.ac = int(std::upper_bound(g.c.begin(), g.c.begin() + N, r) - g.c.begin()) - 1;
      t.wc = (r - g.c[t.ac]) / (g.c[t.ac + 1] - g.c[t.ac]);
    }
    if (r >= g.rho[N]) { t.modeF = 2; t.af = N; }
    else {
      t.af = int(std::upper_bound(g.rho.begin(), g.rho.begin() + N + 1, r) - g.rho.begin()) - 1;
      t.wf = (r - g.rho[t.af]) / (g.rho[t.af + 1] - g.rho[t.af]);
    }
  }
  // zeta stencils: slice n and weight
  std::vector<int> zn(nz);
  std::vector<double> zw(nz);
  for (int k = 0; k < nz; ++k) {
    const double z = k * dzo;
    double wdum = 0.0;
    int n = std::max(0, std::min(Nz - 1, zv.cell(z, wdum)));
    while (n > 0 && zv.node(n) > z) --n;
    while (n < Nz - 1 && zv.node(n + 1) <= z) ++n;
    zn[k] = n;
    zw[k] = std::min(1.0, std::max(0.0, (z - zv.node(n)) / zv.len(n)));
  }

  // complex mode m of quantity q at (slice n, stencil t); q: 0 Er 1 Et 2 Ez 3 Br 4 Bt 5 Bz
  const cd I(0.0, 1.0);
  auto Ephi = [&](int m, int n, int i) { return 0.5 * (fold.Rp[fold.c(m, n, i)] + fold.Lp[fold.c(m, n, i)]); };
  auto Brad = [&](int m, int n, int i) { return 0.5 * (fold.Lp[fold.c(m, n, i)] - fold.Rp[fold.c(m, n, i)]); };
  auto Erad = [&](int m, int n, int i) { return 0.5 * (fold.Rr[fold.f(m, n, i)] + fold.Lr[fold.f(m, n, i)]); };
  auto Bphi = [&](int m, int n, int i) { return 0.5 * (fold.Rr[fold.f(m, n, i)] - fold.Lr[fold.f(m, n, i)]); };
  auto Ezc = [&](int m, int n, int i) { return fold.Ez[fold.c(m, n, i)]; };
  auto Bzf = [&](int m, int n, int i) { return fold.Bz[fold.f(m, n, i)]; };
  auto lin0 = [&](cd f0, cd f1) { return (g.c[1] * f0 - g.c[0] * f1) / (g.c[1] - g.c[0]); };  // extrapolation to r = 0
  auto mode = [&](int q, int m, int n, const Rst& t) -> cd {
    if (q == 1 || q == 2 || q == 3) {  // centres
      auto F = [&](int i) { return q == 1 ? Ephi(m, n, i) : q == 2 ? Ezc(m, n, i) : Brad(m, n, i); };
      if (t.modeC == 0) return (1 - t.wc) * F(t.ac) + t.wc * F(t.ac + 1);
      if (t.modeC == 2) return F(t.ac);
      if (q == 2) return F(0) * ((m == 0) ? 1.0 : std::pow(t.xr, m));
      if (m == 0) return F(0) * t.xr;
      if (m == 1) { const cd a0 = (N > 1) ? lin0(F(0), F(1)) : F(0); return a0 + (F(0) - a0) * t.xr; }
      return F(0) * std::pow(t.xr, m - 1);
    }
    // faces; on the axis E_r = -i E_phi, B_phi = i B_r for m = 1, else 0
    auto F = [&](int i) -> cd {
      if (q == 5) return Bzf(m, n, i);
      if (i == 0) {
        if (m != 1 || N < 2) return cd(0.0);
        return q == 0 ? -I * lin0(Ephi(m, n, 0), Ephi(m, n, 1)) : I * lin0(Brad(m, n, 0), Brad(m, n, 1));
      }
      return q == 0 ? Erad(m, n, i) : Bphi(m, n, i);
    };
    if (t.modeF == 0) return (1 - t.wf) * F(t.af) + t.wf * F(t.af + 1);
    return F(t.af);
  };

  const std::vector<std::string> axes = {"r", "z"};
  auto meshAttrs = [&](hid_t obj, int modes, const std::vector<double>& dim) {
    Writer::attr(obj, "geometry", std::string("thetaMode"));
    Writer::attr(obj, "geometryParameters", "m=" + std::to_string(modes) + ";imag=+");
    Writer::attr(obj, "dataOrder", std::string("C"));
    Writer::attr(obj, "axisLabels", axes);
    Writer::attr(obj, "gridSpacing", std::vector<double>{dr, dzo});
    Writer::attr(obj, "gridGlobalOffset", std::vector<double>{0.0, s0});
    Writer::attr(obj, "gridUnitSI", L);
    Writer::attr(obj, "unitDimension", dim);
    Writer::attr(obj, "timeOffset", 0.0);
  };
  const std::string mp = base + "/meshes";
  { Handle h(W.groups(mp), H5Gclose, "group " + mp); }
  if (O.openpmd_uniform) {
    std::vector<double> out(size_t(NM) * nr * nz);
    std::vector<cd> slab(size_t(Nz + 1) * nr);  // one mode of one quantity on all slices, output radii
    const char* comp[3] = {"r", "t", "z"};
    for (int vec = 0; vec < 2; ++vec) {  // E, B
      const std::string rec = mp + (vec == 0 ? "/E" : "/B");
      {
        Handle h(W.groups(rec), H5Gclose, "group " + rec);
        meshAttrs(h, M + 1, vec == 0 ? unitDim(1, 1, -3, -1) : unitDim(0, 1, -2, -1));
      }
      for (int cpt = 0; cpt < 3; ++cpt) {
        const int q = 3 * vec + cpt;
        for (int m = 0; m <= M; ++m) {
          for (int n = 0; n <= Nz; ++n)
            for (int j = 0; j < nr; ++j) slab[size_t(n) * nr + j] = mode(q, m, n, rs[j]);
          for (int j = 0; j < nr; ++j)
            for (int k = 0; k < nz; ++k) {
              const cd v = (1.0 - zw[k]) * slab[size_t(zn[k]) * nr + j] + zw[k] * slab[size_t(zn[k] + 1) * nr + j];
              if (m == 0) {
                out[(size_t(0) * nr + j) * nz + k] = v.real();
              } else {
                out[(size_t(2 * m - 1) * nr + j) * nz + k] = 2.0 * v.real();
                out[(size_t(2 * m) * nr + j) * nz + k] = -2.0 * v.imag();
              }
            }
        }
        Handle ds(W.plain(rec + "/" + comp[cpt], {hsize_t(NM), hsize_t(nr), hsize_t(nz)}, out.data()), H5Dclose, "E/B");
        Writer::attr(ds, "unitSI", vec == 0 ? E0 : B0);
        Writer::attr(ds, "position", std::vector<double>{0.0, 0.0});
      }
    }
    // densities (m = 0, on half-cells x centres): electrons and the charge of the mobile ions
    auto density = [&](const std::string& nm, bool ions) {
      std::vector<double> d(size_t(nr) * nz);
      std::vector<double> zc(Nz);
      for (int kk = 0; kk < Nz; ++kk) zc[kk] = zv.centre(kk);
      for (int j = 0; j < nr; ++j) {
        const double r = j * dr;
        int a = 0;
        double w = 0.0;
        if (r >= g.c[N - 1]) { a = N - 2; w = 1.0; }
        else if (r > g.c[0]) {
          a = int(std::upper_bound(g.c.begin(), g.c.begin() + N, r) - g.c.begin()) - 1;
          w = (r - g.c[a]) / (g.c[a + 1] - g.c[a]);
        }
        if (N < 2) { a = 0; w = 0.0; }
        for (int k = 0; k < nz; ++k) {
          const double z = k * dzo;
          int b = int(std::upper_bound(zc.begin(), zc.end(), z) - zc.begin()) - 1;
          double v;
          auto val = [&](int kk, int i) {
            const size_t idx = size_t(kk) * N + i;
            return ions ? (c0s.rhoi[idx] + c0b.rhoi[idx]) : -(c0s.rho[idx] + c0b.rho[idx]);
          };
          auto rad = [&](int kk) { return N < 2 ? val(kk, 0) : (1.0 - w) * val(kk, a) + w * val(kk, a + 1); };
          if (b < 0) v = rad(0);
          else if (b >= Nz - 1) v = rad(Nz - 1);
          else {
            const double wz = (z - zc[b]) / (zc[b + 1] - zc[b]);
            v = (1.0 - wz) * rad(b) + wz * rad(b + 1);
          }
          d[size_t(j) * nz + k] = v;
        }
      }
      Handle ds(W.plain(mp + "/" + nm, {1, hsize_t(nr), hsize_t(nz)}, d.data()), H5Dclose, nm);
      meshAttrs(ds, 1, unitDim(-3, 0, 0, 0));
      Writer::attr(ds, "unitSI", nc);
      Writer::attr(ds, "position", std::vector<double>{0.0, 0.0});
    };
    density("n_e", false);
    if (!cfg.ions.empty()) density("Zn_i", true);
  }

  // ---- native records (output.openpmd_grid = native | both): every mode on the code's own staggered,
  // possibly non-uniform r-zeta nodes, without interpolation. Scalar meshes "native_<name>"; the node
  // positions are data sets <basePath>gridCoordinates/{r_faces, r_centers, zeta_slices, zeta_halfcells},
  // referenced by the mesh attribute gridCoordinates (one path per axis), as in the proposal for
  // rectilinear meshes in openPMD 2.0 (docs/openpmd_rectilinear.md). geometry "other" (intended geometry
  // in geometryParameters) while the files are openPMD 1.1, so that readers of regular grids skip them. Fields in m c omega0 / e and
  // m omega0 / e, densities in units of n0 (unitSI converts to SI). E_r and B_phi on the axis (face 0):
  // the regular limit, -i E_phi and i B_r of mode 1 extrapolated to r = 0, else 0.
  if (O.openpmd_native) {
    const std::vector<double> rf(g.rho.begin(), g.rho.begin() + N + 1), rcn(g.c.begin(), g.c.begin() + N);
    std::vector<double> zsl(Nz + 1), zh(Nz);
    for (int n = 0; n <= Nz; ++n) zsl[n] = zv.node(n);
    for (int k = 0; k < Nz; ++k) zh[k] = zv.centre(k);
    auto minStep = [](const std::vector<double>& x) {
      double d = 1e300;
      for (size_t k = 1; k < x.size(); ++k) d = std::min(d, x[k] - x[k - 1]);
      return x.size() > 1 ? d : 1.0;
    };
    const std::string gc = base + "/gridCoordinates";
    { Handle h(W.groups(gc), H5Gclose, gc); }
    auto coords = [&](const char* nm, const std::vector<double>& x) {
      Handle h(W.plain(gc + "/" + nm, {hsize_t(x.size())}, x.data(), true), H5Dclose, nm);
      return gc + "/" + nm;
    };
    const std::string pRf = coords("r_faces", rf), pRc = coords("r_centers", rcn);
    const std::string pZs = coords("zeta_slices", zsl), pZh = coords("zeta_halfcells", zh);
    auto nativeMesh = [&](const std::string& nm, int modes, const std::vector<double>& rr,
                          const std::vector<double>& zz, const std::string& pr, const std::string& pz,
                          const std::vector<double>& data, const std::vector<double>& dim, double unitSI) {
      Handle ds(W.plain(mp + "/" + nm, {hsize_t(2 * modes - 1), hsize_t(rr.size()), hsize_t(zz.size())}, data.data()),
                H5Dclose, nm);
      Writer::attr(ds, "geometry", std::string("other"));
      Writer::attr(ds, "geometryParameters", "thetaMode;m=" + std::to_string(modes) + ";imag=+");
      Writer::attr(ds, "gridType", std::string("rectilinear"));
      Writer::attr(ds, "dataOrder", std::string("C"));
      Writer::attr(ds, "axisLabels", axes);
      Writer::attr(ds, "gridSpacing", std::vector<double>{minStep(rr), minStep(zz)});
      Writer::attr(ds, "gridGlobalOffset", std::vector<double>{0.0, s0});
      Writer::attr(ds, "gridUnitSI", L);
      Writer::attr(ds, "gridCoordinates", std::vector<std::string>{pr, pz});
      Writer::attr(ds, "unitDimension", dim);
      Writer::attr(ds, "timeOffset", 0.0);
      Writer::attr(ds, "position", std::vector<double>{0.0, 0.0});
      Writer::attr(ds, "unitSI", unitSI);
    };
    // q: 0 E_r 1 E_phi 2 E_z 3 B_r 4 B_phi 5 B_z
    auto native = [&](int q, int m, int n, int i) -> cd {
      switch (q) {
        case 0: return i > 0 ? Erad(m, n, i) : (m == 1 && N > 1 ? -I * lin0(Ephi(m, n, 0), Ephi(m, n, 1)) : cd(0.0));
        case 1: return Ephi(m, n, i);
        case 2: return Ezc(m, n, i);
        case 3: return Brad(m, n, i);
        case 4: return i > 0 ? Bphi(m, n, i) : (m == 1 && N > 1 ? I * lin0(Brad(m, n, 0), Brad(m, n, 1)) : cd(0.0));
        default: return Bzf(m, n, i);
      }
    };
    const char* nname[6] = {"native_E_r", "native_E_t", "native_E_z", "native_B_r", "native_B_t", "native_B_z"};
    const bool onFaces[6] = {true, false, false, false, true, true};
    for (int q = 0; q < 6; ++q) {
      const std::vector<double>& rr = onFaces[q] ? rf : rcn;
      const size_t NR = rr.size(), NZ = size_t(Nz) + 1;
      std::vector<double> d(size_t(NM) * NR * NZ);
      for (int m = 0; m <= M; ++m)
        for (size_t i = 0; i < NR; ++i)
          for (size_t n = 0; n < NZ; ++n) {
            const cd v = native(q, m, int(n), int(i));
            if (m == 0) d[i * NZ + n] = v.real();
            else {
              d[(size_t(2 * m - 1) * NR + i) * NZ + n] = 2.0 * v.real();
              d[(size_t(2 * m) * NR + i) * NZ + n] = -2.0 * v.imag();
            }
          }
      nativeMesh(nname[q], M + 1, rr, zsl, onFaces[q] ? pRf : pRc, pZs, d,
                 q < 3 ? unitDim(1, 1, -3, -1) : unitDim(0, 1, -2, -1), q < 3 ? E0 : B0);
    }
    // densities, mode 0, on the radial centres and the half-cells
    std::vector<double> d(size_t(N) * Nz);
    for (int i = 0; i < N; ++i)
      for (int k = 0; k < Nz; ++k) {
        const size_t idx = size_t(k) * N + i;
        d[size_t(i) * Nz + k] = -(c0s.rho[idx] + c0b.rho[idx]) / n_ref;
      }
    nativeMesh("native_n_e", 1, rcn, zh, pRc, pZh, d, unitDim(-3, 0, 0, 0), n_ref * nc);
    if (!cfg.ions.empty()) {
      for (int i = 0; i < N; ++i)
        for (int k = 0; k < Nz; ++k) {
          const size_t idx = size_t(k) * N + i;
          d[size_t(i) * Nz + k] = (c0s.rhoi[idx] + c0b.rhoi[idx]) / n_ref;
        }
      nativeMesh("native_Zn_i", 1, rcn, zh, pRc, pZh, d, unitDim(-3, 0, 0, 0), n_ref * nc);
    }
  }

  // ---- particles
  // particlesPath only when the iteration has particles (the standard requires the path to exist)
  if (!O.particles) { H5Adelete(root, "particlesPath"); return; }
  struct Sp { std::string name; int id; };
  std::vector<Sp> species;
  for (size_t b = 0; b < cfg.bunches.size(); ++b) species.push_back({cfg.bunches[b].name, int(b)});
  species.push_back({"trapped", SPECIES_TRAPPED});
  for (size_t s = 0; s < cfg.ions.size(); ++s)
    species.push_back({"inj_" + cfg.ions[s].name, SPECIES_TRAPPED_IONIZED_BASE - int(s)});
  const std::string pp = base + "/particles";
  bool first = true;
  for (const auto& S : species) {
    std::vector<const Particle*> ps;
    for (const auto& p : beam)
      if (p.sp == S.id) ps.push_back(&p);
    if (ps.empty() && S.id < 0) continue;  // no trapped particles: no species
    if (first) { Handle h(W.groups(pp), H5Gclose, "group " + pp); first = false; }
    const std::string sg = pp + "/" + S.name;
    { Handle h(W.groups(sg), H5Gclose, "group " + sg); }
    const hsize_t np = ps.size();
    std::vector<double> col(np);
    auto record = [&](const std::string& rec, const std::vector<double>& dim, unsigned macroWeighted, double wPower) {
      Handle h(W.groups(sg + "/" + rec), H5Gclose, rec);
      Writer::attr(h, "unitDimension", dim);
      Writer::attr(h, "timeOffset", 0.0);
      Writer::attrU32(h, "macroWeighted", macroWeighted);
      Writer::attr(h, "weightingPower", wPower);
    };
    auto component = [&](const std::string& path, double (*get)(const Particle&), double unitSI) {
      for (size_t k = 0; k < np; ++k) col[k] = get(*ps[k]);
      Handle ds(W.plain(sg + "/" + path, {np}, col.data()), H5Dclose, path);
      Writer::attr(ds, "unitSI", unitSI);
    };
    auto constant = [&](const std::string& path, double value, double unitSI) {
      Handle h(W.groups(sg + "/" + path), H5Gclose, path);
      Writer::attr(h, "value", value);
      Writer::attr(h, "shape", std::vector<unsigned long long>{np});
      Writer::attr(h, "unitSI", unitSI);
    };
    record("position", unitDim(1, 0, 0, 0), 0u, 0.0);
    component("position/x", [](const Particle& p) { return p.x; }, L);
    component("position/y", [](const Particle& p) { return p.y; }, L);
    component("position/z", [](const Particle& p) { return p.zeta; }, L);
    record("positionOffset", unitDim(1, 0, 0, 0), 0u, 0.0);
    constant("positionOffset/x", 0.0, L);
    constant("positionOffset/y", 0.0, L);
    constant("positionOffset/z", s0, L);
    record("momentum", unitDim(1, 1, -1, 0), 0u, 1.0);
    component("momentum/x", [](const Particle& p) { return p.ux; }, ME * C_LIGHT);
    component("momentum/y", [](const Particle& p) { return p.uy; }, ME * C_LIGHT);
    component("momentum/z", [](const Particle& p) { return p.uz; }, ME * C_LIGHT);
    // real particles per macroparticle: |q| Q0 / e (all beam particles are electrons or positrons)
    {
      for (size_t k = 0; k < np; ++k) col[k] = std::fabs(ps[k]->q) * Q0 / E_CHARGE_C;
      Handle ds(W.plain(sg + "/weighting", {np}, col.data()), H5Dclose, "weighting");
      Writer::attr(ds, "unitDimension", unitDim(0, 0, 0, 0));
      Writer::attr(ds, "timeOffset", 0.0);
      Writer::attrU32(ds, "macroWeighted", 1u);
      Writer::attr(ds, "weightingPower", 1.0);
      Writer::attr(ds, "unitSI", 1.0);
    }
    const double qsign = (S.id >= 0) ? (cfg.bunches[size_t(S.id)].qm < 0 ? -1.0 : 1.0) : -1.0;
    {
      Handle h(W.groups(sg + "/charge"), H5Gclose, "charge");
      Writer::attr(h, "unitDimension", unitDim(0, 0, 1, 1));
      Writer::attr(h, "timeOffset", 0.0);
      Writer::attrU32(h, "macroWeighted", 0u);
      Writer::attr(h, "weightingPower", 1.0);
      Writer::attr(h, "value", qsign * E_CHARGE_C);
      Writer::attr(h, "shape", std::vector<unsigned long long>{np});
      Writer::attr(h, "unitSI", 1.0);
    }
    // particlePatches: one patch with the bounding box (lab z = zeta + s)
    {
      double lo[3] = {0, 0, 0}, hi[3] = {0, 0, 0};
      for (size_t k = 0; k < np; ++k) {
        const double v[3] = {ps[k]->x, ps[k]->y, ps[k]->zeta + s0};
        for (int a = 0; a < 3; ++a) {
          if (k == 0 || v[a] < lo[a]) lo[a] = v[a];
          if (k == 0 || v[a] > hi[a]) hi[a] = v[a];
        }
      }
      const std::string pg = sg + "/particlePatches";
      { Handle h(W.groups(pg), H5Gclose, pg); }
      const unsigned long long cnt = np, off0 = 0;
      W.u64(pg + "/numParticles", {1}, &cnt, "particles in the patch");
      W.u64(pg + "/numParticlesOffset", {1}, &off0, "offset of the patch in the particle records");
      for (int which = 0; which < 2; ++which) {
        const std::string rec = pg + (which == 0 ? "/offset" : "/extent");
        {
          Handle h(W.groups(rec), H5Gclose, rec);
          Writer::attr(h, "unitDimension", unitDim(1, 0, 0, 0));
        }
        const char* ax[3] = {"x", "y", "z"};
        for (int a = 0; a < 3; ++a) {
          const double v = which == 0 ? lo[a] : hi[a] - lo[a];
          Handle ds(W.plain(rec + "/" + ax[a], {1}, &v), H5Dclose, rec);
          Writer::attr(ds, "unitSI", L);
        }
      }
    }
    {
      Handle h(W.groups(sg + "/mass"), H5Gclose, "mass");
      Writer::attr(h, "unitDimension", unitDim(0, 1, 0, 0));
      Writer::attr(h, "timeOffset", 0.0);
      Writer::attrU32(h, "macroWeighted", 0u);
      Writer::attr(h, "weightingPower", 1.0);
      Writer::attr(h, "value", ME);
      Writer::attr(h, "shape", std::vector<unsigned long long>{np});
      Writer::attr(h, "unitSI", 1.0);
    }
  }
  if (first) H5Adelete(root, "particlesPath");
}

}  // namespace gem

#endif
