// ============================================================================
//  output_hdf5.cpp -- simulation results in HDF5 (output.format = hdf5 | both)
//
//  One file per output step:  <output.dir>/data_NNNNN.h5  (NNNNN = step)
//
//  /                          attributes: code, format_version, step, z,
//                             lambda0_m, k0inv_m, omega0_per_s, n0_m3,
//                             E0_V_per_m, B0_T, modes, backend, ...
//  /grid/r_faces              [Nr+1]    rho_i   (faces; rho_0 = 0 is the axis)
//  /grid/r_centres            [Nr]      c_i     (centres)
//  /grid/zeta_slices          [Nz+1]    zeta_n (uniform: n * dzeta; slice Nz: ghost ahead of the window)
//  /grid/zeta_halfcells       [Nz]      centres of the half-cells (uniform: (k + 1/2) * dzeta)
//  /fields/Er, Bphi, Bz       [M+1][Nz+1][Nr+1]  complex, on faces
//  /fields/Ephi, Br, Ez       [M+1][Nz+1][Nr]    complex, on centres
//  /densities/ne_over_n0      [Nz][Nr]  electron density (plasma + bunches), m = 0, half-cells x centres
//  /densities/Zni_over_n0     [Nz][Nr]  charge density of the mobile ions / e n0
//  /particles/beam/x, y, z, zeta, ux, uy, uz, charge   [Np]   (output.particles = true)
//  /particles/beam/species    [Np]      int: >= 0 bunch index (attribute bunch_names),
//                                       -1 self-trapped plasma electrons, -30-s trapped
//                                       electrons from ionization of ion block s
//  /input_deck                string: the deck after command-line overrides
//
//  Complex numbers are stored as the compound {r, i} (the h5py convention; a
//  C++ reader can read them directly into std::complex<double>).
//  Mode convention:  F(r, phi) = F_0(r) + 2 Re sum_{m=1}^{M} F_m(r) exp(i m phi).
//  All data are in code units; every dataset carries the attribute unitSI
//  (multiply to get SI) and a human-readable unit and description.
//  Axis values: face 0 of E_r, B_phi holds 0; the fields on the axis follow
//  from regularity as in the particle gather (see tools/h5_extract.cpp).
// ============================================================================
#include "output_hdf5_internal.hpp"

#if defined(GEMPIC_HAVE_HDF5)

namespace gem {
using namespace h5detail;

// ============================================================================
void Simulation::writeHdf5() {
  namespace fs = std::filesystem;
  const auto& O = cfg.output;
  fs::create_directories(O.dir);
  char name[1024];
  std::snprintf(name, sizeof name, "%s/data_%05d.h5", O.dir.c_str(), step_no);
  // written as name.tmp and renamed when complete: a reader (e.g. tools/gempic_gui.py during a run)
  // never sees a partial data_NNNNN.h5
  const std::string tmp = std::string(name) + ".tmp";
  writeHdf5File(tmp);
  fs::rename(tmp, name);
}

void Simulation::writeHdf5File(const std::string& name) {
  const auto& O = cfg.output;
  H5Eset_auto2(H5E_DEFAULT, nullptr, nullptr);  // errors are reported by exceptions
  Writer W(name, O.hdf5_single, O.hdf5_compression);

  const Units& u = cfg.units;
  const double L = u.k0inv_um() * 1e-6;          // length unit [m]
  const double E0 = MC2_MEV * 1e6 / L;           // field unit m c omega0 / e [V/m]
  const double B0 = E0 / C_LIGHT;                // [T]
  const double nc = u.nc_cm3() * 1e6;            // [m^-3]
  const double Q0 = u.codeToPC(1.0) * 1e-12;     // charge unit [C]
  const int N = g.N;
  const std::string um = "1/k0 = " + std::to_string(u.k0inv_um()) + " um";

  // ---- root attributes
  const hid_t root = W.file();
  Writer::attr(root, "code", std::string("GEM-PIC-RZ"));
  Writer::attr(root, "format_version", 1LL);
  Writer::attr(root, "step", (long long)step_no);
  Writer::attr(root, "z", s0);
  Writer::attr(root, "z_m", s0 * L);
  Writer::attr(root, "lambda0_m", u.lambda0_um * 1e-6);
  Writer::attr(root, "k0inv_m", L);
  Writer::attr(root, "omega0_per_s", u.omega0_per_fs() * 1e15);
  Writer::attr(root, "n_c_m3", nc);
  Writer::attr(root, "n0_m3", n_ref * nc);
  Writer::attr(root, "E0_V_per_m", E0);
  Writer::attr(root, "B0_T", B0);
  Writer::attr(root, "charge_unit_C", Q0);
  Writer::attr(root, "modes", (long long)(M + 1));
  Writer::attr(root, "Nr", (long long)N);
  Writer::attr(root, "Nzeta", (long long)Nz);
  Writer::attr(root, "dzeta", cfg.grid.dzeta);  // uniform: the step; table: the smallest step
  Writer::attr(root, "zeta_grid", std::string(zg.uniform ? "uniform" : "table"));
  Writer::attr(root, "dz", cfg.grid.stepAt(s0));  // the next step
  Writer::attr(root, "backend", backendInfo());
  Writer::attr(root, "units", std::string("code units: length 1/k0, fields m c omega0/e, momenta m_e c, "
                                          "charge e n_c/k0^3; each dataset has unitSI"));

  // ---- grid
  W.group("/grid");
  std::vector<double> v(g.rho.begin(), g.rho.begin() + N + 1);
  W.real("/grid/r_faces", {hsize_t(N + 1)}, v.data(), L, um, "radii of the faces rho_i (rho_0 = 0: axis)");
  v.assign(g.c.begin(), g.c.begin() + N);
  W.real("/grid/r_centres", {hsize_t(N)}, v.data(), L, um, "radii of the cell centres c_i");
  v.resize(Nz + 1);
  for (int n = 0; n <= Nz; ++n) v[n] = zv.node(n);
  W.real("/grid/zeta_slices", {hsize_t(Nz + 1)}, v.data(), L, um,
         "zeta of the field slices (zeta = z - c t; slice Nz is the empty ghost ahead of the window)");
  v.resize(Nz);
  for (int k = 0; k < Nz; ++k) v[k] = zv.centre(k);
  W.real("/grid/zeta_halfcells", {hsize_t(Nz)}, v.data(), L, um, "zeta of the half-cells (densities)");

  // ---- fields: all modes on the native staggered grid
  W.group("/fields");
  {
    Handle fg(H5Gopen2(W.file(), "/fields", H5P_DEFAULT), H5Gclose, "open /fields");
    Writer::attr(fg, "convention", std::string("F(r,phi) = F_0(r) + 2 Re sum_{m=1..M} F_m(r) exp(i m phi); "
                                                "first index m, second zeta slice, third radius"));
  }
  const size_t nf = fold.Rr.size(), nc_ = fold.Rp.size();
  const std::vector<hsize_t> df = {hsize_t(M + 1), hsize_t(Nz + 1), hsize_t(N + 1)};
  const std::vector<hsize_t> dc = {hsize_t(M + 1), hsize_t(Nz + 1), hsize_t(N)};
  const std::string eu = "m c omega0 / e", bu = "m omega0 / e";
  std::vector<cd> a(nf);
  for (size_t j = 0; j < nf; ++j) a[j] = 0.5 * (fold.Rr[j] + fold.Lr[j]);
  W.complex("/fields/Er", df, a.data(), E0, eu, "E_r, modes on the faces r_faces");
  for (size_t j = 0; j < nf; ++j) a[j] = 0.5 * (fold.Rr[j] - fold.Lr[j]);
  W.complex("/fields/Bphi", df, a.data(), B0, bu, "B_phi, modes on the faces r_faces");
  W.complex("/fields/Bz", df, fold.Bz.data(), B0, bu, "B_z, modes on the faces r_faces");
  a.resize(nc_);
  for (size_t j = 0; j < nc_; ++j) a[j] = 0.5 * (fold.Rp[j] + fold.Lp[j]);
  W.complex("/fields/Ephi", dc, a.data(), E0, eu, "E_phi, modes on the centres r_centres");
  for (size_t j = 0; j < nc_; ++j) a[j] = 0.5 * (fold.Lp[j] - fold.Rp[j]);
  W.complex("/fields/Br", dc, a.data(), B0, bu, "B_r, modes on the centres r_centres");
  W.complex("/fields/Ez", dc, fold.Ez.data(), E0, eu, "E_z, modes on the centres r_centres");

  // ---- densities (m = 0, half-cells x centres)
  W.group("/densities");
  std::vector<double> d(size_t(Nz) * N);
  for (size_t j = 0; j < d.size(); ++j) d[j] = -(c0s.rho[j] + c0b.rho[j]) / n_ref;
  W.real("/densities/ne_over_n0", {hsize_t(Nz), hsize_t(N)}, d.data(), n_ref * nc, "n0",
         "electron density (plasma, bunches, trapped) / n0, m = 0; [half-cell][centre]");
  for (size_t j = 0; j < d.size(); ++j) d[j] = (c0s.rhoi[j] + c0b.rhoi[j]) / n_ref;
  W.real("/densities/Zni_over_n0", {hsize_t(Nz), hsize_t(N)}, d.data(), n_ref * nc, "n0",
         "charge density of the ion species (ion blocks) / (e n0), m = 0; [half-cell][centre]");

  // ---- beam list
  if (O.particles) {
    W.group("/particles");
    W.group("/particles/beam");
    const hsize_t np = beam.size();
    std::vector<double> x(np);
    auto col = [&](const char* nm, double (*get)(const Particle&), double unitSI, const std::string& unit,
                   const char* desc) {
      for (size_t k = 0; k < np; ++k) x[k] = get(beam[k]);
      W.real(std::string("/particles/beam/") + nm, {np}, x.data(), unitSI, unit, desc);
    };
    col("x", [](const Particle& p) { return p.x; }, L, um, "position x");
    col("y", [](const Particle& p) { return p.y; }, L, um, "position y");
    col("z", [](const Particle& p) { return p.z; }, L, um, "lab position z (= s)");
    col("zeta", [](const Particle& p) { return p.zeta; }, L, um, "zeta = z - c t");
    col("ux", [](const Particle& p) { return p.ux; }, 1.0, "m_e c", "momentum gamma v_x / c");
    col("uy", [](const Particle& p) { return p.uy; }, 1.0, "m_e c", "momentum gamma v_y / c");
    col("uz", [](const Particle& p) { return p.uz; }, 1.0, "m_e c", "momentum gamma v_z / c");
    col("charge", [](const Particle& p) { return p.q; }, Q0, "e n_c / k0^3", "macroparticle charge (signed)");
    std::vector<int> sp(np);
    for (size_t k = 0; k < np; ++k) sp[k] = beam[k].sp;
    W.integer("/particles/beam/species", {np}, sp.data(),
              ">= 0: bunch index (see bunch_names), -1: self-trapped plasma electrons, "
              "-30-s: trapped electrons from ionization of ion block s");
    Handle bg(H5Gopen2(W.file(), "/particles/beam", H5P_DEFAULT), H5Gclose, "open /particles/beam");
    std::string names;
    for (size_t b = 0; b < cfg.bunches.size(); ++b) names += (b ? "," : "") + cfg.bunches[b].name;
    Writer::attr(bg, "bunch_names", names);
    std::string ions;
    for (size_t s = 0; s < cfg.ions.size(); ++s) ions += (s ? "," : "") + cfg.ions[s].name;
    Writer::attr(bg, "ion_names", ions);
  }
  W.text("/input_deck", deckText);
}

}  // namespace gem

#else  // no HDF5 in this build: output.format = hdf5 is rejected when the deck is read

namespace gem {
void Simulation::writeHdf5() { throw std::runtime_error("this build has no HDF5 support"); }
void Simulation::writeHdf5File(const std::string&) { throw std::runtime_error("this build has no HDF5 support"); }
void Simulation::writeOpenPMD() { throw std::runtime_error("this build has no HDF5 support"); }
std::string Simulation::writeCheckpoint() { throw std::runtime_error("checkpoints need HDF5 support"); }
void Simulation::readCheckpoint(const std::string&) { throw std::runtime_error("restarts need HDF5 support"); }
std::string Simulation::latestCheckpoint(const std::string&) { return ""; }
}  // namespace gem

#endif
