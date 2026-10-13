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
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdio>
#include <sstream>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <vector>

#include "simulation.hpp"

#if defined(GEMPIC_HAVE_HDF5)
#include <hdf5.h>

namespace gem {

namespace {

// ---------------------------------------------------------------- HDF5 helpers
void check(herr_t e, const std::string& what) {
  if (e < 0) throw std::runtime_error("HDF5 error: " + what);
}

// owning handle of an HDF5 object
class Handle {
 public:
  Handle(hid_t id, herr_t (*closer)(hid_t), const std::string& what) : id_(id), close_(closer) {
    if (id < 0) throw std::runtime_error("HDF5 error: " + what);
  }
  Handle(const Handle&) = delete;
  Handle& operator=(const Handle&) = delete;
  ~Handle() { if (id_ >= 0) close_(id_); }
  operator hid_t() const { return id_; }

 private:
  hid_t id_;
  herr_t (*close_)(hid_t);
};

class Writer {
 public:
  Writer(const std::string& path, bool single, int level)
      : file_(H5Fcreate(path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT), H5Fclose, "cannot create " + path),
        single_(single), level_(level) {}

  hid_t file() const { return file_; }

  void group(const std::string& path) {
    Handle g(H5Gcreate2(file_, path.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT), H5Gclose, "group " + path);
  }

  // ---- attributes (on the file root or on an open object)
  static void attr(hid_t obj, const char* name, double v) {
    Handle sp(H5Screate(H5S_SCALAR), H5Sclose, "dataspace");
    Handle a(H5Acreate2(obj, name, H5T_IEEE_F64LE, sp, H5P_DEFAULT, H5P_DEFAULT), H5Aclose, std::string("attribute ") + name);
    check(H5Awrite(a, H5T_NATIVE_DOUBLE, &v), name);
  }
  static void attr(hid_t obj, const char* name, long long v) {
    Handle sp(H5Screate(H5S_SCALAR), H5Sclose, "dataspace");
    Handle a(H5Acreate2(obj, name, H5T_STD_I64LE, sp, H5P_DEFAULT, H5P_DEFAULT), H5Aclose, std::string("attribute ") + name);
    check(H5Awrite(a, H5T_NATIVE_LLONG, &v), name);
  }
  static void attr(hid_t obj, const char* name, const std::string& v) {
    Handle t(H5Tcopy(H5T_C_S1), H5Tclose, "string type");
    check(H5Tset_size(t, v.empty() ? 1 : v.size()), "string size");
    check(H5Tset_strpad(t, H5T_STR_NULLTERM), "string pad");
    Handle sp(H5Screate(H5S_SCALAR), H5Sclose, "dataspace");
    Handle a(H5Acreate2(obj, name, t, sp, H5P_DEFAULT, H5P_DEFAULT), H5Aclose, std::string("attribute ") + name);
    const std::string s = v.empty() ? std::string(1, '\0') : v;
    check(H5Awrite(a, t, s.c_str()), name);
  }

  // ---- datasets
  void real(const std::string& path, const std::vector<hsize_t>& dims, const double* data, double unitSI,
            const std::string& unit, const std::string& desc) {
    Handle ds(create(path, dims, single_ ? H5T_IEEE_F32LE : H5T_IEEE_F64LE), H5Dclose, "dataset " + path);
    check(H5Dwrite(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data), "write " + path);
    describe(ds, unitSI, unit, desc);
  }
  void integer(const std::string& path, const std::vector<hsize_t>& dims, const int* data, const std::string& desc) {
    Handle ds(create(path, dims, H5T_STD_I32LE), H5Dclose, "dataset " + path);
    check(H5Dwrite(ds, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data), "write " + path);
    describe(ds, 1.0, "", desc);
  }
  void u64(const std::string& path, const std::vector<hsize_t>& dims, const unsigned long long* data,
           const std::string& desc) {
    Handle ds(create(path, dims, H5T_STD_U64LE), H5Dclose, "dataset " + path);
    check(H5Dwrite(ds, H5T_NATIVE_ULLONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, data), "write " + path);
    describe(ds, 1.0, "", desc);
  }
  void complex(const std::string& path, const std::vector<hsize_t>& dims, const std::complex<double>* data,
               double unitSI, const std::string& unit, const std::string& desc) {
    Handle ft(complexType(single_ ? H5T_IEEE_F32LE : H5T_IEEE_F64LE), H5Tclose, "complex type");
    Handle mt(complexType(H5T_NATIVE_DOUBLE), H5Tclose, "complex type");
    Handle ds(create(path, dims, ft), H5Dclose, "dataset " + path);
    check(H5Dwrite(ds, mt, H5S_ALL, H5S_ALL, H5P_DEFAULT, data), "write " + path);
    describe(ds, unitSI, unit, desc);
  }
  void text(const std::string& path, const std::string& s) {
    Handle t(H5Tcopy(H5T_C_S1), H5Tclose, "string type");
    check(H5Tset_size(t, s.empty() ? 1 : s.size()), "string size");
    Handle sp(H5Screate(H5S_SCALAR), H5Sclose, "dataspace");
    Handle ds(H5Dcreate2(file_, path.c_str(), t, sp, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT), H5Dclose, "dataset " + path);
    const std::string v = s.empty() ? std::string(1, '\0') : s;
    check(H5Dwrite(ds, t, H5S_ALL, H5S_ALL, H5P_DEFAULT, v.c_str()), "write " + path);
  }

 private:
  Handle file_;
  bool single_;
  int level_;

  // {r, i} compound with members of type `base` (std::complex<double> layout for doubles)
  static hid_t complexType(hid_t base) {
    const size_t sz = H5Tget_size(base);
    hid_t t = H5Tcreate(H5T_COMPOUND, 2 * sz);
    H5Tinsert(t, "r", 0, base);
    H5Tinsert(t, "i", sz, base);
    return t;
  }
  hid_t create(const std::string& path, const std::vector<hsize_t>& dims, hid_t type) {
    Handle sp(H5Screate_simple(int(dims.size()), dims.data(), nullptr), H5Sclose, "dataspace " + path);
    Handle pl(H5Pcreate(H5P_DATASET_CREATE), H5Pclose, "property list");
    bool empty = false;
    for (hsize_t d : dims) empty = empty || d == 0;
    if (level_ > 0 && !empty) {
      // chunks: slabs of at most 256 slices (fields) or 65536 particles
      std::vector<hsize_t> ch(dims);
      if (ch.size() == 1) ch[0] = std::min<hsize_t>(ch[0], 65536);
      else {
        for (size_t k = 0; k + 2 < ch.size(); ++k) ch[k] = 1;
        ch[ch.size() - 2] = std::min<hsize_t>(ch[ch.size() - 2], 256);
      }
      check(H5Pset_chunk(pl, int(ch.size()), ch.data()), "chunk " + path);
      check(H5Pset_shuffle(pl), "shuffle " + path);
      check(H5Pset_deflate(pl, unsigned(level_)), "deflate " + path);
    }
    return H5Dcreate2(file_, path.c_str(), type, sp, H5P_DEFAULT, pl, H5P_DEFAULT);
  }
  static void describe(hid_t ds, double unitSI, const std::string& unit, const std::string& desc) {
    attr(ds, "unitSI", unitSI);
    attr(ds, "unit", unit);
    attr(ds, "description", desc);
  }
};

// ---------------------------------------------------------------- reading (checkpoints)
class Reader {
 public:
  explicit Reader(const std::string& path)
      : file_(H5Fopen(path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT), H5Fclose, "cannot open " + path), path_(path) {}
  bool has(const std::string& p) const { return H5Lexists(file_, p.c_str(), H5P_DEFAULT) > 0; }
  double attrD(const char* name) const { double v = 0; readAttr(name, H5T_NATIVE_DOUBLE, &v); return v; }
  long long attrI(const char* name) const { long long v = 0; readAttr(name, H5T_NATIVE_LLONG, &v); return v; }
  std::string attrS(const char* name) const {
    Handle a(H5Aopen_by_name(file_, "/", name, H5P_DEFAULT, H5P_DEFAULT), H5Aclose, path_ + ": attribute " + name);
    Handle t(H5Aget_type(a), H5Tclose, "type");
    std::string s(H5Tget_size(t), '\0');
    check(H5Aread(a, t, &s[0]), name);
    return s.c_str();
  }
  std::string text(const std::string& p) const {
    Handle ds(H5Dopen2(file_, p.c_str(), H5P_DEFAULT), H5Dclose, path_ + ": no dataset " + p);
    Handle t(H5Dget_type(ds), H5Tclose, "type");
    std::string s(H5Tget_size(t), '\0');
    check(H5Dread(ds, t, H5S_ALL, H5S_ALL, H5P_DEFAULT, &s[0]), p);
    return s.c_str();
  }
  // read a dataset of n elements (n = 0: any size, returned in n) into dst via memory type mt
  template <class T>
  std::vector<T> read(const std::string& p, hid_t mt, size_t expect = size_t(-1)) const {
    Handle ds(H5Dopen2(file_, p.c_str(), H5P_DEFAULT), H5Dclose, path_ + ": no dataset " + p);
    Handle sp(H5Dget_space(ds), H5Sclose, "dataspace");
    const hssize_t n = H5Sget_simple_extent_npoints(sp);
    if (expect != size_t(-1) && size_t(n) != expect)
      throw std::runtime_error(path_ + ": " + p + " has " + std::to_string(n) + " elements, expected " +
                               std::to_string(expect) + " (grid of the deck differs from the checkpoint?)");
    std::vector<T> v(static_cast<size_t>(n));
    if (n > 0) check(H5Dread(ds, mt, H5S_ALL, H5S_ALL, H5P_DEFAULT, v.data()), "read " + p);
    return v;
  }
  std::vector<double> reals(const std::string& p, size_t expect = size_t(-1)) const {
    return read<double>(p, H5T_NATIVE_DOUBLE, expect);
  }
  void complexes(const std::string& p, std::vector<std::complex<double>>& dst) const {
    Handle mt(H5Tcreate(H5T_COMPOUND, sizeof(std::complex<double>)), H5Tclose, "complex type");
    H5Tinsert(mt, "r", 0, H5T_NATIVE_DOUBLE);
    H5Tinsert(mt, "i", sizeof(double), H5T_NATIVE_DOUBLE);
    dst = read<std::complex<double>>(p, mt, dst.size());
  }

 private:
  Handle file_;
  std::string path_;
  void readAttr(const char* name, hid_t mt, void* v) const {
    Handle a(H5Aopen_by_name(file_, "/", name, H5P_DEFAULT, H5P_DEFAULT), H5Aclose, path_ + ": attribute " + name);
    check(H5Aread(a, mt, v), name);
  }
};

// complete particle records (every member, double precision)
void writeParticles(Writer& W, const std::string& g, const std::vector<Particle>& ps) {
  W.group(g);
  const hsize_t n = ps.size();
  std::vector<double> d(n);
  auto dbl = [&](const char* nm, double Particle::*m) {
    for (size_t k = 0; k < n; ++k) d[k] = ps[k].*m;
    W.real(g + "/" + nm, {n}, d.data(), 1.0, "code units", nm);
  };
  dbl("x", &Particle::x); dbl("y", &Particle::y); dbl("z", &Particle::z);
  dbl("ux", &Particle::ux); dbl("uy", &Particle::uy); dbl("uz", &Particle::uz);
  dbl("q", &Particle::q); dbl("qm", &Particle::qm); dbl("zeta", &Particle::zeta); dbl("w", &Particle::w);
  std::vector<int> iv(n);
  auto itg = [&](const char* nm, int Particle::*m) {
    for (size_t k = 0; k < n; ++k) iv[k] = ps[k].*m;
    W.integer(g + "/" + nm, {n}, iv.data(), nm);
  };
  itg("status", &Particle::status); itg("sp", &Particle::sp); itg("Zi", &Particle::Zi);
  std::vector<unsigned long long> id(n);
  for (size_t k = 0; k < n; ++k) id[k] = ps[k].id;
  W.u64(g + "/id", {n}, id.data(), "key of the ionization random numbers");
}

std::vector<Particle> readParticles(const Reader& R, const std::string& g) {
  const std::vector<double> x = R.reals(g + "/x");
  const size_t n = x.size();
  std::vector<Particle> ps(n, Particle{});
  auto dbl = [&](const char* nm, double Particle::*m) {
    const std::vector<double> v = R.reals(g + "/" + nm, n);
    for (size_t k = 0; k < n; ++k) ps[k].*m = v[k];
  };
  dbl("x", &Particle::x); dbl("y", &Particle::y); dbl("z", &Particle::z);
  dbl("ux", &Particle::ux); dbl("uy", &Particle::uy); dbl("uz", &Particle::uz);
  dbl("q", &Particle::q); dbl("qm", &Particle::qm); dbl("zeta", &Particle::zeta); dbl("w", &Particle::w);
  auto itg = [&](const char* nm, int Particle::*m) {
    const std::vector<int> v = R.read<int>(g + "/" + nm, H5T_NATIVE_INT, n);
    for (size_t k = 0; k < n; ++k) ps[k].*m = v[k];
  };
  itg("status", &Particle::status); itg("sp", &Particle::sp); itg("Zi", &Particle::Zi);
  const std::vector<unsigned long long> id = R.read<unsigned long long>(g + "/id", H5T_NATIVE_ULLONG, n);
  for (size_t k = 0; k < n; ++k) ps[k].id = id[k];
  return ps;
}

void writeCurrents(Writer& W, const std::string& g, const Currents& C) {
  W.group(g);
  W.complex(g + "/jr", {C.jr.size()}, C.jr.data(), 1.0, "code units", "j_r modes [m][half-cell][face]");
  W.complex(g + "/jp", {C.jp.size()}, C.jp.data(), 1.0, "code units", "j_phi modes [m][half-cell][centre]");
  W.complex(g + "/jz", {C.jz.size()}, C.jz.data(), 1.0, "code units", "j_z modes [m][half-cell][centre]");
  W.real(g + "/rho", {C.rho.size()}, C.rho.data(), 1.0, "code units", "charge density (electrons, bunches), m = 0");
  W.real(g + "/rhoi", {C.rhoi.size()}, C.rhoi.data(), 1.0, "code units", "charge density of the ions, m = 0");
}
void readCurrents(const Reader& R, const std::string& g, Currents& C) {
  R.complexes(g + "/jr", C.jr);
  R.complexes(g + "/jp", C.jp);
  R.complexes(g + "/jz", C.jz);
  C.rho = R.reals(g + "/rho", C.rho.size());
  C.rhoi = R.reals(g + "/rhoi", C.rhoi.size());
}

// step number of a checkpoint file name, -1 if it is none
int checkpointStep(const std::string& fname) {
  int step = -1;
  char tail[8] = {0};
  if (std::sscanf(fname.c_str(), "checkpoint_%d.h%2s", &step, tail) == 2 && std::string(tail) == "5") return step;
  return -1;
}

}  // namespace

// ============================================================================
void Simulation::writeHdf5() {
  namespace fs = std::filesystem;
  const auto& O = cfg.output;
  fs::create_directories(O.dir);
  char name[1024];
  std::snprintf(name, sizeof name, "%s/data_%05d.h5", O.dir.c_str(), step_no);
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

// ============================================================================
//  checkpoints: the complete state, in double precision, so that a restarted
//  run continues exactly as the uninterrupted one
//
//  /                         attributes: kind = "checkpoint", step, z, seed, counters,
//                            grid parameters (checked at the restart)
//  /grid/r_faces, r_centres  radial grid (checked at the restart)
//  /grid/zeta_slices         longitudinal grid (checked at the restart)
//  /state/fields/Rr ... Ez   raw characteristics R, L and E_z, B_z of level s^0
//  /state/beam_currents, /state/stream_currents   j_r, j_phi, j_z, rho, rho_i at s^0
//  /state/beam               beam list (every member of the particle record)
//  /state/pending            particles that continue streaming in the next step
//  /state/rng                state of the random generator of the plasma loading
//  /input_deck               deck of the run that wrote the checkpoint
// ============================================================================
std::string Simulation::writeCheckpoint() {
  namespace fs = std::filesystem;
  if (engine) engine->downloadFull();  // currents and pending particles of the device
  const auto& O = cfg.output;
  fs::create_directories(O.dir);
  char name[1024];
  std::snprintf(name, sizeof name, "%s/checkpoint_%05d.h5", O.dir.c_str(), step_no);
  const std::string path = name, tmp = path + ".tmp";
  H5Eset_auto2(H5E_DEFAULT, nullptr, nullptr);
  {
    Writer W(tmp, /*single=*/false, /*compression=*/0);
    const hid_t root = W.file();
    Writer::attr(root, "code", std::string("GEM-PIC-RZ"));
    Writer::attr(root, "kind", std::string("checkpoint"));
    Writer::attr(root, "format_version", 1LL);
    Writer::attr(root, "step", (long long)step_no);
    Writer::attr(root, "z", s0);
    Writer::attr(root, "seed", (long long)cfg.numerics.seed);
    Writer::attr(root, "n_restreamed", (long long)n_restreamed);
    Writer::attr(root, "n_trapped_total", (long long)n_trapped_total);
    Writer::attr(root, "lambda0_um", cfg.units.lambda0_um);
    Writer::attr(root, "modes", (long long)(M + 1));
    Writer::attr(root, "Nr", (long long)g.N);
    Writer::attr(root, "Nzeta", (long long)Nz);
    Writer::attr(root, "dzeta", cfg.grid.dzeta);
    Writer::attr(root, "dz", cfg.grid.stepAt(s0));  // informative: the z grid may change at a restart
    std::string ions;
    for (size_t s = 0; s < cfg.ions.size(); ++s) ions += (s ? "," : "") + cfg.ions[s].name;
    Writer::attr(root, "ion_names", ions);
    Writer::attr(root, "bunches", (long long)cfg.bunches.size());
    Writer::attr(root, "backend", backendInfo());
    W.group("/grid");
    W.real("/grid/r_faces", {hsize_t(g.N + 1)}, g.rho.data(), 1.0, "1/k0", "faces");
    W.real("/grid/r_centres", {hsize_t(g.N)}, g.c.data(), 1.0, "1/k0", "centres");
    W.real("/grid/zeta_slices", {hsize_t(Nz + 1)}, zg.z.data(), 1.0, "1/k0", "slices");
    W.group("/state");
    W.group("/state/fields");
    const char* desc = "[m][slice][r], level s^0";
    W.complex("/state/fields/Rr", {fold.Rr.size()}, fold.Rr.data(), 1.0, "code units", desc);
    W.complex("/state/fields/Lr", {fold.Lr.size()}, fold.Lr.data(), 1.0, "code units", desc);
    W.complex("/state/fields/Bz", {fold.Bz.size()}, fold.Bz.data(), 1.0, "code units", desc);
    W.complex("/state/fields/Rp", {fold.Rp.size()}, fold.Rp.data(), 1.0, "code units", desc);
    W.complex("/state/fields/Lp", {fold.Lp.size()}, fold.Lp.data(), 1.0, "code units", desc);
    W.complex("/state/fields/Ez", {fold.Ez.size()}, fold.Ez.data(), 1.0, "code units", desc);
    writeCurrents(W, "/state/beam_currents", c0b);
    writeCurrents(W, "/state/stream_currents", c0s);
    writeParticles(W, "/state/beam", beam);
    writeParticles(W, "/state/pending", pendingNext);
    std::ostringstream os;
    os << rng;
    W.text("/state/rng", os.str());
    W.text("/input_deck", deckText);
  }
  fs::rename(tmp, path);  // a crash while writing never leaves a broken checkpoint_NNNNN.h5
  // keep only the newest checkpoint_keep files
  if (O.checkpoint_keep > 0) {
    std::vector<std::pair<int, fs::path>> files;
    for (const auto& e : fs::directory_iterator(O.dir)) {
      const int st = checkpointStep(e.path().filename().string());
      if (st >= 0) files.push_back({st, e.path()});
    }
    std::sort(files.begin(), files.end());
    for (size_t k = 0; k + O.checkpoint_keep < files.size(); ++k) fs::remove(files[k].second);
  }
  return path;
}

void Simulation::readCheckpoint(const std::string& path) {
  H5Eset_auto2(H5E_DEFAULT, nullptr, nullptr);
  Reader R(path);
  auto fail = [&](const std::string& what) {
    throw std::runtime_error("restart from " + path + ": " + what);
  };
  if (!R.has("/state") || R.attrS("kind") != "checkpoint") fail("not a checkpoint file");
  // the grid must be the one of the checkpoint
  auto differs = [](double a, double b) { return std::fabs(a - b) > 1e-12 * std::max(std::fabs(a), std::fabs(b)); };
  if (R.attrI("modes") != M + 1) fail("m_max differs from the checkpoint");
  if (R.attrI("Nr") != g.N) fail("number of radial cells differs from the checkpoint");
  if (R.attrI("Nzeta") != Nz) fail("nzeta differs from the checkpoint");
  if (R.has("/grid/zeta_slices")) {
    const std::vector<double> zs = R.reals("/grid/zeta_slices", zg.z.size());
    for (size_t n = 0; n < zs.size(); ++n)
      if (differs(zs[n], zg.z[n])) fail("the zeta grid differs from the checkpoint");
  } else if (!zg.uniform || differs(R.attrD("dzeta"), cfg.grid.dzeta)) {
    fail("dzeta differs from the checkpoint");
  }
  // (dz and the z grid may change: the scheme takes any step from s^0)
  if (differs(R.attrD("lambda0_um"), cfg.units.lambda0_um)) fail("lambda0 differs from the checkpoint");
  const std::vector<double> rf = R.reals("/grid/r_faces", g.rho.size());
  for (size_t i = 0; i < rf.size(); ++i)
    if (differs(rf[i], g.rho[i])) fail("the radial grid differs from the checkpoint");
  std::string ions;
  for (size_t s = 0; s < cfg.ions.size(); ++s) ions += (s ? "," : "") + cfg.ions[s].name;
  if (R.attrS("ion_names") != ions) fail("the ion blocks differ from the checkpoint ('" + R.attrS("ion_names") + "')");
  if (R.attrI("bunches") != (long long)cfg.bunches.size())
    std::printf("WARNING: the deck has %zu bunch blocks, the checkpoint %lld (bunch names in the output may be wrong)\n",
                cfg.bunches.size(), R.attrI("bunches"));
  if (R.attrI("seed") != (long long)cfg.numerics.seed)
    std::printf("WARNING: numerics.seed differs from the checkpoint: the ionization random numbers differ from the "
                "uninterrupted run\n");
  // state
  R.complexes("/state/fields/Rr", fold.Rr);
  R.complexes("/state/fields/Lr", fold.Lr);
  R.complexes("/state/fields/Bz", fold.Bz);
  R.complexes("/state/fields/Rp", fold.Rp);
  R.complexes("/state/fields/Lp", fold.Lp);
  R.complexes("/state/fields/Ez", fold.Ez);
  readCurrents(R, "/state/beam_currents", c0b);
  readCurrents(R, "/state/stream_currents", c0s);
  beam = readParticles(R, "/state/beam");
  pendingNext = readParticles(R, "/state/pending");
  std::istringstream is(R.text("/state/rng"));
  is >> rng;
  if (!is) fail("cannot restore the random generator");
  step_no = int(R.attrI("step"));
  s0 = R.attrD("z");
  n_restreamed = long(R.attrI("n_restreamed"));
  n_trapped_total = long(R.attrI("n_trapped_total"));
}

std::string Simulation::latestCheckpoint(const std::string& dir) {
  namespace fs = std::filesystem;
  int best = -1;
  std::string path;
  if (fs::is_directory(dir))
    for (const auto& e : fs::directory_iterator(dir)) {
      const int st = checkpointStep(e.path().filename().string());
      if (st > best) { best = st; path = e.path().string(); }
    }
  return path;
}

}  // namespace gem

#else  // no HDF5 in this build: output.format = hdf5 is rejected when the deck is read

namespace gem {
void Simulation::writeHdf5() { throw std::runtime_error("this build has no HDF5 support"); }
std::string Simulation::writeCheckpoint() { throw std::runtime_error("checkpoints need HDF5 support"); }
void Simulation::readCheckpoint(const std::string&) { throw std::runtime_error("restarts need HDF5 support"); }
std::string Simulation::latestCheckpoint(const std::string&) { return ""; }
}  // namespace gem

#endif
