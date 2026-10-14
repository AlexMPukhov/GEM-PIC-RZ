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
#include <ctime>
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

  static void attrU32(hid_t obj, const char* name, unsigned v) {
    Handle sp(H5Screate(H5S_SCALAR), H5Sclose, "dataspace");
    Handle a(H5Acreate2(obj, name, H5T_STD_U32LE, sp, H5P_DEFAULT, H5P_DEFAULT), H5Aclose, std::string("attribute ") + name);
    check(H5Awrite(a, H5T_NATIVE_UINT, &v), name);
  }
  static void attr(hid_t obj, const char* name, const std::vector<double>& v) {
    const hsize_t n = v.size();
    Handle sp(H5Screate_simple(1, &n, nullptr), H5Sclose, "dataspace");
    Handle a(H5Acreate2(obj, name, H5T_IEEE_F64LE, sp, H5P_DEFAULT, H5P_DEFAULT), H5Aclose, std::string("attribute ") + name);
    check(H5Awrite(a, H5T_NATIVE_DOUBLE, v.data()), name);
  }
  static void attr(hid_t obj, const char* name, const std::vector<unsigned long long>& v) {
    const hsize_t n = v.size();
    Handle sp(H5Screate_simple(1, &n, nullptr), H5Sclose, "dataspace");
    Handle a(H5Acreate2(obj, name, H5T_STD_U64LE, sp, H5P_DEFAULT, H5P_DEFAULT), H5Aclose, std::string("attribute ") + name);
    check(H5Awrite(a, H5T_NATIVE_ULLONG, v.data()), name);
  }
  // array of fixed-length strings (openPMD axisLabels)
  static void attr(hid_t obj, const char* name, const std::vector<std::string>& v) {
    size_t len = 1;
    for (const auto& x : v) len = std::max(len, x.size() + 1);
    std::string buf(len * v.size(), '\0');
    for (size_t k = 0; k < v.size(); ++k) std::copy(v[k].begin(), v[k].end(), buf.begin() + long(k * len));
    Handle t(H5Tcopy(H5T_C_S1), H5Tclose, "string type");
    check(H5Tset_size(t, len), "string size");
    check(H5Tset_strpad(t, H5T_STR_NULLTERM), "string pad");
    const hsize_t n = v.size();
    Handle sp(H5Screate_simple(1, &n, nullptr), H5Sclose, "dataspace");
    Handle a(H5Acreate2(obj, name, t, sp, H5P_DEFAULT, H5P_DEFAULT), H5Aclose, std::string("attribute ") + name);
    check(H5Awrite(a, t, buf.data()), name);
  }
  // group with all missing parents, returned open
  hid_t groups(const std::string& path) {
    size_t p = 1;
    for (;;) {
      p = path.find('/', p);
      const std::string sub = path.substr(0, p);
      if (H5Lexists(file_, sub.c_str(), H5P_DEFAULT) <= 0) group(sub);
      if (p == std::string::npos) break;
      ++p;
    }
    return H5Gopen2(file_, path.c_str(), H5P_DEFAULT);
  }
  // dataset without the descriptive attributes (openPMD sets its own), returned open
  hid_t plain(const std::string& path, const std::vector<hsize_t>& dims, const double* data) {
    hid_t ds = create(path, dims, single_ ? H5T_IEEE_F32LE : H5T_IEEE_F64LE);
    if (ds < 0) throw std::runtime_error("HDF5 error: dataset " + path);
    check(H5Dwrite(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data), "write " + path);
    return ds;
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
  Writer W(name, O.hdf5_single, O.hdf5_compression);

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
  Writer::attr(root, "comment", std::string("quasi-static GT2 code: fields of the window at z = zeta + s, "
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
    Writer::attr(obj, "geometryParameters", "m=" + std::to_string(modes - 1) + ";imag=+");
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

  // ---- particles
  if (!O.particles) return;
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
void Simulation::writeOpenPMD() { throw std::runtime_error("this build has no HDF5 support"); }
std::string Simulation::writeCheckpoint() { throw std::runtime_error("checkpoints need HDF5 support"); }
void Simulation::readCheckpoint(const std::string&) { throw std::runtime_error("restarts need HDF5 support"); }
std::string Simulation::latestCheckpoint(const std::string&) { return ""; }
}  // namespace gem

#endif
