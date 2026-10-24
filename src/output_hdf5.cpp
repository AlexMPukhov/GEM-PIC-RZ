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
#include "version.hpp"
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
  // largeAttributes: file format of HDF5 >= 1.8, which stores attributes larger than 64 kB
  // (the node coordinates of the native openPMD records)
  Writer(const std::string& path, bool single, int level, bool largeAttributes = false)
      : file_(create(path, largeAttributes), H5Fclose, "cannot create " + path), single_(single), level_(level) {}

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
  // dataset without the descriptive attributes (openPMD sets its own), returned open;
  // forceDouble: always float64 (grid coordinates)
  hid_t plain(const std::string& path, const std::vector<hsize_t>& dims, const double* data, bool forceDouble = false) {
    hid_t ds = create(path, dims, (single_ && !forceDouble) ? H5T_IEEE_F32LE : H5T_IEEE_F64LE);
    if (ds < 0) throw std::runtime_error("HDF5 error: dataset " + path);
    check(H5Dwrite(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data), "write " + path);
    return ds;
  }

  // ---- datasets written in pieces (block by block): created empty, filled by hyperslabs
  hid_t empty(const std::string& path, const std::vector<hsize_t>& dims, bool complexData, bool forceDouble = false) {
    const hid_t base = (single_ && !forceDouble) ? H5T_IEEE_F32LE : H5T_IEEE_F64LE;
    hid_t ds;
    if (complexData) {
      Handle ct(complexType(base), H5Tclose, "complex type");
      ds = create(path, dims, ct);
    } else {
      ds = create(path, dims, base);
    }
    if (ds < 0) throw std::runtime_error("HDF5 error: dataset " + path);
    return ds;
  }
  // the box start + [0, count) of ds from data (C order of count); complexData: std::complex<double>
  static void slab(hid_t ds, const std::vector<hsize_t>& start, const std::vector<hsize_t>& count, const void* data,
                   bool complexData) {
    for (hsize_t c : count)
      if (c == 0) return;
    Handle fs(H5Dget_space(ds), H5Sclose, "dataspace");
    check(H5Sselect_hyperslab(fs, H5S_SELECT_SET, start.data(), nullptr, count.data(), nullptr), "hyperslab");
    Handle ms(H5Screate_simple(int(count.size()), count.data(), nullptr), H5Sclose, "memory space");
    if (complexData) {
      Handle mt(complexType(H5T_NATIVE_DOUBLE), H5Tclose, "complex type");
      check(H5Dwrite(ds, mt, ms, fs, H5P_DEFAULT, data), "write slab");
    } else {
      check(H5Dwrite(ds, H5T_NATIVE_DOUBLE, ms, fs, H5P_DEFAULT, data), "write slab");
    }
  }
  static void describeDs(hid_t ds, double unitSI, const std::string& unit, const std::string& desc) {
    describe(ds, unitSI, unit, desc);
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

  static hid_t create(const std::string& path, bool largeAttributes) {
    if (!largeAttributes) return H5Fcreate(path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    const hid_t fapl = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_libver_bounds(fapl, H5F_LIBVER_V18, H5F_LIBVER_LATEST);
    const hid_t f = H5Fcreate(path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, fapl);
    H5Pclose(fapl);
    return f;
  }
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
  // the blocks {offset, count} of a 1-D dataset of `total` elements, in that order, into dst (its size: the sum)
  template <class T>
  void readBlocks(const std::string& p, hid_t mt, size_t total, const std::vector<std::pair<hsize_t, hsize_t>>& blocks,
                  std::vector<T>& dst) const {
    Handle ds(H5Dopen2(file_, p.c_str(), H5P_DEFAULT), H5Dclose, path_ + ": no dataset " + p);
    Handle sp(H5Dget_space(ds), H5Sclose, "dataspace");
    if (size_t(H5Sget_simple_extent_npoints(sp)) != total)
      throw std::runtime_error(path_ + ": " + p + " has " + std::to_string(H5Sget_simple_extent_npoints(sp)) +
                               " elements, expected " + std::to_string(total) + " (grid of the deck differs?)");
    hsize_t n = 0;
    check(H5Sselect_none(sp), "select");
    for (const auto& b : blocks) {
      if (b.second == 0) continue;
      check(H5Sselect_hyperslab(sp, H5S_SELECT_OR, &b.first, nullptr, &b.second, nullptr), "hyperslab " + p);
      n += b.second;
    }
    if (n != dst.size()) throw std::runtime_error("readBlocks: size mismatch for " + p);
    if (n == 0) return;
    Handle ms(H5Screate_simple(1, &n, nullptr), H5Sclose, "memory space");
    check(H5Dread(ds, mt, ms, sp, H5P_DEFAULT, dst.data()), "read " + p);
  }
  void complexBlocks(const std::string& p, size_t total, const std::vector<std::pair<hsize_t, hsize_t>>& blocks,
                     std::vector<std::complex<double>>& dst) const {
    Handle mt(H5Tcreate(H5T_COMPOUND, sizeof(std::complex<double>)), H5Tclose, "complex type");
    H5Tinsert(mt, "r", 0, H5T_NATIVE_DOUBLE);
    H5Tinsert(mt, "i", sizeof(double), H5T_NATIVE_DOUBLE);
    readBlocks(p, mt, total, blocks, dst);
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
// the stored half-cells [k0, k0 + nk) of C (the whole window, or an MPI block)
void readCurrents(const Reader& R, const std::string& g, Currents& C) {
  auto blocks = [&](size_t L) {
    std::vector<std::pair<hsize_t, hsize_t>> b;
    for (int m = 0; m <= C.M; ++m) b.push_back({hsize_t((size_t(m) * C.Nz + C.k0) * L), hsize_t(size_t(C.nk) * L)});
    return b;
  };
  const size_t nf = size_t(C.M + 1) * C.Nz * (C.N + 1), nc = size_t(C.M + 1) * C.Nz * C.N;
  R.complexBlocks(g + "/jr", nf, blocks(C.N + 1), C.jr);
  R.complexBlocks(g + "/jp", nc, blocks(C.N), C.jp);
  R.complexBlocks(g + "/jz", nc, blocks(C.N), C.jz);
  const std::vector<std::pair<hsize_t, hsize_t>> d = {{hsize_t(size_t(C.k0) * C.N), hsize_t(size_t(C.nk) * C.N)}};
  R.readBlocks(g + "/rho", H5T_NATIVE_DOUBLE, size_t(C.Nz) * C.N, d, C.rho);
  R.readBlocks(g + "/rhoi", H5T_NATIVE_DOUBLE, size_t(C.Nz) * C.N, d, C.rhoi);
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
  // fields: the stored slices only (the whole window, or the block of an MPI rank)
  auto fblocks = [&](size_t L) {
    std::vector<std::pair<hsize_t, hsize_t>> b;
    for (int m = 0; m <= M; ++m) b.push_back({hsize_t((size_t(m) * (Nz + 1) + fold.n0) * L), hsize_t(size_t(fold.ns) * L)});
    return b;
  };
  const size_t nff = size_t(M + 1) * (Nz + 1) * (g.N + 1), ncf = size_t(M + 1) * (Nz + 1) * g.N;
  R.complexBlocks("/state/fields/Rr", nff, fblocks(g.N + 1), fold.Rr);
  R.complexBlocks("/state/fields/Lr", nff, fblocks(g.N + 1), fold.Lr);
  R.complexBlocks("/state/fields/Bz", nff, fblocks(g.N + 1), fold.Bz);
  R.complexBlocks("/state/fields/Rp", ncf, fblocks(g.N), fold.Rp);
  R.complexBlocks("/state/fields/Lp", ncf, fblocks(g.N), fold.Lp);
  R.complexBlocks("/state/fields/Ez", ncf, fblocks(g.N), fold.Ez);
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

// ============================================================================
//  LASY envelope file (output.lasy): <dir>/lasy/laser_NNNNN.h5
//
//  The level s is the plane z = s, its slices the times t = s - zeta + lasyTref. The forward field R/2 is
//  interpolated (as in the particle gather) to uniform nodes r_j = j dr, t_l = t_min + l dt, turned into E_x,
//  E_y at 2 M + 4 angles, demodulated at omega_L (analytic signal through an FFT, spectrum kept for
//  |omega - omega_L| < lasy_bandwidth omega_L), and projected on the dominant polarization p (eigenvector of
//  the coherence matrix): env = conj(p_x) env_x + conj(p_y) env_y. The azimuthal modes of env, |k| <= M + 1,
//  are written in the openPMD thetaMode form C_0, C_m, S_m (env = C_0 + sum C_m cos + S_m sin), axes (t, r),
//  as LASY writes them: lasy.profiles.FromOpenPMDProfile(path, "laserEnvelope") reads the file.
// ============================================================================

// ============================================================================
//  Block-wise output (Simulation::beginOutput / outputBlock / endOutput): the
//  datasets of the window are created at their full size and filled block by
//  block with hyperslabs, so the writer never needs the whole window in memory.
//  Every block [n0, n1) is passed with fold holding the slices [n0, n1] and the
//  currents the half-cells [n0, n1] (one beyond the block, for the
//  interpolations); the last block (n1 = Nz) also writes the ghost slice Nz.
// ============================================================================

namespace {
// physical units of the files
struct SIunits {
  double L, E0, B0, T0, nc, Q0;
  explicit SIunits(const Units& u) {
    L = u.k0inv_um() * 1e-6;
    E0 = MC2_MEV * 1e6 / L;
    B0 = E0 / C_LIGHT;
    T0 = 1.0 / (u.omega0_per_fs() * 1e15);
    nc = u.nc_cm3() * 1e6;
    Q0 = u.codeToPC(1.0) * 1e-12;
  }
};
constexpr double ME_KG = 9.1093837015e-31;
std::vector<double> unitDim(double Lx, double Mx, double Tx, double Ix) { return {Lx, Mx, Tx, Ix, 0, 0, 0}; }
}  // namespace

// ---------------------------------------------------------------- data_NNNNN.h5
struct Simulation::Hdf5Sink : Simulation::Sink {
  Simulation& S;
  std::string name, tmp;
  std::unique_ptr<Writer> W;
  hid_t dsF[6] = {-1, -1, -1, -1, -1, -1}, dsD[2] = {-1, -1};
  explicit Hdf5Sink(Simulation& s) : S(s) {
    namespace fs = std::filesystem;
    const auto& O = S.cfg.output;
    fs::create_directories(O.dir);
    char nm[1024];
    std::snprintf(nm, sizeof nm, "%s/data_%05d.h5", O.dir.c_str(), S.step_no);
    // written as name.tmp and renamed when complete: a reader (e.g. tools/gempic_gui.py during a run)
    // never sees a partial data_NNNNN.h5
    name = nm;
    tmp = name + ".tmp";
    H5Eset_auto2(H5E_DEFAULT, nullptr, nullptr);  // errors are reported by exceptions
    W = std::make_unique<Writer>(tmp, O.hdf5_single, O.hdf5_compression);
    const SIunits si(S.cfg.units);
    const Units& u = S.cfg.units;
    const int N = S.g.N, M = S.M, Nz = S.Nz;
    const std::string um = "1/k0 = " + std::to_string(u.k0inv_um()) + " um";
    const hid_t root = W->file();
    Writer::attr(root, "code", std::string("GEM-PIC-RZ"));
    Writer::attr(root, "format_version", 1LL);
    Writer::attr(root, "step", (long long)S.step_no);
    Writer::attr(root, "z", S.s0);
    Writer::attr(root, "z_m", S.s0 * si.L);
    Writer::attr(root, "lambda0_m", u.lambda0_um * 1e-6);
    Writer::attr(root, "k0inv_m", si.L);
    Writer::attr(root, "omega0_per_s", u.omega0_per_fs() * 1e15);
    Writer::attr(root, "n_c_m3", si.nc);
    Writer::attr(root, "n0_m3", S.n_ref * si.nc);
    Writer::attr(root, "E0_V_per_m", si.E0);
    Writer::attr(root, "B0_T", si.B0);
    Writer::attr(root, "charge_unit_C", si.Q0);
    Writer::attr(root, "modes", (long long)(M + 1));
    Writer::attr(root, "Nr", (long long)N);
    Writer::attr(root, "Nzeta", (long long)Nz);
    Writer::attr(root, "dzeta", S.cfg.grid.dzeta);  // uniform: the step; table: the smallest step
    Writer::attr(root, "zeta_grid", std::string(S.zg.uniform ? "uniform" : "table"));
    Writer::attr(root, "dz", S.cfg.grid.stepAt(S.s0));  // the next step
    Writer::attr(root, "backend", S.backendInfo());
    Writer::attr(root, "units", std::string("code units: length 1/k0, fields m c omega0/e, momenta m_e c, "
                                            "charge e n_c/k0^3; each dataset has unitSI"));
    W->group("/grid");
    std::vector<double> v(S.g.rho.begin(), S.g.rho.begin() + N + 1);
    W->real("/grid/r_faces", {hsize_t(N + 1)}, v.data(), si.L, um, "radii of the faces rho_i (rho_0 = 0: axis)");
    v.assign(S.g.c.begin(), S.g.c.begin() + N);
    W->real("/grid/r_centres", {hsize_t(N)}, v.data(), si.L, um, "radii of the cell centres c_i");
    v.resize(Nz + 1);
    for (int n = 0; n <= Nz; ++n) v[n] = S.zv.node(n);
    W->real("/grid/zeta_slices", {hsize_t(Nz + 1)}, v.data(), si.L, um,
            "zeta of the field slices (zeta = z - c t; slice Nz is the empty ghost ahead of the window)");
    v.resize(Nz);
    for (int k = 0; k < Nz; ++k) v[k] = S.zv.centre(k);
    W->real("/grid/zeta_halfcells", {hsize_t(Nz)}, v.data(), si.L, um, "zeta of the half-cells (densities)");
    W->group("/fields");
    {
      Handle fg(H5Gopen2(W->file(), "/fields", H5P_DEFAULT), H5Gclose, "open /fields");
      Writer::attr(fg, "convention", std::string("F(r,phi) = F_0(r) + 2 Re sum_{m=1..M} F_m(r) exp(i m phi); "
                                                  "first index m, second zeta slice, third radius"));
    }
    const std::vector<hsize_t> df = {hsize_t(M + 1), hsize_t(Nz + 1), hsize_t(N + 1)};
    const std::vector<hsize_t> dc = {hsize_t(M + 1), hsize_t(Nz + 1), hsize_t(N)};
    const std::string eu = "m c omega0 / e", bu = "m omega0 / e";
    const char* nmF[6] = {"/fields/Er", "/fields/Bphi", "/fields/Bz", "/fields/Ephi", "/fields/Br", "/fields/Ez"};
    const char* dsc[6] = {"E_r, modes on the faces r_faces", "B_phi, modes on the faces r_faces",
                          "B_z, modes on the faces r_faces", "E_phi, modes on the centres r_centres",
                          "B_r, modes on the centres r_centres", "E_z, modes on the centres r_centres"};
    const bool isE[6] = {true, false, false, true, false, true};
    for (int q = 0; q < 6; ++q) {
      dsF[q] = W->empty(nmF[q], q < 3 ? df : dc, true);
      Writer::describeDs(dsF[q], isE[q] ? si.E0 : si.B0, isE[q] ? eu : bu, dsc[q]);
    }
    W->group("/densities");
    dsD[0] = W->empty("/densities/ne_over_n0", {hsize_t(Nz), hsize_t(N)}, false);
    Writer::describeDs(dsD[0], S.n_ref * si.nc, "n0",
                       "electron density (plasma, bunches, trapped) / n0, m = 0; [half-cell][centre]");
    dsD[1] = W->empty("/densities/Zni_over_n0", {hsize_t(Nz), hsize_t(N)}, false);
    Writer::describeDs(dsD[1], S.n_ref * si.nc, "n0",
                       "charge density of the ion species (ion blocks) / (e n0), m = 0; [half-cell][centre]");
  }
  ~Hdf5Sink() override { close(); }
  void close() {
    for (auto& d : dsF) if (d >= 0) { H5Dclose(d); d = -1; }
    for (auto& d : dsD) if (d >= 0) { H5Dclose(d); d = -1; }
  }
  void block(int n0, int n1) override {
    const int N = S.g.N, M = S.M, Nz = S.Nz;
    const int ne = (n1 == Nz) ? Nz + 1 : n1, nb = ne - n0;
    const Fields& F = S.fold;
    for (int q = 0; q < 6; ++q) {
      const bool face = q < 3;
      const int L = face ? N + 1 : N;
      std::vector<cd> a(size_t(M + 1) * nb * L);
      for (int m = 0; m <= M; ++m)
        for (int n = n0; n < ne; ++n)
          for (int i = 0; i < L; ++i) {
            const size_t jf = face ? F.f(m, n, i) : F.c(m, n, i);
            cd v;
            switch (q) {
              case 0: v = 0.5 * (F.Rr[jf] + F.Lr[jf]); break;
              case 1: v = 0.5 * (F.Rr[jf] - F.Lr[jf]); break;
              case 2: v = F.Bz[jf]; break;
              case 3: v = 0.5 * (F.Rp[jf] + F.Lp[jf]); break;
              case 4: v = 0.5 * (F.Lp[jf] - F.Rp[jf]); break;
              default: v = F.Ez[jf];
            }
            a[(size_t(m) * nb + n - n0) * L + i] = v;
          }
      Writer::slab(dsF[q], {0, hsize_t(n0), 0}, {hsize_t(M + 1), hsize_t(nb), hsize_t(L)}, a.data(), true);
    }
    const int kb = n1 - n0;
    std::vector<double> d(size_t(kb) * N);
    for (int k = n0; k < n1; ++k)
      for (int i = 0; i < N; ++i)
        d[size_t(k - n0) * N + i] = -(S.c0s.rho[S.c0s.d(k, i)] + S.c0b.rho[S.c0b.d(k, i)]) / S.n_ref;
    Writer::slab(dsD[0], {hsize_t(n0), 0}, {hsize_t(kb), hsize_t(N)}, d.data(), false);
    for (int k = n0; k < n1; ++k)
      for (int i = 0; i < N; ++i)
        d[size_t(k - n0) * N + i] = (S.c0s.rhoi[S.c0s.d(k, i)] + S.c0b.rhoi[S.c0b.d(k, i)]) / S.n_ref;
    Writer::slab(dsD[1], {hsize_t(n0), 0}, {hsize_t(kb), hsize_t(N)}, d.data(), false);
  }
  std::string finish() override {
    close();
    const auto& O = S.cfg.output;
    const SIunits si(S.cfg.units);
    const std::string um = "1/k0 = " + std::to_string(S.cfg.units.k0inv_um()) + " um";
    if (O.particles) {
      W->group("/particles");
      W->group("/particles/beam");
      const auto& beam = S.beam;
      const hsize_t np = beam.size();
      std::vector<double> x(np);
      auto col = [&](const char* nm, double (*get)(const Particle&), double unitSI, const std::string& unit,
                     const char* desc) {
        for (size_t k = 0; k < np; ++k) x[k] = get(beam[k]);
        W->real(std::string("/particles/beam/") + nm, {np}, x.data(), unitSI, unit, desc);
      };
      col("x", [](const Particle& p) { return p.x; }, si.L, um, "position x");
      col("y", [](const Particle& p) { return p.y; }, si.L, um, "position y");
      col("z", [](const Particle& p) { return p.z; }, si.L, um, "lab position z (= s)");
      col("zeta", [](const Particle& p) { return p.zeta; }, si.L, um, "zeta = z - c t");
      col("ux", [](const Particle& p) { return p.ux; }, 1.0, "m_e c", "momentum gamma v_x / c");
      col("uy", [](const Particle& p) { return p.uy; }, 1.0, "m_e c", "momentum gamma v_y / c");
      col("uz", [](const Particle& p) { return p.uz; }, 1.0, "m_e c", "momentum gamma v_z / c");
      col("charge", [](const Particle& p) { return p.q; }, si.Q0, "e n_c / k0^3", "macroparticle charge (signed)");
      std::vector<int> sp(np);
      for (size_t k = 0; k < np; ++k) sp[k] = beam[k].sp;
      W->integer("/particles/beam/species", {np}, sp.data(),
                 ">= 0: bunch index (see bunch_names), -1: self-trapped plasma electrons, "
                 "-30-s: trapped electrons from ionization of ion block s");
      Handle bg(H5Gopen2(W->file(), "/particles/beam", H5P_DEFAULT), H5Gclose, "open /particles/beam");
      std::string names;
      for (size_t b = 0; b < S.cfg.bunches.size(); ++b) names += (b ? "," : "") + S.cfg.bunches[b].name;
      Writer::attr(bg, "bunch_names", names);
      std::string ions;
      for (size_t s = 0; s < S.cfg.ions.size(); ++s) ions += (s ? "," : "") + S.cfg.ions[s].name;
      Writer::attr(bg, "ion_names", ions);
    }
    W->text("/input_deck", S.deckText);
    W.reset();
    std::filesystem::rename(tmp, name);
    return "";
  }
};

// ---------------------------------------------------------------- openPMD
struct Simulation::OpenPMDSink : Simulation::Sink {
  Simulation& S;
  std::unique_ptr<Writer> W;
  SIunits si;
  std::string base, mp, name;
  int N, M, Nz, NM, nr = 0, nz = 0;
  double dr = 0, dzo = 0;
  struct Rst { int modeC, ac, modeF, af; double wc, wf, xr; };
  std::vector<Rst> rs;
  std::vector<int> zn, zb;        // uniform grid: slice stencil per zeta node; owner half-cell of density rows
  std::vector<double> zw, zc;
  std::vector<int> da;            // density: radial stencil per output radius
  std::vector<double> dw;
  std::vector<hid_t> ds;          // open datasets
  hid_t dsU[6] = {-1, -1, -1, -1, -1, -1}, dsUn[2] = {-1, -1}, dsN[6] = {-1, -1, -1, -1, -1, -1}, dsNn[2] = {-1, -1};
  const cd I = cd(0.0, 1.0);

  hid_t keep(hid_t h) { ds.push_back(h); return h; }
  explicit OpenPMDSink(Simulation& s) : S(s), si(s.cfg.units) {
    namespace fs = std::filesystem;
    const auto& O = S.cfg.output;
    N = S.g.N; M = S.M; Nz = S.Nz; NM = 2 * M + 1;
    const std::string dir = O.dir + "/openpmd";
    fs::create_directories(dir);
    name = dir + "/openpmd_" + std::to_string(S.step_no) + ".h5";
    H5Eset_auto2(H5E_DEFAULT, nullptr, nullptr);
    W = std::make_unique<Writer>(name, O.hdf5_single, O.hdf5_compression, true);
    const hid_t root = W->file();
    Writer::attr(root, "openPMD", std::string("1.1.0"));
    Writer::attrU32(root, "openPMDextension", 0u);
    Writer::attr(root, "basePath", std::string("/data/%T/"));
    Writer::attr(root, "meshesPath", std::string("meshes/"));
    Writer::attr(root, "particlesPath", std::string("particles/"));
    Writer::attr(root, "iterationEncoding", std::string("fileBased"));
    Writer::attr(root, "iterationFormat", std::string("openpmd_%T.h5"));
    Writer::attr(root, "software", std::string("GEM-PIC-RZ"));
    Writer::attr(root, "softwareVersion", std::string(GEMPIC_VERSION));
    Writer::attr(root, "author", std::string("GEM-PIC-RZ user"));
    {
      char buf[64];
      const std::time_t t = std::time(nullptr);
      std::strftime(buf, sizeof buf, "%Y-%m-%d %H:%M:%S %z", std::localtime(&t));
      Writer::attr(root, "date", std::string(buf));
    }
    Writer::attr(root, "comment", std::string("GEM-PIC-RZ (Galilean GT2 scheme): fields of the window at z = zeta + s, "
                                              "t = s / c; input deck in /data/<T>/input_deck"));
    base = "/data/" + std::to_string(S.step_no);
    {
      Handle it(W->groups(base), H5Gclose, "group " + base);
      Writer::attr(it, "time", S.s0);
      Writer::attr(it, "dt", S.Delta);
      Writer::attr(it, "timeUnitSI", si.T0);
    }
    W->text(base + "/input_deck", S.deckText);
    // ---- output grid
    dr = O.openpmd_dr > 0 ? O.openpmd_dr : 1e300;
    dzo = O.openpmd_dzeta > 0 ? O.openpmd_dzeta : 1e300;
    if (!(O.openpmd_dr > 0)) for (int i = 0; i < N; ++i) dr = std::min(dr, S.g.rho[i + 1] - S.g.rho[i]);
    if (!(O.openpmd_dzeta > 0)) for (int n = 0; n < Nz; ++n) dzo = std::min(dzo, S.zv.len(n));
    nr = int(std::floor(S.g.rho[N] / dr + 1e-9)) + 1;
    nz = int(std::floor(S.zv.node(Nz) / dzo + 1e-9)) + 1;
    // radial stencils (the gather of the code, see tools/h5_extract.cpp)
    const RadialGrid& g = S.g;
    rs.resize(nr);
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
    // zeta stencils: slice n and weight; owner half-cell of the density rows
    zn.resize(nz); zw.resize(nz); zb.resize(nz);
    zc.resize(Nz);
    for (int kk = 0; kk < Nz; ++kk) zc[kk] = S.zv.centre(kk);
    for (int k = 0; k < nz; ++k) {
      const double z = k * dzo;
      double wdum = 0.0;
      int n = std::max(0, std::min(Nz - 1, S.zv.cell(z, wdum)));
      while (n > 0 && S.zv.node(n) > z) --n;
      while (n < Nz - 1 && S.zv.node(n + 1) <= z) ++n;
      zn[k] = n;
      zw[k] = std::min(1.0, std::max(0.0, (z - S.zv.node(n)) / S.zv.len(n)));
      const int b = int(std::upper_bound(zc.begin(), zc.end(), z) - zc.begin()) - 1;
      zb[k] = std::min(std::max(b, 0), Nz - 1);
    }
    da.resize(nr); dw.resize(nr);
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
      da[j] = a; dw[j] = w;
    }
    mp = base + "/meshes";
    { Handle h(W->groups(mp), H5Gclose, "group " + mp); }
    const std::vector<std::string> axes = {"r", "z"};
    auto meshAttrs = [&](hid_t obj, int modes, const std::vector<double>& dim) {
      Writer::attr(obj, "geometry", std::string("thetaMode"));
      Writer::attr(obj, "geometryParameters", "m=" + std::to_string(modes) + ";imag=+");
      Writer::attr(obj, "dataOrder", std::string("C"));
      Writer::attr(obj, "axisLabels", axes);
      Writer::attr(obj, "gridSpacing", std::vector<double>{dr, dzo});
      Writer::attr(obj, "gridGlobalOffset", std::vector<double>{0.0, S.s0});
      Writer::attr(obj, "gridUnitSI", si.L);
      Writer::attr(obj, "unitDimension", dim);
      Writer::attr(obj, "timeOffset", 0.0);
    };
    if (O.openpmd_uniform) {
      const char* comp[3] = {"r", "t", "z"};
      for (int vec = 0; vec < 2; ++vec) {
        const std::string rec = mp + (vec == 0 ? "/E" : "/B");
        {
          Handle h(W->groups(rec), H5Gclose, "group " + rec);
          meshAttrs(h, M + 1, vec == 0 ? unitDim(1, 1, -3, -1) : unitDim(0, 1, -2, -1));
        }
        for (int cpt = 0; cpt < 3; ++cpt) {
          const hid_t d = keep(W->empty(rec + "/" + comp[cpt], {hsize_t(NM), hsize_t(nr), hsize_t(nz)}, false));
          Writer::attr(d, "unitSI", vec == 0 ? si.E0 : si.B0);
          Writer::attr(d, "position", std::vector<double>{0.0, 0.0});
          dsU[3 * vec + cpt] = d;
        }
      }
      for (int q = 0; q < (S.cfg.ions.empty() ? 1 : 2); ++q) {
        const hid_t d = keep(W->empty(mp + (q == 0 ? "/n_e" : "/Zn_i"), {1, hsize_t(nr), hsize_t(nz)}, false));
        meshAttrs(d, 1, unitDim(-3, 0, 0, 0));
        Writer::attr(d, "unitSI", S.n_ref * si.nc);
        Writer::attr(d, "position", std::vector<double>{0.0, 0.0});
        dsUn[q] = d;
      }
    }
    if (O.openpmd_native) {
      const std::vector<double> rf(g.rho.begin(), g.rho.begin() + N + 1), rcn(g.c.begin(), g.c.begin() + N);
      std::vector<double> zsl(Nz + 1), zh(Nz);
      for (int n = 0; n <= Nz; ++n) zsl[n] = S.zv.node(n);
      for (int k = 0; k < Nz; ++k) zh[k] = S.zv.centre(k);
      auto minStep = [](const std::vector<double>& x) {
        double d = 1e300;
        for (size_t k = 1; k < x.size(); ++k) d = std::min(d, x[k] - x[k - 1]);
        return x.size() > 1 ? d : 1.0;
      };
      const std::string gc = base + "/gridCoordinates";
      { Handle h(W->groups(gc), H5Gclose, gc); }
      auto coords = [&](const char* nm, const std::vector<double>& x) {
        Handle h(W->plain(gc + "/" + nm, {hsize_t(x.size())}, x.data(), true), H5Dclose, nm);
        return gc + "/" + nm;
      };
      const std::string pRf = coords("r_faces", rf), pRc = coords("r_centers", rcn);
      const std::string pZs = coords("zeta_slices", zsl), pZh = coords("zeta_halfcells", zh);
      auto nativeMesh = [&](const std::string& nm, int modes, const std::vector<double>& rr,
                            const std::vector<double>& zz, const std::string& pr, const std::string& pz,
                            const std::vector<double>& dim, double unitSI) {
        const hid_t d = keep(W->empty(mp + "/" + nm, {hsize_t(2 * modes - 1), hsize_t(rr.size()), hsize_t(zz.size())},
                                      false));
        Writer::attr(d, "geometry", std::string("other"));
        Writer::attr(d, "geometryParameters", "thetaMode;m=" + std::to_string(modes) + ";imag=+");
        Writer::attr(d, "gridType", std::string("rectilinear"));
        Writer::attr(d, "dataOrder", std::string("C"));
        Writer::attr(d, "axisLabels", axes);
        Writer::attr(d, "gridSpacing", std::vector<double>{minStep(rr), minStep(zz)});
        Writer::attr(d, "gridGlobalOffset", std::vector<double>{0.0, S.s0});
        Writer::attr(d, "gridUnitSI", si.L);
        Writer::attr(d, "gridCoordinates", std::vector<std::string>{pr, pz});
        Writer::attr(d, "unitDimension", dim);
        Writer::attr(d, "timeOffset", 0.0);
        Writer::attr(d, "position", std::vector<double>{0.0, 0.0});
        Writer::attr(d, "unitSI", unitSI);
        return d;
      };
      const char* nname[6] = {"native_E_r", "native_E_t", "native_E_z", "native_B_r", "native_B_t", "native_B_z"};
      const bool onFaces[6] = {true, false, false, false, true, true};
      for (int q = 0; q < 6; ++q)
        dsN[q] = nativeMesh(nname[q], M + 1, onFaces[q] ? rf : rcn, zsl, onFaces[q] ? pRf : pRc, pZs,
                            q < 3 ? unitDim(1, 1, -3, -1) : unitDim(0, 1, -2, -1), q < 3 ? si.E0 : si.B0);
      dsNn[0] = nativeMesh("native_n_e", 1, rcn, zh, pRc, pZh, unitDim(-3, 0, 0, 0), S.n_ref * si.nc);
      if (!S.cfg.ions.empty())
        dsNn[1] = nativeMesh("native_Zn_i", 1, rcn, zh, pRc, pZh, unitDim(-3, 0, 0, 0), S.n_ref * si.nc);
    }
  }
  ~OpenPMDSink() override { closeAll(); }
  void closeAll() {
    for (hid_t h : ds) H5Dclose(h);
    ds.clear();
  }
  // complex mode m of quantity q at (slice n, stencil t); q: 0 Er 1 Et 2 Ez 3 Br 4 Bt 5 Bz
  cd Ephi(int m, int n, int i) const { const Fields& F = S.fold; return 0.5 * (F.Rp[F.c(m, n, i)] + F.Lp[F.c(m, n, i)]); }
  cd Brad(int m, int n, int i) const { const Fields& F = S.fold; return 0.5 * (F.Lp[F.c(m, n, i)] - F.Rp[F.c(m, n, i)]); }
  cd Erad(int m, int n, int i) const { const Fields& F = S.fold; return 0.5 * (F.Rr[F.f(m, n, i)] + F.Lr[F.f(m, n, i)]); }
  cd Bphi(int m, int n, int i) const { const Fields& F = S.fold; return 0.5 * (F.Rr[F.f(m, n, i)] - F.Lr[F.f(m, n, i)]); }
  cd Ezc(int m, int n, int i) const { return S.fold.Ez[S.fold.c(m, n, i)]; }
  cd Bzf(int m, int n, int i) const { return S.fold.Bz[S.fold.f(m, n, i)]; }
  cd lin0(cd f0, cd f1) const { return (S.g.c[1] * f0 - S.g.c[0] * f1) / (S.g.c[1] - S.g.c[0]); }
  cd mode(int q, int m, int n, const Rst& t) const {
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
  }
  cd native(int q, int m, int n, int i) const {
    switch (q) {
      case 0: return i > 0 ? Erad(m, n, i) : (m == 1 && N > 1 ? -I * lin0(Ephi(m, n, 0), Ephi(m, n, 1)) : cd(0.0));
      case 1: return Ephi(m, n, i);
      case 2: return Ezc(m, n, i);
      case 3: return Brad(m, n, i);
      case 4: return i > 0 ? Bphi(m, n, i) : (m == 1 && N > 1 ? I * lin0(Brad(m, n, 0), Brad(m, n, 1)) : cd(0.0));
      default: return Bzf(m, n, i);
    }
  }
  void block(int n0, int n1) override {
    const auto& O = S.cfg.output;
    if (O.openpmd_uniform) {
      // field rows k whose slice stencil starts in the block
      const int kb = int(std::lower_bound(zn.begin(), zn.end(), n0) - zn.begin());
      const int ke = int(std::lower_bound(zn.begin(), zn.end(), n1) - zn.begin());
      const int kc = ke - kb;
      if (kc > 0) {
        std::vector<double> out(size_t(NM) * nr * kc);
        std::vector<cd> s0v(nr), s1v(nr);
        for (int q = 0; q < 6; ++q) {
          for (int m = 0; m <= M; ++m) {
            int last = -1;
            for (int k = kb; k < ke; ++k) {
              const int n = zn[k];
              if (n != last) {  // the two slices of the stencil, at all output radii
                for (int j = 0; j < nr; ++j) { s0v[j] = mode(q, m, n, rs[j]); s1v[j] = mode(q, m, n + 1, rs[j]); }
                last = n;
              }
              for (int j = 0; j < nr; ++j) {
                const cd v = (1.0 - zw[k]) * s0v[j] + zw[k] * s1v[j];
                if (m == 0) out[(size_t(0) * nr + j) * kc + k - kb] = v.real();
                else {
                  out[(size_t(2 * m - 1) * nr + j) * kc + k - kb] = 2.0 * v.real();
                  out[(size_t(2 * m) * nr + j) * kc + k - kb] = -2.0 * v.imag();
                }
              }
            }
          }
          Writer::slab(dsU[q], {0, 0, hsize_t(kb)}, {hsize_t(NM), hsize_t(nr), hsize_t(kc)}, out.data(), false);
        }
      }
      // density rows k whose half-cell stencil starts in the block (densities m = 0, half-cells x centres)
      const int db = int(std::lower_bound(zb.begin(), zb.end(), n0) - zb.begin());
      const int de = int(std::lower_bound(zb.begin(), zb.end(), n1) - zb.begin());
      const int dc = de - db;
      for (int q = 0; q < 2 && dc > 0; ++q) {
        if (dsUn[q] < 0) continue;
        const bool ions = (q == 1);
        std::vector<double> d(size_t(nr) * dc);
        auto val = [&](int kk, int i) {
          return ions ? (S.c0s.rhoi[S.c0s.d(kk, i)] + S.c0b.rhoi[S.c0b.d(kk, i)])
                      : -(S.c0s.rho[S.c0s.d(kk, i)] + S.c0b.rho[S.c0b.d(kk, i)]);
        };
        for (int j = 0; j < nr; ++j) {
          const int a = da[j];
          const double w = dw[j];
          auto rad = [&](int kk) { return N < 2 ? val(kk, 0) : (1.0 - w) * val(kk, a) + w * val(kk, a + 1); };
          for (int k = db; k < de; ++k) {
            const double z = k * dzo;
            const int b = int(std::upper_bound(zc.begin(), zc.end(), z) - zc.begin()) - 1;
            double v;
            if (b < 0) v = rad(0);
            else if (b >= Nz - 1) v = rad(Nz - 1);
            else {
              const double wz = (z - zc[b]) / (zc[b + 1] - zc[b]);
              v = (1.0 - wz) * rad(b) + wz * rad(b + 1);
            }
            d[size_t(j) * dc + k - db] = v;
          }
        }
        Writer::slab(dsUn[q], {0, 0, hsize_t(db)}, {1, hsize_t(nr), hsize_t(dc)}, d.data(), false);
      }
    }
    if (O.openpmd_native) {
      const int ne = (n1 == Nz) ? Nz + 1 : n1, nb = ne - n0;
      const bool onFaces[6] = {true, false, false, false, true, true};
      for (int q = 0; q < 6; ++q) {
        const size_t NR = onFaces[q] ? size_t(N) + 1 : size_t(N);
        std::vector<double> d(size_t(NM) * NR * nb);
        for (int m = 0; m <= M; ++m)
          for (size_t i = 0; i < NR; ++i)
            for (int n = n0; n < ne; ++n) {
              const cd v = native(q, m, n, int(i));
              const size_t c = size_t(n - n0);
              if (m == 0) d[i * nb + c] = v.real();
              else {
                d[(size_t(2 * m - 1) * NR + i) * nb + c] = 2.0 * v.real();
                d[(size_t(2 * m) * NR + i) * nb + c] = -2.0 * v.imag();
              }
            }
        Writer::slab(dsN[q], {0, 0, hsize_t(n0)}, {hsize_t(NM), hsize_t(NR), hsize_t(nb)}, d.data(), false);
      }
      const int kb = n1 - n0;
      std::vector<double> d(size_t(N) * kb);
      for (int q = 0; q < 2; ++q) {
        if (dsNn[q] < 0) continue;
        for (int i = 0; i < N; ++i)
          for (int k = n0; k < n1; ++k)
            d[size_t(i) * kb + k - n0] = q == 0 ? -(S.c0s.rho[S.c0s.d(k, i)] + S.c0b.rho[S.c0b.d(k, i)]) / S.n_ref
                                                : (S.c0s.rhoi[S.c0s.d(k, i)] + S.c0b.rhoi[S.c0b.d(k, i)]) / S.n_ref;
        Writer::slab(dsNn[q], {0, 0, hsize_t(n0)}, {1, hsize_t(N), hsize_t(kb)}, d.data(), false);
      }
    }
  }
  void particles() {
    const auto& O = S.cfg.output;
    Writer& Wr = *W;
    const hid_t root = Wr.file();
    const double L = si.L, Q0 = si.Q0, ME = ME_KG;
    (void)Q0;
    if (!O.particles) { H5Adelete(root, "particlesPath"); return; }
    struct Sp { std::string name; int id; };
    std::vector<Sp> species;
    for (size_t b = 0; b < S.cfg.bunches.size(); ++b) species.push_back({S.cfg.bunches[b].name, int(b)});
    species.push_back({"trapped", SPECIES_TRAPPED});
    for (size_t s = 0; s < S.cfg.ions.size(); ++s)
      species.push_back({"inj_" + S.cfg.ions[s].name, SPECIES_TRAPPED_IONIZED_BASE - int(s)});
    const std::string pp = base + "/particles";
    bool first = true;
    for (const auto& Sp_ : species) {
      std::vector<const Particle*> ps;
      for (const auto& p : S.beam)
        if (p.sp == Sp_.id) ps.push_back(&p);
      if (ps.empty() && Sp_.id < 0) continue;  // no trapped particles: no species
      if (first) { Handle h(Wr.groups(pp), H5Gclose, "group " + pp); first = false; }
      const std::string sg = pp + "/" + Sp_.name;
      { Handle h(Wr.groups(sg), H5Gclose, "group " + sg); }
      const hsize_t np = ps.size();
      std::vector<double> col(np);
      auto record = [&](const std::string& rec, const std::vector<double>& dim, unsigned macroWeighted, double wPower) {
        Handle h(Wr.groups(sg + "/" + rec), H5Gclose, rec);
        Writer::attr(h, "unitDimension", dim);
        Writer::attr(h, "timeOffset", 0.0);
        Writer::attrU32(h, "macroWeighted", macroWeighted);
        Writer::attr(h, "weightingPower", wPower);
      };
      auto component = [&](const std::string& path, double (*get)(const Particle&), double unitSI) {
        for (size_t k = 0; k < np; ++k) col[k] = get(*ps[k]);
        Handle ds(Wr.plain(sg + "/" + path, {np}, col.data()), H5Dclose, path);
        Writer::attr(ds, "unitSI", unitSI);
      };
      auto constant = [&](const std::string& path, double value, double unitSI) {
        Handle h(Wr.groups(sg + "/" + path), H5Gclose, path);
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
      constant("positionOffset/z", S.s0, L);
      record("momentum", unitDim(1, 1, -1, 0), 0u, 1.0);
      component("momentum/x", [](const Particle& p) { return p.ux; }, ME * C_LIGHT);
      component("momentum/y", [](const Particle& p) { return p.uy; }, ME * C_LIGHT);
      component("momentum/z", [](const Particle& p) { return p.uz; }, ME * C_LIGHT);
      // real particles per macroparticle: |q| Q0 / e (all S.beam particles are electrons or positrons)
      {
        for (size_t k = 0; k < np; ++k) col[k] = std::fabs(ps[k]->q) * Q0 / E_CHARGE_C;
        Handle ds(Wr.plain(sg + "/weighting", {np}, col.data()), H5Dclose, "weighting");
        Writer::attr(ds, "unitDimension", unitDim(0, 0, 0, 0));
        Writer::attr(ds, "timeOffset", 0.0);
        Writer::attrU32(ds, "macroWeighted", 1u);
        Writer::attr(ds, "weightingPower", 1.0);
        Writer::attr(ds, "unitSI", 1.0);
      }
      const double qsign = (Sp_.id >= 0) ? (S.cfg.bunches[size_t(Sp_.id)].qm < 0 ? -1.0 : 1.0) : -1.0;
      {
        Handle h(Wr.groups(sg + "/charge"), H5Gclose, "charge");
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
          const double v[3] = {ps[k]->x, ps[k]->y, ps[k]->zeta + S.s0};
          for (int a = 0; a < 3; ++a) {
            if (k == 0 || v[a] < lo[a]) lo[a] = v[a];
            if (k == 0 || v[a] > hi[a]) hi[a] = v[a];
          }
        }
        const std::string pg = sg + "/particlePatches";
        { Handle h(Wr.groups(pg), H5Gclose, pg); }
        const unsigned long long cnt = np, off0 = 0;
        Wr.u64(pg + "/numParticles", {1}, &cnt, "particles in the patch");
        Wr.u64(pg + "/numParticlesOffset", {1}, &off0, "offset of the patch in the particle records");
        for (int which = 0; which < 2; ++which) {
          const std::string rec = pg + (which == 0 ? "/offset" : "/extent");
          {
            Handle h(Wr.groups(rec), H5Gclose, rec);
            Writer::attr(h, "unitDimension", unitDim(1, 0, 0, 0));
          }
          const char* ax[3] = {"x", "y", "z"};
          for (int a = 0; a < 3; ++a) {
            const double v = which == 0 ? lo[a] : hi[a] - lo[a];
            Handle ds(Wr.plain(rec + "/" + ax[a], {1}, &v), H5Dclose, rec);
            Writer::attr(ds, "unitSI", L);
          }
        }
      }
      {
        Handle h(Wr.groups(sg + "/mass"), H5Gclose, "mass");
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
  std::string finish() override {
    closeAll();
    particles();
    W.reset();
    return "";
  }
};

// ---------------------------------------------------------------- LASY envelope file (see laser_file.hpp)
// The forward field R/2 at the output radii is collected block by block (R_r,m and R_phi,m at every time sample
// and output radius: the demodulation needs the whole time series of a radius); the rest at the end.
struct Simulation::LasySink : Simulation::Sink {
  Simulation& S;
  std::string path;
  double Lm, E0, T0, wL, dr = 0, dt = 0, tmin = 0;
  int nr = 0, nt = 0;
  std::vector<int> zn;
  std::vector<double> zw;
  std::vector<cd> RR, RP;  // [(m * nt + l) * nr + j]
  explicit LasySink(Simulation& s) : S(s) {
    namespace fs = std::filesystem;
    const auto& O = S.cfg.output;
    const std::string dir = O.dir + "/lasy";
    fs::create_directories(dir);
    char nm[64];
    std::snprintf(nm, sizeof nm, "laser_%05d.h5", S.step_no);
    path = dir + "/" + nm;
    const Units& u = S.cfg.units;
    Lm = u.k0inv_um() * 1e-6;
    E0 = MC2_MEV * 1e6 / Lm;
    T0 = 1.0 / (u.omega0_per_fs() * 1e15);
    wL = O.lasy_omega;
    const int N = S.g.N, Nz = S.Nz, M = S.M;
    const RadialGrid& g = S.g;
    dr = O.lasy_dr;
    dt = O.lasy_dt;
    if (!(dr > 0)) { dr = 1e300; for (int i = 0; i < N; ++i) dr = std::min(dr, g.rho[i + 1] - g.rho[i]); }
    if (!(dt > 0)) { dt = 1e300; for (int n = 0; n < Nz; ++n) dt = std::min(dt, S.zv.len(n)); }
    double rmax = O.lasy_r_max > 0 ? O.lasy_r_max : (S.cfg.grid.r_absorber < g.rho[N] ? S.cfg.grid.r_absorber : g.rho[N]);
    rmax = std::min(rmax, g.rho[N]);
    nr = int(std::floor(rmax / dr + 1e-9)) + 1;
    const double Z = S.zv.node(Nz - 1);
    nt = int(std::floor(Z / dt + 1e-9)) + 1;
    tmin = S.s0 + S.lasyTref - Z;  // t_l = tmin + l dt  <->  zeta_l = Z - l dt
    zn.resize(nt);
    zw.resize(nt);
    for (int l = 0; l < nt; ++l) {
      const double z = std::max(0.0, Z - l * dt);
      double wdum = 0.0;
      int n = std::max(0, std::min(Nz - 2, S.zv.cell(z, wdum)));
      while (n > 0 && S.zv.node(n) > z) --n;
      while (n < Nz - 2 && S.zv.node(n + 1) <= z) ++n;
      zn[l] = n;
      zw[l] = std::min(1.0, std::max(0.0, (z - S.zv.node(n)) / S.zv.len(n)));
    }
    RR.assign(size_t(M + 1) * nt * nr, cd(0.0));
    RP.assign(size_t(M + 1) * nt * nr, cd(0.0));
  }
  // R_r,m (faces) and R_phi,m (centres) of slice n at radius r, with the regular forms on the axis
  cd lin0(cd f0, cd f1) const { return (S.g.c[1] * f0 - S.g.c[0] * f1) / (S.g.c[1] - S.g.c[0]); }
  cd Rphi(int m, int n, double r) const {
    const RadialGrid& g = S.g;
    const int N = g.N;
    auto F = [&](int i) { return S.fold.Rp[S.fold.c(m, n, i)]; };
    if (r >= g.c[N - 1]) return F(N - 1);
    if (r < g.c[0]) {
      const double x = r / g.c[0];
      if (m == 0) return F(0) * x;
      if (m == 1) { const cd a0 = lin0(F(0), F(1)); return a0 + (F(0) - a0) * x; }
      return F(0) * std::pow(x, m - 1);
    }
    const int a = int(std::upper_bound(g.c.begin(), g.c.begin() + N, r) - g.c.begin()) - 1;
    const double w = (r - g.c[a]) / (g.c[a + 1] - g.c[a]);
    return (1 - w) * F(a) + w * F(a + 1);
  }
  cd Rrad(int m, int n, double r) const {
    const RadialGrid& g = S.g;
    const int N = g.N;
    const cd I(0.0, 1.0);
    auto F = [&](int i) -> cd {
      if (i > 0) return S.fold.Rr[S.fold.f(m, n, i)];
      return m == 1 ? -I * lin0(S.fold.Rp[S.fold.c(1, n, 0)], S.fold.Rp[S.fold.c(1, n, 1)]) : cd(0.0);
    };
    if (r >= g.rho[N]) return F(N);
    const int a = int(std::upper_bound(g.rho.begin(), g.rho.begin() + N + 1, r) - g.rho.begin()) - 1;
    const double w = (r - g.rho[a]) / (g.rho[a + 1] - g.rho[a]);
    return (1 - w) * F(a) + w * F(a + 1);
  }
  void block(int n0, int n1) override {
    const int M = S.M;
#pragma omp parallel for schedule(static)
    for (int l = 0; l < nt; ++l) {
      const int n = zn[l];
      if (n < n0 || n >= n1) continue;
      const double w = zw[l];
      for (int j = 0; j < nr; ++j) {
        const double r = j * dr;
        for (int m = 0; m <= M; ++m) {
          RR[(size_t(m) * nt + l) * nr + j] = (1 - w) * Rrad(m, n, r) + w * Rrad(m, n + 1, r);
          RP[(size_t(m) * nt + l) * nr + j] = (1 - w) * Rphi(m, n, r) + w * Rphi(m, n + 1, r);
        }
      }
    }
  }
  std::string finish() override {
    const auto& O = S.cfg.output;
    const int M = S.M, K = M + 1, NK = 2 * K + 1, nth = 2 * M + 4;
    const cd I(0.0, 1.0);
    H5Eset_auto2(H5E_DEFAULT, nullptr, nullptr);
    // ---- per radius: E_x, E_y at the angles, envelopes, azimuthal modes b_k (x and y)
    int n2 = 1;
    while (n2 < 2 * nt) n2 <<= 1;
    const double dOm = 2.0 * PI / (n2 * dt);
    std::vector<cd> bx(size_t(NK) * nt * nr), by(size_t(NK) * nt * nr);  // [k + K][l][j]
  #pragma omp parallel for schedule(dynamic)
    for (int j = 0; j < nr; ++j) {
      const double r = j * dr;
      std::vector<cd> rr(size_t(M + 1) * nt), rp(size_t(M + 1) * nt), X(n2), Y(n2);
      (void)r;
      for (int m = 0; m <= M; ++m)
        for (int l = 0; l < nt; ++l) {
          rr[size_t(m) * nt + l] = RR[(size_t(m) * nt + l) * nr + j];
          rp[size_t(m) * nt + l] = RP[(size_t(m) * nt + l) * nr + j];
        }
      for (int q = 0; q < nth; ++q) {
        const double th = 2.0 * PI * q / nth, ct = std::cos(th), st = std::sin(th);
        std::fill(X.begin(), X.end(), cd(0.0));
        std::fill(Y.begin(), Y.end(), cd(0.0));
        for (int l = 0; l < nt; ++l) {
          // E = R_0 / 2 + Re sum_m R_m e^{i m theta} (forward field E = R / 2)
          double er = 0.5 * rr[l].real(), ep = 0.5 * rp[l].real();
          for (int m = 1; m <= M; ++m) {
            const cd e = std::polar(1.0, m * th);
            er += (rr[size_t(m) * nt + l] * e).real();
            ep += (rp[size_t(m) * nt + l] * e).real();
          }
          X[l] = er * ct - ep * st;
          Y[l] = er * st + ep * ct;
        }
        // analytic signal in the band around omega_L (e^{-i omega t}, omega > 0), then the envelope
        for (auto* V : {&X, &Y}) {
          fft(*V, +1);
          for (int k = 0; k < n2; ++k) {
            const double om = (k < n2 / 2 ? k : k - n2) * dOm;
            if (!(om > 0 && std::fabs(om - wL) < O.lasy_band * wL)) (*V)[k] = 0.0;
          }
          fft(*V, -1);
        }
        for (int l = 0; l < nt; ++l) {
          const cd c = std::polar(2.0 / n2, wL * (tmin + l * dt));
          const cd ex = X[l] * c, ey = Y[l] * c;
          for (int k = -K; k <= K; ++k) {  // env = sum_k b_k e^{-i k theta}
            const cd ph = std::polar(1.0 / nth, k * th);
            bx[(size_t(k + K) * nt + l) * nr + j] += ex * ph;
            by[(size_t(k + K) * nt + l) * nr + j] += ey * ph;
          }
        }
      }
    }
    // ---- dominant polarization: coherence matrix (Parseval over the modes, weight r)
    double jxx = 0, jyy = 0;
    cd jxy = 0.0;
    for (int k = 0; k < NK; ++k)
      for (int l = 0; l < nt; ++l)
        for (int j = 0; j < nr; ++j) {
          const size_t id = (size_t(k) * nt + l) * nr + j;
          const double w = std::max(j * dr, 0.25 * dr);
          jxx += std::norm(bx[id]) * w;
          jyy += std::norm(by[id]) * w;
          jxy += std::conj(bx[id]) * by[id] * w;
        }
    const double tr = jxx + jyy, det = jxx * jyy - std::norm(jxy);
    const double lmax = 0.5 * tr + std::sqrt(std::max(0.0, 0.25 * tr * tr - det));
    cd px = 1.0, py = 0.0;
    if (tr > 0) {
      // eigenvector of [[jxx, jxy], [conj jxy, jyy]] for lmax
      cd vx = jxy, vy = lmax - jxx;
      if (std::norm(vx) + std::norm(vy) < 1e-24 * tr * tr) { vx = lmax - jyy; vy = std::conj(jxy); }
      if (std::norm(vx) + std::norm(vy) < 1e-24 * tr * tr) { vx = jxx >= jyy ? 1.0 : 0.0; vy = jxx >= jyy ? 0.0 : 1.0; }
      const double nn = std::sqrt(std::norm(vx) + std::norm(vy));
      px = vx / nn;
      py = vy / nn;
      const cd ref = std::abs(px) > 1e-6 ? px : py;  // phase convention: the first nonzero component real > 0
      const cd ph = std::conj(ref) / std::abs(ref);
      px *= ph;
      py *= ph;
      if (std::fabs(px.imag()) < 1e-12) px = px.real();
      if (std::fabs(py.imag()) < 1e-12 * std::max(1.0, std::abs(py))) py = py.real();
    }
    const double captured = tr > 0 ? lmax / tr : 1.0;
    // ---- envelope modes -> openPMD components [2K + 1][t][r] in V/m
    std::vector<cd> out(size_t(NK) * nt * nr);
    for (int l = 0; l < nt; ++l)
      for (int j = 0; j < nr; ++j) {
        auto e = [&](int k) {
          const size_t id = (size_t(k + K) * nt + l) * nr + j;
          return (std::conj(px) * bx[id] + std::conj(py) * by[id]) * E0;
        };
        out[(size_t(0) * nt + l) * nr + j] = e(0);
        for (int m = 1; m <= K; ++m) {
          out[(size_t(2 * m - 1) * nt + l) * nr + j] = e(m) + e(-m);
          out[(size_t(2 * m) * nt + l) * nr + j] = -I * e(m) + I * e(-m);
        }
      }

    // ---- file (layout of LASY's writer)
    Writer W(path, false, O.hdf5_compression);
    const hid_t root = W.file();
    Writer::attr(root, "openPMD", std::string("1.1.0"));
    Writer::attrU32(root, "openPMDextension", 0u);
    Writer::attr(root, "basePath", std::string("/data/%T/"));
    Writer::attr(root, "meshesPath", std::string("meshes/"));
    Writer::attr(root, "iterationEncoding", std::string("fileBased"));
    Writer::attr(root, "iterationFormat", std::string("laser_%05T"));
    Writer::attr(root, "software", std::string("GEM-PIC-RZ"));
    Writer::attr(root, "softwareVersion", std::string(GEMPIC_VERSION));
    {
      char buf[64];
      const std::time_t t = std::time(nullptr);
      std::strftime(buf, sizeof buf, "%Y-%m-%d %H:%M:%S %z", std::localtime(&t));
      Writer::attr(root, "date", std::string(buf));
    }
    char cm[512];
    std::snprintf(cm, sizeof cm,
                  "GEM-PIC-RZ laser envelope at the plane z = %.9g m of the simulation (step %d): forward field R/2, "
                  "polarization captures %.6f of the energy",
                  S.s0 * Lm, S.step_no, captured);
    Writer::attr(root, "comment", std::string(cm));
    const std::string base = "/data/" + std::to_string(S.step_no);
    {
      Handle it(W.groups(base), H5Gclose, "group " + base);
      Writer::attr(it, "time", 0.0);
      Writer::attr(it, "dt", 1.0);
      Writer::attr(it, "timeUnitSI", 1.0);
    }
    { Handle h(W.groups(base + "/meshes"), H5Gclose, "group meshes"); }
    const std::string mp = base + "/meshes/laserEnvelope";
    W.complex(mp, {hsize_t(NK), hsize_t(nt), hsize_t(nr)}, out.data(), 1.0, "V/m",
              "complex envelope of the transverse electric field, openPMD thetaMode components");
    Handle ds(H5Dopen2(W.file(), mp.c_str(), H5P_DEFAULT), H5Dclose, "laserEnvelope");
    Writer::attr(ds, "geometry", std::string("thetaMode"));
    Writer::attr(ds, "geometryParameters", "m=" + std::to_string(K + 1) + ";imag=+");
    Writer::attr(ds, "dataOrder", std::string("C"));
    Writer::attr(ds, "axisLabels", std::vector<std::string>{"t", "r"});
    Writer::attr(ds, "gridSpacing", std::vector<double>{dt * T0, dr * Lm});
    Writer::attr(ds, "gridGlobalOffset", std::vector<double>{tmin * T0, 0.0});
    Writer::attr(ds, "gridUnitSI", 1.0);
    Writer::attr(ds, "position", std::vector<double>{0.0, 0.0});
    Writer::attr(ds, "timeOffset", 0.0);
    Writer::attr(ds, "unitDimension", std::vector<double>{1, 1, -3, -1, 0, 0, 0});
    Writer::attr(ds, "angularFrequency", wL / T0);
    Writer::attr(ds, "envelopeField", std::string("electric_field"));
    if (px.imag() == 0.0 && py.imag() == 0.0) {
      Writer::attr(ds, "polarization", std::vector<double>{px.real(), py.real()});
    } else {
      const hsize_t two = 2;
      Handle sp(H5Screate_simple(1, &two, nullptr), H5Sclose, "dataspace");
      Handle ct(H5Tcreate(H5T_COMPOUND, 2 * sizeof(double)), H5Tclose, "complex type");
      H5Tinsert(ct, "r", 0, H5T_NATIVE_DOUBLE);
      H5Tinsert(ct, "i", sizeof(double), H5T_NATIVE_DOUBLE);
      Handle a(H5Acreate2(ds, "polarization", ct, sp, H5P_DEFAULT, H5P_DEFAULT), H5Aclose, "polarization");
      const cd pv[2] = {px, py};
      check(H5Awrite(a, ct, pv), "polarization");
    }
    return "LASY file: " + path + "  (polarization (" + std::to_string(px.real()).substr(0, 6) + (px.imag() != 0 ? "..." : "") + ", " +
           std::to_string(py.real()).substr(0, 6) + (py.imag() >= 0 ? "+" : "") + std::to_string(py.imag()).substr(0, 6) +
           "i), captures " + std::to_string(captured).substr(0, 8) + " of the energy)";
  }
};

// ---------------------------------------------------------------- checkpoint_NNNNN.h5
// the complete state in double precision; fields and currents as 1-D datasets [m][slice / half-cell][r]
struct Simulation::CheckpointSink : Simulation::Sink {
  Simulation& S;
  std::string path, tmp;
  std::unique_ptr<Writer> W;
  hid_t dsF[6] = {-1, -1, -1, -1, -1, -1};
  hid_t dsC[2][5] = {{-1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1}};  // beam / stream: jr, jp, jz, rho, rhoi
  explicit CheckpointSink(Simulation& s) : S(s) {
    namespace fs = std::filesystem;
    const auto& O = S.cfg.output;
    fs::create_directories(O.dir);
    char name[1024];
    std::snprintf(name, sizeof name, "%s/checkpoint_%05d.h5", O.dir.c_str(), S.step_no);
    path = name;
    tmp = path + ".tmp";
    H5Eset_auto2(H5E_DEFAULT, nullptr, nullptr);
    W = std::make_unique<Writer>(tmp, /*single=*/false, /*compression=*/0);
    const hid_t root = W->file();
    const int N = S.g.N, M = S.M, Nz = S.Nz;
    Writer::attr(root, "code", std::string("GEM-PIC-RZ"));
    Writer::attr(root, "kind", std::string("checkpoint"));
    Writer::attr(root, "format_version", 1LL);
    Writer::attr(root, "step", (long long)S.step_no);
    Writer::attr(root, "z", S.s0);
    Writer::attr(root, "seed", (long long)S.cfg.numerics.seed);
    Writer::attr(root, "n_restreamed", (long long)S.n_restreamed);
    Writer::attr(root, "n_trapped_total", (long long)S.n_trapped_total);
    Writer::attr(root, "lambda0_um", S.cfg.units.lambda0_um);
    Writer::attr(root, "modes", (long long)(M + 1));
    Writer::attr(root, "Nr", (long long)N);
    Writer::attr(root, "Nzeta", (long long)Nz);
    Writer::attr(root, "dzeta", S.cfg.grid.dzeta);
    Writer::attr(root, "dz", S.cfg.grid.stepAt(S.s0));  // informative: the z grid may change at a restart
    std::string ions;
    for (size_t k = 0; k < S.cfg.ions.size(); ++k) ions += (k ? "," : "") + S.cfg.ions[k].name;
    Writer::attr(root, "ion_names", ions);
    Writer::attr(root, "bunches", (long long)S.cfg.bunches.size());
    Writer::attr(root, "backend", S.backendInfo());
    W->group("/grid");
    W->real("/grid/r_faces", {hsize_t(N + 1)}, S.g.rho.data(), 1.0, "1/k0", "faces");
    W->real("/grid/r_centres", {hsize_t(N)}, S.g.c.data(), 1.0, "1/k0", "centres");
    W->real("/grid/zeta_slices", {hsize_t(Nz + 1)}, S.zg.z.data(), 1.0, "1/k0", "slices");
    W->group("/state");
    W->group("/state/fields");
    const char* desc = "[m][slice][r], level s^0";
    const hsize_t nf = hsize_t(M + 1) * (Nz + 1) * (N + 1), nc = hsize_t(M + 1) * (Nz + 1) * N;
    const char* fn[6] = {"Rr", "Lr", "Bz", "Rp", "Lp", "Ez"};
    for (int q = 0; q < 6; ++q) {
      dsF[q] = W->empty(std::string("/state/fields/") + fn[q], {q < 3 ? nf : nc}, true, true);
      Writer::describeDs(dsF[q], 1.0, "code units", desc);
    }
    const hsize_t cf = hsize_t(M + 1) * Nz * (N + 1), cc = hsize_t(M + 1) * Nz * N, cd_ = hsize_t(Nz) * N;
    const char* gname[2] = {"/state/beam_currents", "/state/stream_currents"};
    for (int b = 0; b < 2; ++b) {
      const std::string g = gname[b];
      W->group(g);
      dsC[b][0] = W->empty(g + "/jr", {cf}, true, true);
      Writer::describeDs(dsC[b][0], 1.0, "code units", "j_r modes [m][half-cell][face]");
      dsC[b][1] = W->empty(g + "/jp", {cc}, true, true);
      Writer::describeDs(dsC[b][1], 1.0, "code units", "j_phi modes [m][half-cell][centre]");
      dsC[b][2] = W->empty(g + "/jz", {cc}, true, true);
      Writer::describeDs(dsC[b][2], 1.0, "code units", "j_z modes [m][half-cell][centre]");
      dsC[b][3] = W->empty(g + "/rho", {cd_}, false, true);
      Writer::describeDs(dsC[b][3], 1.0, "code units", "charge density (electrons, bunches), m = 0");
      dsC[b][4] = W->empty(g + "/rhoi", {cd_}, false, true);
      Writer::describeDs(dsC[b][4], 1.0, "code units", "charge density of the ions, m = 0");
    }
  }
  ~CheckpointSink() override { close(); }
  void close() {
    for (auto& d : dsF) if (d >= 0) { H5Dclose(d); d = -1; }
    for (auto& r : dsC) for (auto& d : r) if (d >= 0) { H5Dclose(d); d = -1; }
  }
  void block(int n0, int n1) override {
    const int N = S.g.N, M = S.M, Nz = S.Nz;
    const int ne = (n1 == Nz) ? Nz + 1 : n1;
    const Fields& F = S.fold;
    const std::vector<cd>* fv[6] = {&F.Rr, &F.Lr, &F.Bz, &F.Rp, &F.Lp, &F.Ez};
    for (int q = 0; q < 6; ++q) {
      const int L = q < 3 ? N + 1 : N;
      for (int m = 0; m <= M; ++m) {
        const size_t j0 = q < 3 ? F.f(m, n0, 0) : F.c(m, n0, 0);
        Writer::slab(dsF[q], {hsize_t((size_t(m) * (Nz + 1) + n0) * L)}, {hsize_t(size_t(ne - n0) * L)},
                     &(*fv[q])[j0], true);
      }
    }
    const Currents* cs[2] = {&S.c0b, &S.c0s};
    for (int b = 0; b < 2; ++b) {
      const Currents& C = *cs[b];
      const std::vector<cd>* cv[3] = {&C.jr, &C.jp, &C.jz};
      for (int q = 0; q < 3; ++q) {
        const int L = q == 0 ? N + 1 : N;
        for (int m = 0; m <= M; ++m) {
          const size_t j0 = q == 0 ? C.f(m, n0, 0) : C.c(m, n0, 0);
          Writer::slab(dsC[b][q], {hsize_t((size_t(m) * Nz + n0) * L)}, {hsize_t(size_t(n1 - n0) * L)},
                       &(*cv[q])[j0], true);
        }
      }
      Writer::slab(dsC[b][3], {hsize_t(size_t(n0) * N)}, {hsize_t(size_t(n1 - n0) * N)}, &C.rho[C.d(n0, 0)], false);
      Writer::slab(dsC[b][4], {hsize_t(size_t(n0) * N)}, {hsize_t(size_t(n1 - n0) * N)}, &C.rhoi[C.d(n0, 0)], false);
    }
  }
  std::string finish() override {
    namespace fs = std::filesystem;
    close();
    writeParticles(*W, "/state/beam", S.beam);
    writeParticles(*W, "/state/pending", S.pendingNext);
    std::ostringstream os;
    os << S.rng;
    W->text("/state/rng", os.str());
    W->text("/input_deck", S.deckText);
    W.reset();
    fs::rename(tmp, path);  // a crash while writing never leaves a broken checkpoint_NNNNN.h5
    const auto& O = S.cfg.output;
    if (O.checkpoint_keep > 0) {  // keep only the newest checkpoint_keep files
      std::vector<std::pair<int, fs::path>> files;
      for (const auto& e : fs::directory_iterator(O.dir)) {
        const int st = checkpointStep(e.path().filename().string());
        if (st >= 0) files.push_back({st, e.path()});
      }
      std::sort(files.begin(), files.end());
      for (size_t k = 0; k + O.checkpoint_keep < files.size(); ++k) fs::remove(files[k].second);
    }
    return "checkpoint: " + path;
  }
};

std::unique_ptr<Simulation::Sink> Simulation::makeHdf5Sink(Simulation& S, int kind) {
  switch (kind) {
    case 0: return std::make_unique<Hdf5Sink>(S);
    case 1: return std::make_unique<OpenPMDSink>(S);
    case 2: return std::make_unique<LasySink>(S);
    default: return std::make_unique<CheckpointSink>(S);
  }
}

// serial run / engine: the whole window as one block
std::string Simulation::writeCheckpoint() {
  if (engine) engine->downloadFull();  // currents and pending particles of the device
  OutputPlan plan;
  plan.checkpoint = true;
  auto out = beginOutput(plan);
  outputBlock(*out, 0, Nz);
  const auto msg = endOutput(*out);
  return msg.empty() ? std::string() : msg.back().substr(12);
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
std::unique_ptr<Simulation::Sink> Simulation::makeHdf5Sink(Simulation&, int) {
  throw std::runtime_error("this build has no HDF5 support");
}
std::string Simulation::writeCheckpoint() { throw std::runtime_error("checkpoints need HDF5 support"); }
void Simulation::readCheckpoint(const std::string&) { throw std::runtime_error("restarts need HDF5 support"); }
std::string Simulation::latestCheckpoint(const std::string&) { return ""; }
}  // namespace gem

#endif
