// ============================================================================
//  output_hdf5_internal.hpp -- HDF5 helpers of the output (src/output_*.cpp): handles, Writer, Reader,
//  particles and currents to and from groups. The file layout is described in src/output_hdf5.cpp.
// ============================================================================
#pragma once
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
namespace h5detail {

// ---------------------------------------------------------------- HDF5 helpers
inline void check(herr_t e, const std::string& what) {
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
inline void writeParticles(Writer& W, const std::string& g, const std::vector<Particle>& ps) {
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

inline std::vector<Particle> readParticles(const Reader& R, const std::string& g) {
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

inline void writeCurrents(Writer& W, const std::string& g, const Currents& C) {
  W.group(g);
  W.complex(g + "/jr", {C.jr.size()}, C.jr.data(), 1.0, "code units", "j_r modes [m][half-cell][face]");
  W.complex(g + "/jp", {C.jp.size()}, C.jp.data(), 1.0, "code units", "j_phi modes [m][half-cell][centre]");
  W.complex(g + "/jz", {C.jz.size()}, C.jz.data(), 1.0, "code units", "j_z modes [m][half-cell][centre]");
  W.real(g + "/rho", {C.rho.size()}, C.rho.data(), 1.0, "code units", "charge density (electrons, bunches), m = 0");
  W.real(g + "/rhoi", {C.rhoi.size()}, C.rhoi.data(), 1.0, "code units", "charge density of the ions, m = 0");
}
inline void readCurrents(const Reader& R, const std::string& g, Currents& C) {
  R.complexes(g + "/jr", C.jr);
  R.complexes(g + "/jp", C.jp);
  R.complexes(g + "/jz", C.jz);
  C.rho = R.reals(g + "/rho", C.rho.size());
  C.rhoi = R.reals(g + "/rhoi", C.rhoi.size());
}

// step number of a checkpoint file name, -1 if it is none
inline int checkpointStep(const std::string& fname) {
  int step = -1;
  char tail[8] = {0};
  if (std::sscanf(fname.c_str(), "checkpoint_%d.h%2s", &step, tail) == 2 && std::string(tail) == "5") return step;
  return -1;
}

}  // namespace h5detail
}  // namespace gem

#endif
