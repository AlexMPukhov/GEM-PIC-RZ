// ============================================================================
//  beam_import.cpp -- read a particle beam from openPMD or native HDF5 files
// ============================================================================
#include "beam_import.hpp"

#include <cmath>
#include <stdexcept>

#if defined(GEMPIC_HAVE_HDF5)
#include <hdf5.h>

#include <algorithm>
#include <map>
#include <sstream>

namespace gem {

namespace {

constexpr double QE = 1.602176634e-19, ME = 9.1093837015e-31, CL = 299792458.0;

struct H {  // owning handle
  hid_t id;
  herr_t (*close)(hid_t);
  H(hid_t i, herr_t (*c)(hid_t)) : id(i), close(c) {}
  H(const H&) = delete;
  ~H() { if (id >= 0) close(id); }
  operator hid_t() const { return id; }
};

class File {
 public:
  explicit File(const std::string& p) : path(p), f(H5Fopen(p.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT), H5Fclose) {
    if (f.id < 0) fail("cannot open the file (not HDF5?)");
  }
  std::string path;
  H f;

  [[noreturn]] void fail(const std::string& m) const { throw std::runtime_error(path + ": " + m); }
  bool exists(const std::string& p) const {
    // check every level (H5Lexists fails on a missing parent)
    size_t q = 1;
    for (;;) {
      q = p.find('/', q);
      const std::string sub = p.substr(0, q);
      if (H5Lexists(f, sub.c_str(), H5P_DEFAULT) <= 0) return false;
      if (q == std::string::npos) return true;
      ++q;
    }
  }
  bool isGroup(const std::string& p) const {
    H5O_info_t info;
#if H5_VERSION_GE(1, 12, 0)
    if (H5Oget_info_by_name3(f, p.c_str(), &info, H5O_INFO_BASIC, H5P_DEFAULT) < 0) return false;
#else
    if (H5Oget_info_by_name(f, p.c_str(), &info, H5P_DEFAULT) < 0) return false;
#endif
    return info.type == H5O_TYPE_GROUP;
  }
  bool hasAttr(const std::string& obj, const char* name) const {
    return H5Aexists_by_name(f, obj.c_str(), name, H5P_DEFAULT) > 0;
  }
  std::vector<double> attrD(const std::string& obj, const char* name) const {
    H a(H5Aopen_by_name(f, obj.c_str(), name, H5P_DEFAULT, H5P_DEFAULT), H5Aclose);
    if (a.id < 0) fail("missing attribute " + obj + "@" + name);
    H sp(H5Aget_space(a), H5Sclose);
    std::vector<double> v(size_t(std::max<hssize_t>(1, H5Sget_simple_extent_npoints(sp))));
    if (H5Aread(a, H5T_NATIVE_DOUBLE, v.data()) < 0) fail("cannot read attribute " + obj + "@" + name);
    return v;
  }
  double attr1(const std::string& obj, const char* name, double def) const {
    return hasAttr(obj, name) ? attrD(obj, name)[0] : def;
  }
  std::string attrS(const std::string& obj, const char* name) const {
    H a(H5Aopen_by_name(f, obj.c_str(), name, H5P_DEFAULT, H5P_DEFAULT), H5Aclose);
    if (a.id < 0) return "";
    H t(H5Aget_type(a), H5Tclose);
    if (H5Tget_class(t) != H5T_STRING) return "";
    if (H5Tis_variable_str(t) > 0) {
      char* s = nullptr;
      H mt(H5Tcopy(H5T_C_S1), H5Tclose);
      H5Tset_size(mt, H5T_VARIABLE);
      if (H5Aread(a, mt, &s) < 0 || !s) return "";
      std::string r(s);
      H5free_memory(s);
      return r;
    }
    std::string s(H5Tget_size(t), '\0');
    H5Aread(a, t, &s[0]);
    return s.c_str();
  }
  std::vector<double> data(const std::string& p) const {
    H d(H5Dopen2(f, p.c_str(), H5P_DEFAULT), H5Dclose);
    if (d.id < 0) fail("no dataset " + p);
    H sp(H5Dget_space(d), H5Sclose);
    std::vector<double> v(size_t(H5Sget_simple_extent_npoints(sp)));
    if (!v.empty() && H5Dread(d, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, v.data()) < 0)
      fail("cannot read " + p);
    return v;
  }
  std::vector<std::string> children(const std::string& p) const {
    std::vector<std::string> out;
    H g(H5Gopen2(f, p.c_str(), H5P_DEFAULT), H5Gclose);
    if (g.id < 0) return out;
    H5G_info_t gi;
    H5Gget_info(g, &gi);
    for (hsize_t k = 0; k < gi.nlinks; ++k) {
      char buf[512];
      if (H5Lget_name_by_idx(g, ".", H5_INDEX_NAME, H5_ITER_INC, k, buf, sizeof buf, H5P_DEFAULT) > 0)
        out.push_back(buf);
    }
    return out;
  }
  // openPMD record component: dataset or constant (group with "value"), times unitSI
  std::vector<double> component(const std::string& p, size_t n, bool required = true) const {
    if (!exists(p)) {
      if (required) fail("missing record component " + p);
      return std::vector<double>(n, 0.0);
    }
    const double u = attr1(p, "unitSI", 1.0);
    std::vector<double> v;
    if (isGroup(p)) {
      v.assign(n, attrD(p, "value")[0]);
    } else {
      v = data(p);
      if (v.size() != n) fail(p + " has " + std::to_string(v.size()) + " entries, expected " + std::to_string(n));
    }
    for (auto& x : v) x *= u;
    return v;
  }
  // number of particles of a record (dataset size or constant "shape")
  size_t count(const std::string& p) const {
    if (isGroup(p)) return size_t(attrD(p, "shape")[0]);
    H d(H5Dopen2(f, p.c_str(), H5P_DEFAULT), H5Dclose);
    if (d.id < 0) fail("no dataset " + p);
    H sp(H5Dget_space(d), H5Sclose);
    return size_t(H5Sget_simple_extent_npoints(sp));
  }
};

std::string join(const std::vector<std::string>& v) {
  std::string s;
  for (const auto& x : v) s += (s.empty() ? "" : ", ") + x;
  return s.empty() ? "none" : s;
}

ImportedBeam readOpenPMD(const File& F, const std::string& species, long iteration) {
  ImportedBeam B;
  B.format = "openPMD " + F.attrS("/", "openPMD");
  std::string base = F.attrS("/", "basePath");
  if (base.empty()) base = "/data/%T/";
  const std::string itRoot = base.substr(0, base.find("%T"));
  std::vector<long> its;
  for (const auto& c : F.children(itRoot.empty() ? "/" : itRoot)) {
    char* end = nullptr;
    const long v = std::strtol(c.c_str(), &end, 10);
    if (end != c.c_str() && *end == '\0') its.push_back(v);
  }
  if (its.empty()) F.fail("no iteration found under " + itRoot);
  std::sort(its.begin(), its.end());
  if (iteration < 0) iteration = its.back();
  if (!std::binary_search(its.begin(), its.end(), iteration)) {
    std::string l;
    for (long v : its) l += (l.empty() ? "" : ", ") + std::to_string(v);
    F.fail("iteration " + std::to_string(iteration) + " not in the file (iterations: " + l + ")");
  }
  B.iteration = iteration;
  std::string pp = F.attrS("/", "particlesPath");
  if (pp.empty()) pp = "particles/";
  if (pp.back() != '/') pp += '/';
  const std::string ip = itRoot + std::to_string(iteration) + "/" + pp;
  const std::vector<std::string> sps = F.children(ip);
  std::string sp = species;
  if (sp.empty()) {
    if (sps.size() != 1) F.fail("give file_species: the file holds the species " + join(sps));
    sp = sps[0];
  } else if (std::find(sps.begin(), sps.end(), sp) == sps.end()) {
    F.fail("species '" + sp + "' not found (species: " + join(sps) + ")");
  }
  B.species = sp;
  const std::string g = ip + sp + "/";
  const size_t n = F.count(g + "position/x");
  // charge and mass of one particle: electrons or positrons only
  const double q = F.exists(g + "charge") ? F.component(g + "charge", n)[0] : -QE;
  const double m = F.exists(g + "mass") ? F.component(g + "mass", n)[0] : ME;
  if (std::fabs(std::fabs(q) / QE - 1.0) > 1e-3 || std::fabs(m / ME - 1.0) > 1e-3)
    F.fail("species '" + sp + "' is not electrons or positrons (charge " + std::to_string(q) + " C, mass " +
           std::to_string(m) + " kg)");
  B.charge_sign = q < 0 ? -1.0 : 1.0;
  const std::vector<double> w = F.exists(g + "weighting") ? F.component(g + "weighting", n)
                                                          : std::vector<double>(n, 1.0);
  std::vector<double> pos[3], mom[3];
  const char* ax[3] = {"x", "y", "z"};
  for (int a = 0; a < 3; ++a) {
    pos[a] = F.component(g + "position/" + ax[a], n, a != 1 || F.exists(g + "position/y"));
    const std::vector<double> off = F.component(g + "positionOffset/" + ax[a], n, false);
    for (size_t k = 0; k < n; ++k) pos[a][k] += off[k];
    mom[a] = F.component(g + "momentum/" + ax[a], n, a != 1 || F.exists(g + "momentum/y"));
  }
  // momentum per real particle: macroWeighted = 1 records hold the sum over the macroparticle
  const bool mw = F.hasAttr(g + "momentum", "macroWeighted") && F.attrD(g + "momentum", "macroWeighted")[0] != 0.0;
  const double wp = F.attr1(g + "momentum", "weightingPower", 1.0);
  B.p.resize(n);
  for (size_t k = 0; k < n; ++k) {
    const double f = mw ? std::pow(w[k], -wp) : 1.0;
    B.p[k] = ImportedParticle{pos[0][k], pos[1][k], pos[2][k], mom[0][k] * f / (ME * CL), mom[1][k] * f / (ME * CL),
                              mom[2][k] * f / (ME * CL), w[k]};
  }
  return B;
}

ImportedBeam readNative(const File& F, const std::string& species) {
  ImportedBeam B;
  B.format = "GEM-PIC-RZ native";
  if (!F.exists("/particles/beam/x")) F.fail("no beam list in the file (written with output.particles = false?)");
  const double L = F.attrD("/", "k0inv_m")[0], Q = F.attrD("/", "charge_unit_C")[0];
  std::vector<std::string> names;
  {
    std::stringstream ss(F.attrS("/particles/beam", "bunch_names"));
    std::string t;
    while (std::getline(ss, t, ',')) if (!t.empty()) names.push_back(t);
  }
  std::map<std::string, int> ids;
  for (size_t b = 0; b < names.size(); ++b) ids[names[b]] = int(b);
  ids["trapped"] = -1;
  {
    std::stringstream ss(F.attrS("/particles/beam", "ion_names"));
    std::string t;
    int s = 0;
    while (std::getline(ss, t, ',')) if (!t.empty()) ids["inj_" + t] = -30 - s++;
  }
  std::vector<std::string> all;
  for (const auto& kv : ids) all.push_back(kv.first);
  std::string sp = species;
  if (sp.empty()) {
    if (names.size() != 1) F.fail("give file_species: the file holds " + join(all));
    sp = names[0];
  }
  if (!ids.count(sp)) F.fail("species '" + sp + "' not found (species: " + join(all) + ")");
  B.species = sp;
  const int id = ids[sp];
  const auto x = F.data("/particles/beam/x"), y = F.data("/particles/beam/y"), zeta = F.data("/particles/beam/zeta");
  const auto ux = F.data("/particles/beam/ux"), uy = F.data("/particles/beam/uy"), uz = F.data("/particles/beam/uz");
  const auto q = F.data("/particles/beam/charge"), s = F.data("/particles/beam/species");
  const double z0 = F.attrD("/", "z")[0];
  bool signSet = false;
  for (size_t k = 0; k < x.size(); ++k) {
    if (int(std::lround(s[k])) != id) continue;
    if (!signSet) { B.charge_sign = q[k] < 0 ? -1.0 : 1.0; signSet = true; }
    B.p.push_back(ImportedParticle{x[k] * L, y[k] * L, (zeta[k] + z0) * L, ux[k], uy[k], uz[k],
                                   std::fabs(q[k]) * Q / QE});
  }
  return B;
}

}  // namespace

ImportedBeam importBeam(const std::string& path, const std::string& species, long iteration) {
  H5Eset_auto2(H5E_DEFAULT, nullptr, nullptr);
  File F(path);
  ImportedBeam B;
  if (F.hasAttr("/", "openPMD")) B = readOpenPMD(F, species, iteration);
  else if (F.attrS("/", "code") == "GEM-PIC-RZ") B = readNative(F, species);
  else F.fail("neither an openPMD file nor a GEM-PIC-RZ data file");
  if (B.p.empty()) F.fail("species '" + B.species + "' holds no particles");
  return B;
}

}  // namespace gem

#else

namespace gem {
ImportedBeam importBeam(const std::string& path, const std::string&, long) {
  throw std::runtime_error(path + ": reading beams from files needs HDF5 support (install libhdf5-dev and rebuild)");
}
}  // namespace gem

#endif
