// ============================================================================
//  beam_import.cpp -- read a particle beam from openPMD or native HDF5 files
// ============================================================================
#include "beam_import.hpp"

#include <cmath>
#include <stdexcept>

#if defined(GEMPIC_HAVE_HDF5)
#include "h5read.hpp"

#include <algorithm>
#include <map>
#include <sstream>

namespace gem {

namespace {

constexpr double QE = 1.602176634e-19, ME = 9.1093837015e-31, CL = 299792458.0;

using h5r::H;

// openPMD record helpers on top of the shared reader
class File : public h5r::File {
 public:
  using h5r::File::File;
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
  const std::vector<double> w = F.exists(g + "weighting") ? F.component(g + "weighting", n)
                                                          : std::vector<double>(n, 1.0);
  // charge and mass of one real particle (macroWeighted records: divided by w^weightingPower); one value for the
  // whole species, the charge a non-zero multiple of e (no electrons when absent: the default species)
  auto perParticle = [&](const std::string& rec, double dflt) {
    if (!F.exists(g + rec)) return std::vector<double>(n, dflt);
    std::vector<double> v = F.component(g + rec, n);
    if (F.hasAttr(g + rec, "macroWeighted") && F.attrD(g + rec, "macroWeighted")[0] != 0.0) {
      const double pw = F.attr1(g + rec, "weightingPower", 1.0);
      for (size_t k = 0; k < n; ++k) v[k] *= std::pow(w[k], -pw);
    }
    for (size_t k = 1; k < n; ++k)
      if (std::fabs(v[k] - v[0]) > 1e-6 * std::fabs(v[0]))
        F.fail("species '" + sp + "': the " + rec + " differs between particles (one species with one " + rec +
               " is needed; select the particles of one charge state in a separate species)");
    return v;
  };
  const double q = perParticle("charge", -QE)[0], m = perParticle("mass", ME)[0];
  const double Z = q / QE;
  if (!(std::fabs(Z) >= 0.999) || std::fabs(Z - std::round(Z)) > 1e-3)
    F.fail("species '" + sp + "': the charge " + std::to_string(q) + " C is not a non-zero multiple of e");
  if (!(m > 0.0)) F.fail("species '" + sp + "': the mass " + std::to_string(m) + " kg is not positive");
  B.charge_e = std::round(Z);
  B.mass = std::fabs(m / ME - 1.0) < 1e-6 ? 1.0 : m / ME;  // electrons and positrons: exactly m_e
  const double MC = B.mass * ME * CL;                     // momenta in units of m c of the species
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
    B.p[k] = ImportedParticle{pos[0][k], pos[1][k], pos[2][k], mom[0][k] * f / MC, mom[1][k] * f / MC,
                              mom[2][k] * f / MC, w[k]};
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
  // mass [m_e] and charge [e] of one particle of a bunch (attributes bunch_masses, bunch_charges; files written
  // before they existed hold electron or positron bunches only); trapped species are electrons
  auto listAttr = [&](const char* a, size_t b, double dflt) {
    if (id < 0 || !F.hasAttr("/particles/beam", a)) return dflt;
    std::stringstream ss(F.attrS("/particles/beam", a));
    std::string t;
    for (size_t i = 0; std::getline(ss, t, ','); ++i)
      if (i == b) return std::stod(t);
    F.fail(std::string("attribute ") + a + " has no entry for bunch '" + sp + "'");
    return dflt;
  };
  B.mass = listAttr("bunch_masses", size_t(std::max(id, 0)), 1.0);
  const double Zb = listAttr("bunch_charges", size_t(std::max(id, 0)), 0.0);  // 0: sign from the particle charge
  const auto x = F.data("/particles/beam/x"), y = F.data("/particles/beam/y"), zeta = F.data("/particles/beam/zeta");
  const auto ux = F.data("/particles/beam/ux"), uy = F.data("/particles/beam/uy"), uz = F.data("/particles/beam/uz");
  const auto q = F.data("/particles/beam/charge"), s = F.data("/particles/beam/species");
  const double z0 = F.attrD("/", "z")[0];
  bool signSet = false;
  for (size_t k = 0; k < x.size(); ++k) {
    if (int(std::lround(s[k])) != id) continue;
    if (!signSet) { B.charge_e = Zb != 0.0 ? Zb : (q[k] < 0 ? -1.0 : 1.0); signSet = true; }
    B.p.push_back(ImportedParticle{x[k] * L, y[k] * L, (zeta[k] + z0) * L, ux[k], uy[k], uz[k],
                                   std::fabs(q[k]) * Q / (QE * std::fabs(B.charge_e))});
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
