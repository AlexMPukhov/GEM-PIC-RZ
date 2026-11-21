// ============================================================================
//  density_file.cpp -- density profiles from text tables and openPMD meshes
//  (see density_file.hpp)
// ============================================================================
#include "density_file.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <stdexcept>

#if defined(GEMPIC_HAVE_HDF5)
#include "h5read.hpp"
#endif

namespace gem {

namespace {

[[noreturn]] void failAt(const std::string& path, const std::string& m) { throw std::runtime_error(path + ": " + m); }

const char* HEADER_HELP =
    "a header line must name the columns with their units, e.g. '# z[mm] n[cm^-3]' for n(z) or "
    "'# z[mm] r[um] n[cm^-3]' for n(r, z); a density without unit, n or n[1], is a factor on 'density'";

struct Column {
  char kind = '?';  // 'z', 'r', 'n'
  std::string name, unit;
};

bool isSep(char c) { return c == ' ' || c == '\t' || c == ',' || c == ';' || c == '\r'; }

// "z[mm]  r [um], n(cm^-3)" -> columns; false if a name is not z, r or a density
bool parseHeader(const std::string& s, std::vector<Column>& cols) {
  cols.clear();
  size_t p = 0;
  while (p < s.size()) {
    while (p < s.size() && isSep(s[p])) ++p;
    if (p >= s.size()) break;
    Column c;
    while (p < s.size() && !isSep(s[p]) && s[p] != '[' && s[p] != '(') c.name += s[p++];
    size_t q = p;
    while (q < s.size() && (s[q] == ' ' || s[q] == '\t')) ++q;
    if (q < s.size() && (s[q] == '[' || s[q] == '(')) {
      const char close = s[q] == '[' ? ']' : ')';
      const size_t e = s.find(close, q);
      if (e == std::string::npos) return false;
      c.unit = s.substr(q + 1, e - q - 1);
      while (!c.unit.empty() && c.unit.back() == ' ') c.unit.pop_back();
      while (!c.unit.empty() && c.unit.front() == ' ') c.unit.erase(0, 1);
      p = e + 1;
    }
    std::string low;
    for (char ch : c.name) low += char(std::tolower(static_cast<unsigned char>(ch)));
    if (low == "z") c.kind = 'z';
    else if (low == "r") c.kind = 'r';
    else if (low == "density" || (!low.empty() && low[0] == 'n')) c.kind = 'n';
    else return false;
    cols.push_back(c);
  }
  int nz = 0, nr = 0, nn = 0;
  for (const auto& c : cols) { nz += c.kind == 'z'; nr += c.kind == 'r'; nn += c.kind == 'n'; }
  return nz == 1 && nn == 1 && nr <= 1 && cols.size() == size_t(nz + nr + nn);
}

bool dimensionless(const std::string& unit) { return unit.empty() || unit == "1" || unit == "-"; }

// peak and its position; rejects negative and non-finite values
void finish(DensityFile& D) {
  const size_t nr = D.nr();
  for (size_t i = 0; i < D.z.size(); ++i)
    for (size_t j = 0; j < nr; ++j) {
      const double x = D.v[i * nr + j];
      if (x > D.peak) { D.peak = x; D.zPeak = D.z[i]; D.rPeak = D.r.empty() ? 0.0 : D.r[j]; }
    }
  if (!(D.peak > 0.0)) failAt(D.path, "the density is zero everywhere");
}

std::string num(double x) {
  char b[64];
  std::snprintf(b, sizeof b, "%.10g", x);
  return b;
}

std::shared_ptr<const DensityFile> readText(const std::string& path, const Units& u) {
  std::ifstream in(path);
  if (!in) failAt(path, "cannot open the file");
  auto D = std::make_shared<DensityFile>();
  D->path = path;
  D->format = "text";
  std::vector<Column> cols, tmp;
  bool haveHeader = false, data = false;
  struct Row { double z, r, n; int line; };
  std::vector<Row> rows;
  std::string line;
  int ln = 0;
  int iz = -1, ir = -1, in_ = -1;
  while (std::getline(in, line)) {
    ++ln;
    size_t a = 0;
    while (a < line.size() && (line[a] == ' ' || line[a] == '\t' || line[a] == '\r')) ++a;
    if (a >= line.size()) continue;
    const auto where = [&]() { return "line " + std::to_string(ln) + ": "; };
    if (line[a] == '#') {
      const size_t b = line.find_first_not_of("# \t\r", a);
      if (!data && b != std::string::npos && parseHeader(line.substr(b), tmp)) {
        cols = tmp;  // the last comment line that names the columns
        haveHeader = true;
        D->detail = line.substr(a);
      }
      continue;
    }
    // a first line that does not start with a number: the header of a CSV file
    char* end = nullptr;
    std::strtod(line.c_str() + a, &end);
    if (end == line.c_str() + a) {
      if (data) failAt(path, where() + "not a number: '" + line.substr(a, 30) + "'");
      if (!parseHeader(line.substr(a), tmp)) failAt(path, where() + "cannot read the header; " + HEADER_HELP);
      cols = tmp;
      haveHeader = true;
      D->detail = line.substr(a);
      continue;
    }
    if (!haveHeader) failAt(path, std::string("no header before the data; ") + HEADER_HELP);
    if (!data) {
      for (size_t k = 0; k < cols.size(); ++k) {
        if (cols[k].kind == 'z') iz = int(k);
        if (cols[k].kind == 'r') ir = int(k);
        if (cols[k].kind == 'n') in_ = int(k);
      }
      data = true;
    }
    // numbers of the row
    double x[3] = {0, 0, 0};
    size_t n = 0;
    const char* p = line.c_str() + a;
    for (;;) {
      while (*p && isSep(*p)) ++p;
      if (!*p || *p == '#') break;
      const double val = std::strtod(p, &end);
      if (end == p) failAt(path, where() + "not a number: '" + std::string(p).substr(0, 30) + "'");
      if (n >= cols.size())
        failAt(path, where() + "more than " + std::to_string(cols.size()) + " columns (the header names " +
                         std::to_string(cols.size()) + ")");
      x[n++] = val;
      p = end;
    }
    if (n != cols.size())
      failAt(path, where() + std::to_string(n) + " columns, the header names " + std::to_string(cols.size()));
    for (size_t k = 0; k < n; ++k)
      if (!std::isfinite(x[k])) failAt(path, where() + "a value is not finite (nan or inf)");
    if (x[in_] < 0.0) failAt(path, where() + "negative density " + num(x[in_]));
    rows.push_back({x[iz], ir >= 0 ? x[ir] : 0.0, x[in_], ln});
  }
  if (!haveHeader) failAt(path, std::string("no header; ") + HEADER_HELP);
  D->zUnit = cols[size_t(iz < 0 ? 0 : iz)].unit;
  if (rows.empty()) failAt(path, "no data rows");
  // units (each value converted as a deck value in this unit would be)
  D->nUnit = cols[size_t(in_)].unit;
  D->absolute = !dimensionless(D->nUnit);
  if (ir >= 0) D->rUnit = cols[size_t(ir)].unit;
  try {
    for (auto& w : rows) {
      w.z = toCode(w.z, D->zUnit, Dim::Length, u);
      if (ir >= 0) w.r = toCode(w.r, D->rUnit, Dim::Length, u);
      if (D->absolute) w.n = toCode(w.n, D->nUnit, Dim::Density, u);
    }
  } catch (const std::invalid_argument& e) {
    failAt(path, std::string("header '") + D->detail + "': " + e.what());
  }
  if (ir < 0) {  // n(z)
    if (rows.size() < 2) failAt(path, "a profile n(z) needs at least 2 points");
    for (size_t k = 0; k < rows.size(); ++k) {
      if (k > 0 && !(rows[k].z > rows[k - 1].z))
        failAt(path, "line " + std::to_string(rows[k].line) + ": z must be strictly increasing");
      D->z.push_back(rows[k].z);
      D->v.push_back(rows[k].n);
    }
  } else {  // n(r, z): the rows fill the grid of the distinct z and r values
    for (const auto& w : rows) {
      if (w.r < 0.0) failAt(path, "line " + std::to_string(w.line) + ": negative radius");
      D->z.push_back(w.z);
      D->r.push_back(w.r);
    }
    std::sort(D->z.begin(), D->z.end());
    D->z.erase(std::unique(D->z.begin(), D->z.end()), D->z.end());
    std::sort(D->r.begin(), D->r.end());
    D->r.erase(std::unique(D->r.begin(), D->r.end()), D->r.end());
    const size_t nz = D->z.size(), nr = D->r.size();
    if (nz < 2 || nr < 2)
      failAt(path, "a profile n(r, z) needs at least 2 values of z and of r (found " + std::to_string(nz) + " and " +
                       std::to_string(nr) + "); for n(z) leave the r column out");
    if (rows.size() != nz * nr)
      failAt(path, "the rows do not fill a rectilinear grid: " + std::to_string(rows.size()) + " rows, but " +
                       std::to_string(nz) + " distinct z x " + std::to_string(nr) + " distinct r = " +
                       std::to_string(nz * nr) + " grid points (every (z, r) pair once)");
    D->v.assign(nz * nr, -1.0);
    for (const auto& w : rows) {
      const size_t i = size_t(std::lower_bound(D->z.begin(), D->z.end(), w.z) - D->z.begin());
      const size_t j = size_t(std::lower_bound(D->r.begin(), D->r.end(), w.r) - D->r.begin());
      if (D->v[i * nr + j] >= 0.0)
        failAt(path, "line " + std::to_string(w.line) + ": the grid point (z, r) is given twice");
      D->v[i * nr + j] = w.n;
    }
  }
  finish(*D);
  return D;
}

#if defined(GEMPIC_HAVE_HDF5)
std::shared_ptr<const DensityFile> readOpenPMD(const std::string& path, const std::string& mesh, long iteration,
                                               const Units& u) {
  H5Eset_auto2(H5E_DEFAULT, nullptr, nullptr);
  h5r::File F(path);
  auto D = std::make_shared<DensityFile>();
  D->path = path;
  D->format = "openPMD";
  if (!F.hasAttr("/", "openPMD")) F.fail("not an openPMD file (no attribute 'openPMD')");
  std::string itRoot;
  const long it = F.iteration(iteration, itRoot);
  std::string mp = F.attrS("/", "meshesPath");
  if (mp.empty()) mp = "meshes/";
  if (mp.back() != '/') mp += '/';
  const std::string mroot = itRoot + std::to_string(it) + "/" + mp;
  const std::vector<std::string> all = F.children(mroot);
  std::string name = mesh;
  if (name.empty()) {
    if (F.exists(mroot + "density")) name = "density";
    else if (all.size() == 1) name = all[0];
  }
  if (name.empty() || !F.exists(mroot + name)) {
    std::string l;
    for (const auto& c : all) l += (l.empty() ? "" : ", ") + c;
    F.fail((name.empty() ? std::string("no mesh 'density'") : "no mesh '" + name + "'") + " in iteration " +
           std::to_string(it) + " (meshes: " + (l.empty() ? "none" : l) + "); give file_mesh");
  }
  const std::string p = mroot + name;
  if (F.isGroup(p)) F.fail("mesh '" + name + "' is not a scalar record (a density is one dataset)");
  D->detail = "mesh '" + name + "', iteration " + std::to_string(it);
  const std::string geo = F.attrS(p, "geometry");
  const bool theta = geo == "thetaMode" || F.attrS(p, "geometryParameters").rfind("thetaMode", 0) == 0;
  std::vector<std::string> lab = F.attrSV(p, "axisLabels");
  const size_t na = lab.size();
  std::vector<double> sp(na, 1.0), off(na, 0.0), pos(na, 0.0);
  std::vector<std::string> gc(na);
  auto take = [&](const char* a, std::vector<double>& v, bool required) {
    if (!F.hasAttr(p, a)) {
      if (required) F.fail("mesh '" + name + "' has no attribute " + a);
      return;
    }
    v = F.attrD(p, a);
    if (v.size() != na) F.fail("mesh '" + name + "': " + a + " must have one entry per axis of axisLabels");
  };
  if (F.hasAttr(p, "gridCoordinates")) {
    gc = F.attrSV(p, "gridCoordinates");
    if (gc.size() != na) F.fail("mesh '" + name + "': gridCoordinates must have one entry per axis of axisLabels");
  }
  bool anyRegular = false;
  for (const auto& s : gc) anyRegular = anyRegular || s.empty();
  take("gridSpacing", sp, anyRegular);
  take("gridGlobalOffset", off, false);
  take("position", pos, false);
  if (F.attrS(p, "dataOrder") == "F") {
    std::reverse(lab.begin(), lab.end());
    std::reverse(sp.begin(), sp.end());
    std::reverse(off.begin(), off.end());
    std::reverse(pos.begin(), pos.end());
    std::reverse(gc.begin(), gc.end());
  }
  const std::vector<hsize_t> dims = F.dims(p);
  const size_t lead = (theta && dims.size() == na + 1) ? 1 : 0;  // thetaMode: the first index is the mode
  if (dims.size() != na + lead || na < 1 || na > 2) {
    std::string l;
    for (const auto& c : lab) l += (l.empty() ? "" : ", ") + c;
    F.fail("mesh '" + name + "': need a dataset over the axes z, or r and z (axisLabels: " + l + "; " +
           std::to_string(dims.size()) + " dimensions)");
  }
  int az = -1, ar = -1;
  for (size_t k = 0; k < na; ++k) {
    if (lab[k] == "z") az = int(k);
    else if (lab[k] == "r") ar = int(k);
  }
  if (az < 0 || (na == 2 && ar < 0)) {
    std::string l;
    for (const auto& c : lab) l += (l.empty() ? "" : ", ") + c;
    F.fail("mesh '" + name + "': the axes must be 'z', or 'r' and 'z' (axisLabels: " + l + ")");
  }
  const double gu = F.attr1(p, "gridUnitSI", 1.0), Lm = u.k0inv_um() * 1e-6;
  auto coords = [&](int k) {
    const size_t n = size_t(dims[lead + size_t(k)]);
    std::vector<double> x(n);
    if (!gc[size_t(k)].empty()) {  // explicit coordinates, relative to gridGlobalOffset (docs/openpmd_rectilinear.md)
      const std::vector<double> c = F.data(gc[size_t(k)]);
      if (c.size() != n)
        F.fail("coordinates " + gc[size_t(k)] + ": " + std::to_string(c.size()) + " values for " + std::to_string(n) +
               " data points along " + lab[size_t(k)]);
      for (size_t j = 0; j < n; ++j) x[j] = (off[size_t(k)] + c[j]) * gu / Lm;
    } else {
      for (size_t j = 0; j < n; ++j) x[j] = (off[size_t(k)] + (double(j) + pos[size_t(k)]) * sp[size_t(k)]) * gu / Lm;
    }
    for (size_t j = 0; j < n; ++j) {
      if (!std::isfinite(x[j])) F.fail("coordinates along " + lab[size_t(k)] + " are not finite");
      if (j > 0 && !(x[j] > x[j - 1])) F.fail("coordinates along " + lab[size_t(k)] + " must be strictly increasing");
    }
    return x;
  };
  D->z = coords(az);
  if (ar >= 0) D->r = coords(ar);
  D->zUnit = D->rUnit = "m";
  if (D->z.size() < 2 || (ar >= 0 && D->r.size() < 2)) F.fail("mesh '" + name + "': at least 2 points per axis");
  if (ar >= 0 && D->r.front() < 0.0) F.fail("mesh '" + name + "': negative radius");
  // unit: number density (L^-3) or dimensionless
  const double uSI = F.attr1(p, "unitSI", 1.0);
  std::vector<double> ud(7, 0.0);
  if (F.hasAttr(p, "unitDimension")) ud = F.attrD(p, "unitDimension");
  ud.resize(7, 0.0);
  bool zero = true, rest = true;
  for (size_t k = 0; k < 7; ++k) {
    zero = zero && ud[k] == 0.0;
    if (k > 0) rest = rest && ud[k] == 0.0;
  }
  if (zero) {
    D->absolute = false;
    D->nUnit = "1";
  } else if (ud[0] == -3.0 && rest) {
    D->absolute = true;
    D->nUnit = "m^-3";
  } else {
    F.fail("mesh '" + name + "': unitDimension is neither a number density (L^-3) nor dimensionless "
           "(a charge or mass density must be converted to a number density first)");
  }
  const double fac = D->absolute ? uSI * 1e-6 / u.nc_cm3() : uSI;
  const std::vector<double> raw = F.data(p);
  const size_t nz = D->z.size(), nr = D->nr();
  D->v.resize(nz * nr);
  for (size_t i = 0; i < nz; ++i)
    for (size_t j = 0; j < nr; ++j) {
      // mode 0 of a thetaMode record: the leading index 0
      const size_t idx = (ar < 0) ? i : (ar < az ? j * nz + i : i * nr + j);
      const double x = raw[idx] * fac;
      if (!std::isfinite(x)) F.fail("mesh '" + name + "': a value is not finite (nan or inf)");
      if (x < 0.0) F.fail("mesh '" + name + "': negative density " + num(raw[idx]));
      D->v[i * nr + j] = x;
    }
  if (lead && dims[0] > 1) D->detail += ", mode m = 0 of " + std::to_string((dims[0] + 1) / 2) + " modes";
  finish(*D);
  return D;
}
#endif

}  // namespace

double DensityFile::at(double rr, double zz, bool hold) const {
  if (zz < z.front() || zz > z.back()) {
    if (!hold) return 0.0;
    zz = zz < z.front() ? z.front() : z.back();
  }
  // the interval [z_i, z_i+1] and the weight, as for the tables of the deck
  size_t i = z.size() - 2;
  double tz = 1.0;
  if (zz < z.back()) {
    i = size_t(std::upper_bound(z.begin(), z.end(), zz) - z.begin()) - 1;
    tz = (zz - z[i]) / (z[i + 1] - z[i]);
  }
  if (r.empty()) {
    if (zz >= z.back()) return v.back();
    return v[i] + (v[i + 1] - v[i]) * tz;
  }
  const size_t n = r.size();
  if (rr > r.back()) {
    if (!hold) return 0.0;
    rr = r.back();
  }
  size_t j = n - 2;
  double tr = 1.0;
  if (rr <= r.front()) { j = 0; tr = 0.0; }  // inside the first radius: its value
  else if (rr < r.back()) {
    j = size_t(std::upper_bound(r.begin(), r.end(), rr) - r.begin()) - 1;
    tr = (rr - r[j]) / (r[j + 1] - r[j]);
  }
  const double a = v[i * n + j] + (v[i * n + j + 1] - v[i * n + j]) * tr;
  const double b = v[(i + 1) * n + j] + (v[(i + 1) * n + j + 1] - v[(i + 1) * n + j]) * tr;
  return a + (b - a) * tz;
}

std::shared_ptr<const DensityFile> readDensityFile(const std::string& path, const std::string& mesh, long iteration,
                                                   const Units& u) {
  // HDF5 signature \211HDF\r\n\032\n
  bool hdf5 = false;
  {
    std::ifstream in(path, std::ios::binary);
    if (!in) failAt(path, "cannot open the file");
    char sig[8] = {0};
    in.read(sig, 8);
    hdf5 = in.gcount() == 8 && std::string(sig, 8) == std::string("\211HDF\r\n\032\n", 8);
  }
  if (!hdf5) {
    if (!mesh.empty() || iteration >= 0)
      failAt(path, "file_mesh / file_iteration are for openPMD (HDF5) files; this is a text table");
    return readText(path, u);
  }
#if defined(GEMPIC_HAVE_HDF5)
  return readOpenPMD(path, mesh, iteration, u);
#else
  (void)mesh; (void)iteration;
  failAt(path, "reading openPMD density files needs HDF5 support (install libhdf5-dev and rebuild), "
               "or give the profile as a text table");
#endif
}

}  // namespace gem
