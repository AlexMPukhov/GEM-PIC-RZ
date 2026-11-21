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
#include <map>
#include <stdexcept>

#if defined(GEMPIC_HAVE_HDF5)
#include "h5read.hpp"
#endif

namespace gem {

namespace {

[[noreturn]] void failAt(const std::string& path, const std::string& m) { throw std::runtime_error(path + ": " + m); }

const char* HEADER_HELP =
    "the header (comment lines '# key = value') must state the geometry and the data layout, e.g.\n"
    "    # geometry = z                        or    # geometry = rz\n"
    "    # columns  = z[mm] n[cm^-3]                 # columns  = z[mm] r[um] n[cm^-3]\n"
    "(see docs/density_file.md; a density without unit, n or n[1], is a factor on 'density')";

struct Column {
  char kind = '?';  // 'z', 'r', 'x', 'y', 'n'
  std::string name, unit;
};

bool isSep(char c) { return c == ' ' || c == '\t' || c == ',' || c == ';' || c == '\r'; }

std::string trim(std::string s) {
  while (!s.empty() && (s.back() == ' ' || s.back() == '\t' || s.back() == '\r')) s.pop_back();
  size_t a = 0;
  while (a < s.size() && (s[a] == ' ' || s[a] == '\t')) ++a;
  return s.substr(a);
}

// "z[mm]  r [um], n(cm^-3)" -> columns; false if a name is not a coordinate (z, r, x, y) or a density
bool parseColumns(const std::string& s, std::vector<Column>& cols) {
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
      c.unit = trim(s.substr(q + 1, e - q - 1));
      p = e + 1;
    }
    std::string low;
    for (char ch : c.name) low += char(std::tolower(static_cast<unsigned char>(ch)));
    if (low == "z" || low == "r" || low == "x" || low == "y") c.kind = low[0];
    else if (low == "density" || (!low.empty() && low[0] == 'n')) c.kind = 'n';
    else return false;
    cols.push_back(c);
  }
  return !cols.empty();
}

bool dimensionless(const std::string& unit) { return unit.empty() || unit == "1" || unit == "-"; }

// peak and its position
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

// the transverse function of a file with geometry z: variables r, x, y; numbers with units allowed
void setTransverse(DensityFile& D, std::string text, const Units& u) {
  text = trim(text);
  if (text.size() >= 2 && text.front() == '"' && text.back() == '"') text = text.substr(1, text.size() - 2);
  if (D.twoD()) failAt(D.path, "a transverse function is only for geometry z (the file holds n(r, z))");
  try {
    D.transverse = Expr(text, {"r", "x", "y"}, std::map<std::string, double>(), u);
    const double t = D.transverse({0.0, 0.0, 0.0});
    if (!std::isfinite(t)) failAt(D.path, "transverse function '" + text + "' is not finite on the axis");
  } catch (const std::invalid_argument& e) {
    failAt(D.path, "transverse function '" + text + "': " + e.what());
  }
}

// 3D and Cartesian files of the common format: not for a quasi-cylindrical code
void needAxial(const std::string& path, const std::string& geo) {
  if (geo == "z" || geo == "rz") return;
  if (geo == "xyz")
    failAt(path, "geometry xyz: a 3D profile n(x, y, z) is for the 3D codes; GEM-PIC-RZ reads geometry z and rz "
                 "(the file is not averaged over the azimuth: give n(r, z))");
  if (geo == "xz")
    failAt(path, "geometry xz: a Cartesian 2D profile n(x, z); GEM-PIC-RZ reads geometry z and rz");
  failAt(path, "unknown geometry '" + geo + "' (z | rz | xz | xyz; GEM-PIC-RZ reads z and rz)");
}

std::shared_ptr<const DensityFile> readText(const std::string& path, const Units& u) {
  std::ifstream in(path);
  if (!in) failAt(path, "cannot open the file");
  auto D = std::make_shared<DensityFile>();
  D->path = path;
  D->format = "text";
  // ---- header: comment lines "# key = value" before the data
  std::string geometry, columns, densityUnit, order, transverse;
  bool haveDensityUnit = false;
  struct AxisLine { char name; std::string text; int line; };
  std::vector<AxisLine> axes;
  std::vector<std::string> dataLines;
  std::vector<int> dataNo;
  std::string line;
  int ln = 0;
  while (std::getline(in, line)) {
    ++ln;
    const std::string t = trim(line);
    if (t.empty()) continue;
    if (t[0] != '#') {
      dataLines.push_back(t);
      dataNo.push_back(ln);
      continue;
    }
    if (!dataLines.empty()) continue;  // comments between the data rows
    if (t.size() > 1 && t[1] == '#') continue;  // "## ...": free comment
    const size_t eq = t.find('=');
    if (eq == std::string::npos) continue;
    const std::string key = trim(t.substr(1, eq - 1)), val = trim(t.substr(eq + 1));
    bool ident = !key.empty() && (std::isalpha(static_cast<unsigned char>(key[0])) || key[0] == '_');
    for (char c : key) ident = ident && (std::isalnum(static_cast<unsigned char>(c)) || c == '_');
    if (!ident) continue;  // free text with an '=' in it
    const std::string where = "line " + std::to_string(ln) + ": ";
    if (key == "density_profile") {
      if (val != "1") failAt(path, where + "format version '" + val + "' (this reader knows version 1)");
    } else if (key == "geometry") geometry = val;
    else if (key == "columns") columns = val;
    else if (key == "density_unit") { densityUnit = val; haveDensityUnit = true; }
    else if (key == "order") order = val;
    else if (key == "transverse") transverse = val;
    else if (key == "z" || key == "r" || key == "x" || key == "y") axes.push_back({key[0], val, ln});
    else
      failAt(path, where + "unknown header key '" + key + "' (keys: density_profile, geometry, columns, "
                       "density_unit, z, r, x, y, order, transverse; write free comments without '=' or start "
                       "them with '##')");
  }
  if (geometry.empty()) failAt(path, std::string("no line '# geometry = ...' in the header; ") + HEADER_HELP);
  needAxial(path, geometry);
  const bool two = geometry == "rz";
  if (dataLines.empty()) failAt(path, "no data");
  if (columns.empty() == axes.empty())
    failAt(path, columns.empty() ? std::string("the header gives neither '# columns = ...' nor the axes of a grid; ") +
                                       HEADER_HELP
                                 : "the header gives both '# columns' and axis lines ('# z = ...'): choose one layout");
  D->detail = "geometry " + geometry + (columns.empty() ? ", grid layout" : ", columns");
  // numbers of a data line
  auto numbers = [&](size_t k, std::vector<double>& out) {
    const char* p = dataLines[k].c_str();
    for (;;) {
      while (*p && isSep(*p)) ++p;
      if (!*p || *p == '#') break;
      char* end = nullptr;
      const double val = std::strtod(p, &end);
      if (end == p)
        failAt(path, "line " + std::to_string(dataNo[k]) + ": not a number: '" + std::string(p).substr(0, 30) + "'");
      if (!std::isfinite(val)) failAt(path, "line " + std::to_string(dataNo[k]) + ": a value is not finite (nan or inf)");
      out.push_back(val);
      p = end;
    }
  };
  try {
    if (!columns.empty()) {
      // ---- columns layout: coordinates and density per row; the grid is taken from the rows
      if (haveDensityUnit || !order.empty())
        failAt(path, "density_unit / order belong to the grid layout; with '# columns' the units stand in that line");
      std::vector<Column> cols;
      if (!parseColumns(columns, cols)) failAt(path, "cannot read '# columns = " + columns + "'; " + HEADER_HELP);
      int iz = -1, ir = -1, in_ = -1;
      for (size_t k = 0; k < cols.size(); ++k) {
        int& slot = cols[k].kind == 'z' ? iz : cols[k].kind == 'r' ? ir : in_;
        if (cols[k].kind == 'x' || cols[k].kind == 'y' || slot >= 0)
          failAt(path, "'# columns = " + columns + "' does not fit geometry " + geometry + " (" +
                           (two ? "z, r and the density" : "z and the density") + ")");
        slot = int(k);
      }
      if (iz < 0 || in_ < 0 || (ir >= 0) != two)
        failAt(path, "'# columns = " + columns + "' does not fit geometry " + geometry + " (" +
                         (two ? "z, r and the density" : "z and the density") + ")");
      D->zUnit = cols[size_t(iz)].unit;
      D->nUnit = cols[size_t(in_)].unit;
      if (two) D->rUnit = cols[size_t(ir)].unit;
      D->absolute = !dimensionless(D->nUnit);
      struct Row { double z, r, n; int line; };
      std::vector<Row> rows;
      std::vector<double> x;
      for (size_t k = 0; k < dataLines.size(); ++k) {
        x.clear();
        numbers(k, x);
        if (x.size() != cols.size())
          failAt(path, "line " + std::to_string(dataNo[k]) + ": " + std::to_string(x.size()) + " numbers, '# columns' names " +
                           std::to_string(cols.size()));
        if (x[size_t(in_)] < 0.0)
          failAt(path, "line " + std::to_string(dataNo[k]) + ": negative density " + num(x[size_t(in_)]));
        // each value converted as a deck value in this unit would be
        rows.push_back({toCode(x[size_t(iz)], D->zUnit, Dim::Length, u),
                        two ? toCode(x[size_t(ir)], D->rUnit, Dim::Length, u) : 0.0,
                        D->absolute ? toCode(x[size_t(in_)], D->nUnit, Dim::Density, u) : x[size_t(in_)], dataNo[k]});
      }
      if (!two) {
        if (rows.size() < 2) failAt(path, "a profile n(z) needs at least 2 points");
        for (size_t k = 0; k < rows.size(); ++k) {
          if (k > 0 && !(rows[k].z > rows[k - 1].z))
            failAt(path, "line " + std::to_string(rows[k].line) + ": z must be strictly increasing");
          D->z.push_back(rows[k].z);
          D->v.push_back(rows[k].n);
        }
      } else {  // the rows fill the grid of the distinct z and r values
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
          failAt(path, "a profile n(r, z) needs at least 2 values of z and of r (found " + std::to_string(nz) +
                           " and " + std::to_string(nr) + ")");
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
    } else {
      // ---- grid layout: the axes in the header, the values alone as data
      if (!haveDensityUnit)
        failAt(path, "the grid layout needs '# density_unit = ...' (cm^-3, m^-3, n_c; 1: a factor on 'density')");
      D->nUnit = densityUnit;
      D->absolute = !dimensionless(densityUnit);
      std::string axisOrder;  // names, outer to inner: the last one varies fastest
      for (const auto& a : axes) {
        const std::string where = "line " + std::to_string(a.line) + ": ";
        if ((a.name != 'z' && a.name != 'r') || (a.name == 'r' && !two))
          failAt(path, where + "axis '" + std::string(1, a.name) + "' does not fit geometry " + geometry);
        if (axisOrder.find(a.name) != std::string::npos) failAt(path, where + "axis given twice");
        axisOrder += a.name;
        // "uniform <origin> <step> <count> <unit>"  or  "list <x1> <x2> ... <unit>"
        std::vector<std::string> tok;
        std::string cur;
        for (char c : a.text + " ") {
          if (isSep(c)) { if (!cur.empty()) tok.push_back(cur); cur.clear(); }
          else cur += c;
        }
        const char* form = "'uniform <origin> <step> <count> <unit>' or 'list <x1> <x2> ... <unit>'";
        if (tok.size() < 4) failAt(path, where + "an axis is " + form);
        const std::string unit = tok.back();
        std::vector<double> xs;
        for (size_t k = 1; k + 1 < tok.size(); ++k) {
          char* end = nullptr;
          xs.push_back(std::strtod(tok[k].c_str(), &end));
          if (end == tok[k].c_str() || *end || !std::isfinite(xs.back()))
            failAt(path, where + "not a number: '" + tok[k] + "' (an axis is " + form + ")");
        }
        std::vector<double>& ax = a.name == 'z' ? D->z : D->r;
        (a.name == 'z' ? D->zUnit : D->rUnit) = unit;
        if (tok[0] == "uniform") {
          if (xs.size() != 3 || xs[2] < 2 || xs[2] != std::floor(xs[2]) || !(xs[1] > 0))
            failAt(path, where + "uniform axis: origin, step > 0 and a count >= 2, then the unit");
          for (long k = 0; k < long(xs[2]); ++k) ax.push_back(toCode(xs[0] + double(k) * xs[1], unit, Dim::Length, u));
        } else if (tok[0] == "list") {
          if (xs.size() < 2) failAt(path, where + "a list axis needs at least 2 coordinates");
          for (double x : xs) ax.push_back(toCode(x, unit, Dim::Length, u));
        } else {
          failAt(path, where + "an axis is " + form);
        }
        for (size_t k = 1; k < ax.size(); ++k)
          if (!(ax[k] > ax[k - 1])) failAt(path, where + "the coordinates must be strictly increasing");
        if (a.name == 'r' && ax.front() < 0.0) failAt(path, where + "negative radius");
      }
      if (axisOrder.size() != (two ? 2u : 1u))
        failAt(path, "geometry " + geometry + " needs the ax" + (two ? "es '# z = ...' and '# r = ...'" : "is '# z = ...'"));
      if (!order.empty()) {
        std::string o;
        for (char c : order)
          if (!isSep(c)) o += c;
        if (o != "zr" && o != "rz" && o != "z") failAt(path, "'# order = " + order + "': the axis names, e.g. 'z r'");
        if (o.size() != axisOrder.size()) failAt(path, "'# order = " + order + "' does not fit geometry " + geometry);
        axisOrder = o;
      }
      std::vector<double> x;
      for (size_t k = 0; k < dataLines.size(); ++k) numbers(k, x);
      const size_t nz = D->z.size(), nr = D->nr();
      if (x.size() != nz * nr)
        failAt(path, std::to_string(x.size()) + " values, but the axes of the header give " +
                         (two ? std::to_string(nz) + " x " + std::to_string(nr) + " = " : std::string()) +
                         std::to_string(nz * nr) + " grid points");
      D->v.resize(nz * nr);
      for (size_t i = 0; i < nz; ++i)
        for (size_t j = 0; j < nr; ++j) {
          const double val = x[axisOrder == "rz" ? j * nz + i : i * nr + j];
          if (val < 0.0) failAt(path, "negative density " + num(val));
          D->v[i * nr + j] = D->absolute ? toCode(val, densityUnit, Dim::Density, u) : val;
        }
    }
  } catch (const std::invalid_argument& e) {
    failAt(path, std::string("units of the header: ") + e.what());
  }
  if (!transverse.empty()) setTransverse(*D, transverse, u);
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
  {  // the geometry of the common format (docs/density_file.md) from the axes: z | rz | xz | xyz
    std::string g;
    for (const char* a : {"r", "x", "y", "z"})
      if (std::find(lab.begin(), lab.end(), a) != lab.end()) g += a;
    if (g == "xz" || g == "xyz") needAxial(path, g);
  }
  const std::vector<hsize_t> dims = F.dims(p);
  const size_t lead =(theta && dims.size() == na + 1) ? 1 : 0;  // thetaMode: the first index is the mode
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
  D->detail += std::string(", geometry ") + (ar >= 0 ? "rz" : "z");
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
  // transverse function of a mesh n(z): string attribute transverseFunction (variables r, x, y)
  if (F.hasAttr(p, "transverseFunction")) setTransverse(*D, F.attrS(p, "transverseFunction"), u);
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
