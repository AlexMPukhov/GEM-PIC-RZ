// ============================================================================
//  density_file.hpp -- plasma / gas density profiles from external files
//  (plasma.density_file, ion.density_file with profile = own)
//
//  The format is specified in docs/density_file.md (one format for several
//  codes, geometry z | rz | xz | xyz; this code reads z and rz). In short:
//
//  1. Text table. The header, comment lines "# key = value", defines the
//     geometry and the data layout:
//         # geometry = z                 # geometry = rz
//         # columns  = z[mm] n[cm^-3]    # columns  = z[mm] r[um] n[cm^-3]
//     columns layout: coordinates and density per row (z: increasing; rz: one
//     row per grid point in any order, the rows fill the grid z x r);
//     or the grid layout, the axes in the header and the values alone as data:
//         # geometry = rz
//         # density_unit = cm^-3
//         # z = uniform 0 0.05 41 mm        (origin, step, count, unit)
//         # r = list 0 5 10 20 40 um
//         # order = z r                     (the last axis varies fastest)
//     A file with geometry z may carry a transverse function of r, x, y:
//         # transverse = "exp(-(r / (50 um))^2)"
//     Density: a physical unit (cm^-3, m^-3, n_c) = the absolute density;
//     no unit, 1 or - = a factor on the `density` of the block.
//  2. openPMD (HDF5): a scalar mesh record (default name "density") with the
//     axes z, or r and z (also thetaMode: the mode m = 0 is taken). Regular
//     axes (gridSpacing, gridGlobalOffset, position) or explicit coordinates
//     (attribute gridCoordinates, docs/openpmd_rectilinear.md). unitDimension
//     L^-3: absolute density (unitSI -> m^-3); all zero: a factor. Optional
//     string attribute transverseFunction of a mesh n(z).
//
//  Between the points the density is linear in z and r (bilinear). Inside the
//  first radius of the file it is the value there (the axis side); outside the
//  file it is zero, or the edge value with `hold`.
// ============================================================================
#pragma once
#include <cmath>
#include <memory>
#include <string>
#include <vector>

#include "expr.hpp"
#include "units.hpp"

namespace gem {

struct DensityFile {
  std::string path;
  std::string format;         // "text" or "openPMD"
  std::string detail;         // geometry and layout (text), mesh and iteration (openPMD)
  std::string zUnit, rUnit, nUnit;  // units as stated in the file
  bool absolute = false;      // values are densities [n_c]; else a factor on the block's density
  std::vector<double> z, r;   // strictly increasing coordinates [1/k0]; r empty: n(z) only
  std::vector<double> v;      // values [iz * nr + ir], nr = max(1, r.size())
  double peak = 0.0, zPeak = 0.0, rPeak = 0.0;
  Expr transverse;            // geometry z: transverse factor T(r, x, y) (empty: 1)

  bool twoD() const { return !r.empty(); }
  size_t nr() const { return r.empty() ? 1 : r.size(); }
  // value at (r, z) of the file's coordinates, without the transverse function; outside the file 0, or the
  // edge value (hold)
  double at(double rr, double zz, bool hold) const;
  // transverse factor at (x, y); values < 0 count as 0
  double trans(double x, double y) const {
    if (transverse.empty()) return 1.0;
    const double t = transverse({std::sqrt(x * x + y * y), x, y});
    return t > 0.0 ? t : 0.0;
  }
  // its mean over the azimuth at the radius rr (16 angles): the axisymmetric part; exact for a function of r
  double transMean(double rr) const {
    if (transverse.empty()) return 1.0;
    double s = 0.0;
    for (int k = 0; k < 16; ++k) {
      const double ph = 0.39269908169872414 * k;  // 2 pi / 16
      const double t = transverse({rr, rr * std::cos(ph), rr * std::sin(ph)});
      s += t > 0.0 ? t : 0.0;
    }
    return s / 16.0;
  }
};

// Reads a text table or an openPMD file (recognized by the HDF5 signature). `mesh` and `iteration`
// (< 0: the last) select the openPMD record. Throws std::runtime_error with the file name in the message.
std::shared_ptr<const DensityFile> readDensityFile(const std::string& path, const std::string& mesh, long iteration,
                                                   const Units& u);

}  // namespace gem
