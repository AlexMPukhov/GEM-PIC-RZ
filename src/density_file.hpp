// ============================================================================
//  density_file.hpp -- plasma / gas density profiles from external files
//  (plasma.density_file, ion.density_file with profile = own)
//
//  A file holds the density on the axis, n(z), or on a rectilinear grid of the
//  r-z plane, n(r, z): measured gas-jet profiles (interferometry, Abel
//  inversion), fluid simulations of a nozzle, the hand-off of a hydro code.
//
//  1. Text / CSV table. Lines starting with '#' are comments; one of them (or
//     the first line of a CSV file) is the header that names the columns and
//     states their units in brackets:
//         # z[mm]  n[cm^-3]              n(z):    one row per point, z increasing
//         # z[mm]  r[um]  n[cm^-3]       n(r, z): one row per grid point, any order;
//                                                 the rows must fill the grid z x r
//     Columns are separated by blanks, tabs, commas or semicolons. Length units
//     as in the deck (nm, um, mm, cm, m, k0^-1, ...). Density: a physical unit
//     (cm^-3, m^-3, n_c) = the absolute density; no unit, [1] or [-] = a factor
//     on the `density` of the block.
//  2. openPMD (HDF5): a scalar mesh record (default name "density") with the
//     axes z, or r and z (also thetaMode: the mode m = 0 is taken). Regular
//     axes (gridSpacing, gridGlobalOffset, position) or explicit coordinates
//     (attribute gridCoordinates, docs/openpmd_rectilinear.md). unitDimension
//     L^-3: absolute density (unitSI -> m^-3); all zero: a factor.
//
//  Between the points the density is linear in z and r (bilinear). Inside the
//  first radius of the file it is the value there (the axis side); outside the
//  file it is zero, or the edge value with `hold`.
// ============================================================================
#pragma once
#include <memory>
#include <string>
#include <vector>

#include "units.hpp"

namespace gem {

struct DensityFile {
  std::string path;
  std::string format;         // "text" or "openPMD"
  std::string detail;         // openPMD: mesh and iteration; text: the header found
  std::string zUnit, rUnit, nUnit;  // units as stated in the file
  bool absolute = false;      // values are densities [n_c]; else a factor on the block's density
  std::vector<double> z, r;   // strictly increasing coordinates [1/k0]; r empty: n(z) only
  std::vector<double> v;      // values [iz * nr + ir], nr = max(1, r.size())
  double peak = 0.0, zPeak = 0.0, rPeak = 0.0;

  bool twoD() const { return !r.empty(); }
  size_t nr() const { return r.empty() ? 1 : r.size(); }
  // value at (r, z) of the file's coordinates; outside the file 0, or the edge value (hold)
  double at(double rr, double zz, bool hold) const;
};

// Reads a text table or an openPMD file (recognized by the HDF5 signature). `mesh` and `iteration`
// (< 0: the last) select the openPMD record. Throws std::runtime_error with the file name in the message.
std::shared_ptr<const DensityFile> readDensityFile(const std::string& path, const std::string& mesh, long iteration,
                                                   const Units& u);

}  // namespace gem
