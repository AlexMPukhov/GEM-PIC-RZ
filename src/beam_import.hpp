// ============================================================================
//  beam_import.hpp -- read a particle beam from a file (deck: bunch.file)
//
//  Formats (detected from the file):
//    openPMD 1.x (HDF5, any code: GEM-PIC-RZ, FBPIC, HiPACE++, WarpX, ...):
//        /data/<iteration>/<particlesPath>/<species>/ position, positionOffset,
//        momentum, weighting, charge, mass (datasets or constant records)
//    GEM-PIC-RZ native output data_NNNNN.h5: /particles/beam (species by the
//        bunch names, "trapped", or "inj_<ion>")
//  The particles are returned in SI-free form: positions in metres (lab z),
//  momenta u = p / (m c) with the mass m of the species, weight = number of
//  real particles, and the charge and mass of one particle. Any species with
//  one charge (a non-zero multiple of e) and one mass is accepted: electrons,
//  positrons, protons, ions in a given charge state.
// ============================================================================
#pragma once
#include <string>
#include <vector>

namespace gem {

struct ImportedParticle {
  double x, y, z;     // [m]; z: lab position
  double ux, uy, uz;  // p / (m c), m = mass of the species
  double w;           // real particles per macroparticle
};

struct ImportedBeam {
  std::vector<ImportedParticle> p;
  double charge_e = -1.0;     // charge of one particle [e]: -1 electrons, +1 positrons and protons, +Z ions
  double mass = 1.0;          // mass of one particle [m_e]
  std::string format;         // "openPMD 1.1.0", "GEM-PIC-RZ native"
  std::string species;        // species read
  long iteration = -1;        // openPMD iteration read
};

// species: empty = the only species of the file (error if there are several);
// iteration < 0: the only (or the last) iteration of the file.
// Throws std::runtime_error with a message naming the file.
ImportedBeam importBeam(const std::string& path, const std::string& species, long iteration);

}  // namespace gem
