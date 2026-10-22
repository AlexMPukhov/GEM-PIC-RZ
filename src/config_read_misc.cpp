// ============================================================================
//  config_read_misc.cpp -- block readers: output, undulator, magnet, chicane, numerics
// ============================================================================
#include "config_internal.hpp"

namespace gem {
namespace cfgdetail {

void readOutput(const Block& b, const Deck& d, const Units& u, OutputConfig& O) {
  BlockReader r(b, d, u);
  O.dir = r.str("dir", O.dir);
  O.every = int(r.integer("every", O.every));
  O.map_stride = int(r.integer("map_stride", O.map_stride));
  O.axis = r.flag("axis", O.axis);
  O.maps = r.flag("maps", O.maps);
  O.particles = r.flag("particles", O.particles);
  // format: one or more of text, hdf5, openpmd (both = text, hdf5)
  if (r.has("format")) {
    O.text = O.hdf5 = O.openpmd = false;
    const Entry* e = b.find("format");
    for (std::string w : e->tokens) {
      for (auto& ch : w) ch = char(std::tolower(static_cast<unsigned char>(ch)));
      if (w == "text") O.text = true;
      else if (w == "hdf5") O.hdf5 = true;
      else if (w == "openpmd") O.openpmd = true;
      else if (w == "both") O.text = O.hdf5 = true;
      else r.fail(e, "'" + w + "' is not allowed (choose one or more of: text, hdf5, openpmd, both)");
    }
    e->used = true;
  }
#if !defined(GEMPIC_HAVE_HDF5)
  if (O.hdf5 || O.openpmd)
    r.failKey("format", "this build has no HDF5 support (install the HDF5 library, e.g. libhdf5-dev, and rebuild)");
#endif
  // openPMD: fields interpolated to a uniform r-z grid (default: the smallest steps)
  O.openpmd_dr = r.num("openpmd_dr", Dim::Length, 0.0);
  O.openpmd_dzeta = r.num("openpmd_dzeta", Dim::Length, 0.0);
  if (O.openpmd_dr < 0) r.failKey("openpmd_dr", "must be positive");
  if (O.openpmd_dzeta < 0) r.failKey("openpmd_dzeta", "must be positive");
  const std::string og = r.word("openpmd_grid", "both", {"uniform", "native", "both"});
  O.openpmd_uniform = (og != "native");
  O.openpmd_native = (og != "uniform");
  if (!O.openpmd && r.has("openpmd_grid")) r.failKey("openpmd_grid", "needs format = openpmd");
  if (!O.openpmd && (r.has("openpmd_dr") || r.has("openpmd_dzeta")))
    r.failKey(r.has("openpmd_dr") ? "openpmd_dr" : "openpmd_dzeta", "needs format = openpmd");
  O.hdf5_single = r.word("hdf5_precision", "double", {"double", "single"}) == "single";
  O.hdf5_compression = int(r.integer("hdf5_compression", O.hdf5_compression));
  if (O.hdf5_compression < 0 || O.hdf5_compression > 9) r.failKey("hdf5_compression", "must be 0 ... 9");
  O.checkpoint_every = int(r.integer("checkpoint_every", O.checkpoint_every));
  O.checkpoint_keep = int(r.integer("checkpoint_keep", O.checkpoint_keep));
  O.walltime_hours = r.num("walltime_hours", Dim::None, O.walltime_hours);
  if (O.checkpoint_every < 0) r.failKey("checkpoint_every", "must be >= 0");
  if (O.checkpoint_keep < 0) r.failKey("checkpoint_keep", "must be >= 0");
  if (O.walltime_hours < 0) r.failKey("walltime_hours", "must be >= 0");
#if !defined(GEMPIC_HAVE_HDF5)
  if (O.checkpoint_every > 0) r.failKey("checkpoint_every", "checkpoints need HDF5 support (install libhdf5-dev and rebuild)");
  if (O.walltime_hours > 0) r.failKey("walltime_hours", "the checkpoint at the wall-time limit needs HDF5 support");
#endif
  if (O.every < 0) r.failKey("every", "must be >= 0");
  const std::string fel = r.word("fel", "auto", {"auto", "true", "false", "yes", "no", "on", "off"});
  O.fel_mode = (fel == "auto") ? -1 : (fel == "true" || fel == "yes" || fel == "on") ? 1 : 0;
  O.fel_wavelength = r.num("fel_wavelength", Dim::Length, 0.0);
  if (O.fel_wavelength < 0) r.failKey("fel_wavelength", "must be positive");
  O.fel_slice = r.num("fel_slice", Dim::Length, 0.0);
  if (O.fel_slice < 0) r.failKey("fel_slice", "must be positive");
  O.fel_harmonics = int(r.integer("fel_harmonics", 1));
  if (O.fel_harmonics < 1 || O.fel_harmonics > 9) r.failKey("fel_harmonics", "must be between 1 and 9");
  if (O.map_stride < 1) r.failKey("map_stride", "must be >= 1");
}

void readUndulator(const Block& b, const Deck& d, const Units& u, UndulatorConfig& U) {
  BlockReader r(b, d, u);
  const std::string t = r.word("type", "planar", {"planar", "helical"});
  U.type = (t == "planar") ? UndulatorConfig::Type::Planar : UndulatorConfig::Type::Helical;
  U.period = r.num("period", Dim::Length);
  if (!(U.period > 0)) r.failKey("period", "must be positive");
  U.z_start = r.num("z_start", Dim::Length, 0.0);
  if (r.oneOf("periods", "length") == 0) {
    const long n = r.integer("periods");
    if (n < 1) r.failKey("periods", "must be >= 1");
    U.length = double(n) * U.period;
  } else {
    U.length = r.num("length", Dim::Length);
    if (!(U.length > 0)) r.failKey("length", "must be positive");
  }
  // keys are case-insensitive: K, K_function
  if (r.oneOf("k", "k_function") == 0) {
    U.K = r.num("k", Dim::None);
    if (!(U.K >= 0)) r.failKey("k", "must not be negative");
  } else {
    // taper: K(z) of the lab position z (the ramps multiply it as usual)
    U.K_fn = r.expr("k_function", {"z"});
    for (int k = 0; k <= 1000; ++k) {
      const double v = U.K_fn({U.z_start + U.length * k / 1000.0});
      if (!(v >= 0) || !std::isfinite(v)) r.failKey("k_function", "must be finite and >= 0 inside the undulator");
    }
    U.K = U.K_fn({U.z_start});  // entrance value: resonance and messages
  }
  // whole periods: then the linear ramps leave no net transverse momentum (no drift of the orbit)
  const long rp = r.integer("ramp_periods", 1);
  if (rp < 0 || 2.0 * double(rp) * U.period > U.length * (1.0 + 1e-12))
    r.failKey("ramp_periods", "must be >= 0 and fit twice into the length");
  U.ramp = double(rp) * U.period;
  U.focusing = r.flag("natural_focusing", true);
  if (U.type == UndulatorConfig::Type::Helical && rp < 1)
    r.failKey("ramp_periods", "a helical undulator needs ramp_periods >= 1 (a hard edge deflects the orbit by K/gamma)");
  if (U.type == UndulatorConfig::Type::Helical && r.has("natural_focusing"))
    r.failKey("natural_focusing", "applies to planar undulators only (the helical field is the on-axis field)");
}

void readMagnet(const Block& b, const Deck& d, const Units& u, MagnetConfig& G, std::vector<std::string>& warn) {
  BlockReader r(b, d, u);
  G.label = b.label();
  G.kind = b.name == "quadrupole" ? MagnetConfig::Kind::Quadrupole
           : b.name == "plasma_lens" ? MagnetConfig::Kind::PlasmaLens : MagnetConfig::Kind::Solenoid;
  G.z_start = r.num("z_start", Dim::Length);
  G.length = r.num("length", Dim::Length);
  if (!(G.length > 0)) r.failKey("length", "must be positive");
  if (G.kind == MagnetConfig::Kind::Solenoid) {
    G.field = r.num("field", Dim::Field);
    G.ramp = r.num("ramp", Dim::Length, 0.0);
    if (G.ramp < 0 || 2.0 * G.ramp > G.length) r.failKey("ramp", "must be >= 0 and fit twice into the length");
  } else {
    G.gradient = r.num("gradient", Dim::Gradient);
    if (G.kind == MagnetConfig::Kind::PlasmaLens) {
      G.radius = r.num("radius", Dim::Length, 1e300);
      if (!(G.radius > 0)) r.failKey("radius", "must be positive");
    }
  }
  (void)warn;
}

void readChicane(const Block& b, const Deck& d, const Units& u, ChicaneConfig& C) {
  BlockReader r(b, d, u);
  C.z = r.num("z", Dim::Length);
  if (!(C.z > 0)) r.failKey("z", "must be > 0 (the chicane acts when the beam passes z)");
  C.r56 = r.num("r56", Dim::Length);
  C.delay = r.num("delay", Dim::Length, 0.5 * std::fabs(C.r56));
  if (r.oneOf("energy", "gamma") == 0) C.gamma_ref = 1.0 + r.num("energy", Dim::Energy) / MC2_MEV;
  else C.gamma_ref = r.num("gamma", Dim::None);
  if (!(C.gamma_ref > 1.0)) r.failKey(r.has("gamma") ? "gamma" : "energy", "the reference must be relativistic");
}

void readNumerics(const Block& b, const Deck& d, const Units& u, NumericsConfig& N) {
  BlockReader r(b, d, u);
#if defined(GEMPIC_HAVE_CUDA)
  const std::string bdef = "cuda";
#else
  const std::string bdef = "serial";
#endif
  const std::string be = r.word("backend", bdef, {"serial", "openmp", "cuda"});
  N.backend = be == "cuda" ? NumericsConfig::Backend::Cuda
              : be == "openmp" ? NumericsConfig::Backend::OpenMP
                               : NumericsConfig::Backend::Serial;
#if defined(GEMPIC_HAVE_CUDA)
  if (N.backend == NumericsConfig::Backend::OpenMP)
    r.failKey("backend", "this is a CUDA build: use 'cuda' or 'serial' (build without -DGEMPIC_CUDA=ON for 'openmp')");
#else
  if (N.backend == NumericsConfig::Backend::Cuda)
    r.failKey("backend", "this build has no GPU support: configure with cmake -DGEMPIC_CUDA=ON (or use 'openmp')");
#endif
  N.iterations = int(r.integer("iterations", N.iterations));
  N.beam_filter = int(r.integer("beam_filter", N.beam_filter));
  if (N.beam_filter < 0 || N.beam_filter > 8) r.failKey("beam_filter", "must be 0 ... 8");
  N.beam_vay = r.word("beam_pusher", "vay", {"vay", "boris"}) == "vay";
  N.seed = static_cast<unsigned long>(r.integer("seed", long(N.seed)));
  N.beam_dz_max = r.num("beam_dz_max", Dim::Length, N.beam_dz_max);
  N.vz_min = r.num("vz_min", Dim::None, N.vz_min);
  N.vz_beam = r.num("vz_beam", Dim::None, N.vz_beam);
  N.trap_advance = r.num("trap_advance", Dim::None, N.trap_advance);
  if (!(N.trap_advance >= 1.0)) r.failKey("trap_advance", "must be >= 1 (in steps dz)");
  N.theta = r.num("theta", Dim::None, N.theta);
  if (!(N.theta >= 0.5 && N.theta <= 1.0)) r.failKey("theta", "must be in [0.5, 1]");
  if (!(N.vz_min > 0) || !(N.vz_beam < 1) || !(N.vz_min < N.vz_beam))
    r.failKey("vz_beam", "need 0 < vz_min < vz_beam < 1");
  if (N.iterations < 1) r.failKey("iterations", "must be >= 1");
  if (!(N.beam_dz_max > 0)) r.failKey("beam_dz_max", "must be positive");
}

}  // namespace cfgdetail
}  // namespace gem
