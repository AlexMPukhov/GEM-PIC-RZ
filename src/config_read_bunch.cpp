// ============================================================================
//  config_read_bunch.cpp -- block readers: shot noise, bunch from a file, bunch
// ============================================================================
#include "config_internal.hpp"

namespace gem {
namespace cfgdetail {

// shot-noise model (all bunch kinds)
void readShotNoise(const BlockReader& r, BunchConfig& B) {
  B.fawley = r.word("shot_noise", "random", {"random", "fawley"}) == "fawley";
  if (B.fawley) {
    B.noise_wavelength = r.num("noise_wavelength", Dim::Length);
    if (!(B.noise_wavelength > 0)) r.failKey("noise_wavelength", "must be positive");
    B.beamlet = int(r.integer("beamlet", 4));
    if (B.beamlet < 3) r.failKey("beamlet", "must be >= 3 (macroparticles per beamlet)");
    if (B.macroparticles % B.beamlet != 0)
      r.failKey("macroparticles", "must be a multiple of beamlet = " + std::to_string(B.beamlet) +
                                      " with shot_noise = fawley");
  } else {
    if (r.has("noise_wavelength")) r.failKey("noise_wavelength", "is used with shot_noise = fawley");
    if (r.has("beamlet")) r.failKey("beamlet", "is used with shot_noise = fawley");
  }
}

// bunch read from a file (bunch.file): select, place at zeta_center, rescale charge
void readBunchFile(const Block& b, const Deck& d, const Units& u, const GridConfig& g, BunchConfig& B,
                   std::vector<std::string>& warn) {
  BlockReader r(b, d, u);
  B.shape = BunchConfig::Shape::File;
  std::string path = r.str("file", "");
  if (!std::filesystem::exists(path)) {  // relative to the directory of the deck
    const std::filesystem::path alt = std::filesystem::path(d.source).parent_path() / path;
    if (std::filesystem::exists(alt)) path = alt.string();
  }
  ImportedBeam imp;
  try {
    imp = importBeam(path, r.str("file_species", ""), r.integer("file_iteration", -1));
  } catch (const std::exception& e) {
    r.failKey("file", e.what());
  }
  const double L = u.k0inv_um() * 1e-6;  // code length [m]
  const double Qc = 1.602176634e-19 / (u.codeToPC(1.0) * 1e-12);  // one electron in code charge units
  // longitudinal coordinate relative to the charge centroid of the species
  double sw = 0.0, sz = 0.0;
  for (const auto& p : imp.p) { sw += p.w; sz += p.w * p.z; }
  const double zc = sz / sw;
  Expr sel;
  if (r.has("select")) sel = r.expr("select", {"x", "y", "z", "ux", "uy", "uz", "gamma"});
  std::vector<BunchConfig::FileParticle> v;
  v.reserve(imp.p.size());
  for (const auto& p : imp.p) {
    const BunchConfig::FileParticle f{p.x / L, p.y / L, (p.z - zc) / L, p.ux, p.uy, p.uz, p.w * Qc};
    if (!sel.empty()) {
      const double vars[7] = {f.x, f.y, f.zeta, f.ux, f.uy, f.uz, std::sqrt(1.0 + f.ux * f.ux + f.uy * f.uy + f.uz * f.uz)};
      if (sel(vars) == 0.0) continue;
    }
    v.push_back(f);
  }
  if (v.empty()) r.failKey("select", "no particle of the file passes the selection");
  // place the centroid of the selected particles at zeta_center (plus x0, y0)
  B.x0 = r.num("x0", Dim::Length, 0.0);
  B.y0 = r.num("y0", Dim::Length, 0.0);
  const double zcen = r.num("zeta_center", Dim::Length);
  double Q = 0.0, mz = 0.0;
  for (const auto& f : v) { Q += f.q; mz += f.q * f.zeta; }
  mz /= Q;
  double scale = 1.0;
  if (r.has("charge")) {
    const double q = std::fabs(r.num("charge", Dim::Charge));
    if (!(q > 0)) r.failKey("charge", "must not be zero");
    scale = q / Q;
  }
  for (auto& f : v) { f.zeta += zcen - mz; f.x += B.x0; f.y += B.y0; f.q *= scale; }
  Q *= scale;
  // statistics for the summary, the checks and the self-fields
  double g1 = 0, x1 = 0, x2 = 0, y1 = 0, y2 = 0, z2 = 0, zmin = 1e300, zmax = -1e300;
  double xu = 0, ux1 = 0, ux2 = 0, yu = 0, uy1 = 0, uy2 = 0;
  for (const auto& f : v) {
    const double gm = std::sqrt(1.0 + f.ux * f.ux + f.uy * f.uy + f.uz * f.uz), w = f.q / Q;
    g1 += w * gm; x1 += w * f.x; x2 += w * f.x * f.x; y1 += w * f.y; y2 += w * f.y * f.y;
    z2 += w * (f.zeta - zcen) * (f.zeta - zcen); ux1 += w * f.ux; ux2 += w * f.ux * f.ux; xu += w * f.x * f.ux;
    uy1 += w * f.uy; uy2 += w * f.uy * f.uy; yu += w * f.y * f.uy;
    zmin = std::min(zmin, f.zeta); zmax = std::max(zmax, f.zeta);
  }
  double g2 = 0;
  for (const auto& f : v) {
    const double gm = std::sqrt(1.0 + f.ux * f.ux + f.uy * f.uy + f.uz * f.uz);
    g2 += f.q / Q * (gm - g1) * (gm - g1);
  }
  B.qm = imp.charge_sign;
  B.charge = imp.charge_sign * Q;
  B.charge_pC = u.codeToPC(Q);
  B.gamma = g1;
  if (!(B.gamma > 1.0)) r.failKey("file", "the mean gamma of the beam is not > 1");
  B.energy_spread = std::sqrt(g2) / g1;
  B.sigma_x = std::max(1e-12, std::sqrt(std::max(0.0, x2 - x1 * x1)));
  B.sigma_y = std::max(1e-12, std::sqrt(std::max(0.0, y2 - y1 * y1)));
  B.sigma_z = std::max(1e-12, std::sqrt(z2));
  const double cx = x2 - x1 * x1, cu = ux2 - ux1 * ux1, cxu = xu - x1 * ux1;
  B.emittance_x = std::sqrt(std::max(0.0, cx * cu - cxu * cxu));
  const double cy = y2 - y1 * y1, cv = uy2 - uy1 * uy1, cyv = yu - y1 * uy1;
  B.emittance_y = std::sqrt(std::max(0.0, cy * cv - cyv * cyv));
  B.file_center = zcen;
  B.zeta_min = zmin;
  B.zeta_max = zmax;
  B.length = std::max(zmax - zmin, 1e-12);
  B.zeta_center = 0.5 * (zmin + zmax);  // the checks use the extent
  B.macroparticles = int(r.integer("macroparticles", long(v.size())));
  if (B.macroparticles < 1) r.failKey("macroparticles", "must be >= 1");
  B.resample_smoothing = r.num("resample_smoothing", Dim::None, 0.0);
  if (B.resample_smoothing < 0) r.failKey("resample_smoothing", "must be >= 0");
  readShotNoise(r, B);
  B.file_particles.swap(v);
  char buf[512];
  std::snprintf(buf, sizeof buf, "%s, %s '%s'%s: %zu of %zu particles%s", path.c_str(), imp.format.c_str(),
                imp.species.c_str(), imp.iteration >= 0 ? (" iteration " + std::to_string(imp.iteration)).c_str() : "",
                B.file_particles.size(), imp.p.size(), sel.empty() ? "" : (" (select: " + sel.text() + ")").c_str());
  B.file_info = buf;
  if (B.macroparticles != int(B.file_particles.size()) || B.fawley)
    B.file_info += ", resampled to " + std::to_string(B.macroparticles) + " macroparticles";
  if (B.zeta_min < 0 || B.zeta_max > g.window())
    warn.push_back(b.label() + " (" + B.name + ") extends beyond the window (particles outside are removed)");
  if (std::max(B.sigma_x, B.sigma_y) < 2.0 * g.dr_min)
    warn.push_back(b.label() + " (" + B.name + "): sigma < 2 dr_min, the bunch is not resolved radially");
}

void readBunch(const Block& b, const Deck& d, const Units& u, const GridConfig& g, BunchConfig& B,
               std::vector<std::string>& warn) {
  BlockReader r(b, d, u);
  B.name = r.str("name", "bunch" + std::to_string(b.index));
  if (r.has("file")) {
    readBunchFile(b, d, u, g, B, warn);
    return;
  }
  const std::string sp = r.word("species", "electron", {"electron", "positron"});
  B.qm = (sp == "electron") ? -1.0 : 1.0;
  // charge: magnitude from the deck, sign from the species
  const double q = std::fabs(r.num("charge", Dim::Charge));
  if (!(q > 0)) r.failKey("charge", "must not be zero");
  B.charge = (B.qm < 0) ? -q : q;
  B.charge_pC = u.codeToPC(q);
  const int eg = r.oneOf("energy", "gamma", !r.has("gamma_function"));
  if (eg == 0) {
    const double E = r.num("energy", Dim::Energy);
    B.gamma = 1.0 + E / MC2_MEV;
  } else if (eg == 1) {
    B.gamma = r.num("gamma", Dim::None);
  } else {
    B.gamma = 2.0;  // gamma_function: replaced below by the mean over the profile
  }
  if (!(B.gamma > 1.0)) r.failKey(r.has("gamma") ? "gamma" : "energy", "the bunch must be relativistic (gamma > 1)");
  if (eg >= 0 && B.gamma < 10.0)
    warn.push_back(b.label() + ": gamma < 10 - the beam push along s assumes v_z close to c");
  B.energy_spread = r.num("energy_spread", Dim::None, 0.0);
  B.sigma_x = r.num("sigma_x", Dim::Length);
  B.sigma_y = r.num("sigma_y", Dim::Length, B.sigma_x);
  B.emittance_x = r.num("emittance_x", Dim::Length, 0.0);
  B.emittance_y = r.num("emittance_y", Dim::Length, B.emittance_x);
  if (!(B.sigma_x > 0) || !(B.sigma_y > 0)) r.failKey("sigma_x", "must be positive");
  const std::string shape = r.word("distribution_z", "gaussian", {"gaussian", "flattop", "function"});
  if (shape == "gaussian") {
    B.shape = BunchConfig::Shape::Gaussian;
    B.sigma_z = r.num("sigma_z", Dim::Length);
    if (!(B.sigma_z > 0)) r.failKey("sigma_z", "must be positive");
  } else if (shape == "flattop") {
    B.shape = BunchConfig::Shape::Flattop;
    B.length = r.num("length", Dim::Length);
    if (!(B.length > 0)) r.failKey("length", "must be positive");
  } else {
    // current profile f(zeta) (any normalization) on [zeta_min, zeta_max]
    B.shape = BunchConfig::Shape::Function;
    B.profile_z = r.expr("profile_z", {"zeta"});
    B.zeta_min = r.num("zeta_min", Dim::Length);
    B.zeta_max = r.num("zeta_max", Dim::Length);
    if (!(B.zeta_max > B.zeta_min)) r.failKey("zeta_max", "must be larger than zeta_min");
    double sum = 0.0;
    for (int k = 0; k <= 1000; ++k) {
      const double v = B.profile_z({B.zeta_min + (B.zeta_max - B.zeta_min) * k / 1000.0});
      if (v < 0 || !std::isfinite(v)) r.failKey("profile_z", "must be finite and >= 0 on [zeta_min, zeta_max]");
      sum += v;
    }
    if (!(sum > 0)) r.failKey("profile_z", "is zero everywhere on [zeta_min, zeta_max]");
    B.zeta_center = 0.5 * (B.zeta_min + B.zeta_max);  // for messages and the box advice
    B.length = B.zeta_max - B.zeta_min;
  }
  if (r.has("gamma_function")) {
    if (shape != "function") r.failKey("gamma_function", "needs distribution_z = function");
    B.gamma_fn = r.expr("gamma_function", {"zeta"});
    for (int k = 0; k <= 100; ++k)
      if (!(B.gamma_fn({B.zeta_min + (B.zeta_max - B.zeta_min) * k / 100.0}) > 1.0))
        r.failKey("gamma_function", "must be > 1 on [zeta_min, zeta_max]");
    if (r.has("energy") || r.has("gamma"))
      r.failKey("gamma_function", "replaces 'energy' / 'gamma': give only one of them");
    // the bunch's gamma (self-fields, box advice, messages): mean over the current profile
    double sw = 0.0, sg = 0.0;
    for (int k = 0; k <= 1000; ++k) {
      const double z = B.zeta_min + (B.zeta_max - B.zeta_min) * k / 1000.0, w = B.profile_z({z});
      sw += w;
      sg += w * B.gamma_fn({z});
    }
    B.gamma = sg / sw;
    if (B.gamma < 10.0)
      warn.push_back(b.label() + ": gamma < 10 - the beam push along s assumes v_z close to c");
  }
  if (shape != "function") B.zeta_center = r.num("zeta_center", Dim::Length);
  B.x0 = r.num("x0", Dim::Length, 0.0);
  B.y0 = r.num("y0", Dim::Length, 0.0);
  B.z_waist = r.num("z_waist", Dim::Length, 0.0);
  B.macroparticles = int(r.integer("macroparticles", 10000));
  if (B.macroparticles < 1) r.failKey("macroparticles", "must be >= 1");
  readShotNoise(r, B);

  const double half = (B.shape == BunchConfig::Shape::Gaussian) ? 3.0 * B.sigma_z : 0.5 * B.length;
  if (B.zeta_center - half < 0 || B.zeta_center + half > g.window())
    warn.push_back(b.label() + " (" + B.name + ") extends beyond the window (particles outside are removed)");
  if (3.0 * std::max(B.sigma_x, B.sigma_y) > g.r_max)
    warn.push_back(b.label() + " (" + B.name + ") is wider than the radial box");
  if (std::max(B.sigma_x, B.sigma_y) < 2.0 * g.dr_min)
    warn.push_back(b.label() + " (" + B.name + "): sigma < 2 dr_min, the bunch is not resolved radially");
}

}  // namespace cfgdetail
}  // namespace gem
