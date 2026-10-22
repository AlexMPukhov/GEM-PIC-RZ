// ============================================================================
//  config_summary.cpp -- Config::summary, the run summary of --check
// ============================================================================
#include "config_internal.hpp"

namespace gem {
using namespace cfgdetail;

// ============================================================================
std::string Config::summary() const {
  const double um = units.k0inv_um();
  const double fs = 1.0 / units.omega0_per_fs();
  char buf[512];
  std::ostringstream os;
  auto line = [&](const char* fmt, auto... args) {
    std::snprintf(buf, sizeof buf, fmt, args...);
    os << buf << '\n';
  };
  line("units: lambda0 = %.4g um, 1/k0 = %.5g um, 1/omega0 = %.5g fs, n_c = %.4g cm^-3", units.lambda0_um, um, fs,
       units.nc_cm3());
  const auto& g = grid;
  line("grid : modes m = 0..%d", g.M);
  if (g.radial == GridConfig::Radial::Table) {
    std::string pts;
    for (size_t k = 0; k < g.rtable.points().size(); ++k) {
      std::snprintf(buf, sizeof buf, "%s(%.4g, %.4g)", k ? " " : "", g.rtable.points()[k] * um,
                    g.rtable.steps()[k] * um);
      pts += buf;
    }
    line("       radial table: Nr = %d, dr = %.4g .. %.4g um, r_max = %.4g um, max cell-to-cell change %.1f %%",
         g.Nr, g.dr_min * um, g.dr_max * um, g.r_max * um, 100 * g.rtable.maxRelativeChange());
    line("       reference points (r, dr) [um]: %s", pts.c_str());
  } else if (g.stretched)
    line("       radial stretched: Nr = %d, dr = %.4g .. %.4g um, transition at r = %.4g um (cell %.2f, width %.3g "
         "cells), r_max = %.4g um",
         g.Nr, g.dr_min * um, g.dr_max * um, g.mapping(g.xi0) * um, g.xi0, g.wcells, g.r_max * um);
  else
    line("       radial uniform: Nr = %d, dr = %.4g um, r_max = %.4g um", g.Nr, g.dr_min * um, g.r_max * um);
  auto points = [&](const StepMapping& m, double f) {
    std::string pts;
    for (size_t k = 0; k < m.points().size(); ++k) {
      std::snprintf(buf, sizeof buf, "%s(%.4g, %.4g)", k ? " " : "", m.points()[k] * f, m.steps()[k] * um);
      pts += buf;
    }
    return pts;
  };
  if (g.zetaTable) {
    line("       window table: Nz = %d, dzeta = %.4g .. %.4g um, length = %.4g um, max step-to-step change %.1f %%",
         g.Nz, g.dzeta * um, g.dzeta_max * um, g.window() * um, 100 * g.ztable.maxRelativeChange());
    line("       reference points (zeta, dzeta) [um]: %s", points(g.ztable, um).c_str());
  } else {
    line("       window: Nz = %d, dzeta = %.4g um, length = %.4g um", g.Nz, g.dzeta * um, g.window() * um);
  }
  if (g.zTable) {
    line("       z table: dz = %.4g .. %.4g um, z_end = %.4g mm (%d steps), max step-to-step change %.1f %%",
         g.dz * um, g.dz_max * um, g.z_end * um * 1e-3, g.stepsFrom(0.0), 100 * g.stable.maxRelativeChange());
    line("       reference points (z [mm], dz [um]): %s", points(g.stable, um * 1e-3).c_str());
    line("       absorber r > %.4g um, strength %.3g", g.r_absorber * um, g.absorber_strength);
  } else {
    line("       dz = %.4g um, z_end = %.4g mm (%d steps); absorber r > %.4g um, strength %.3g", g.dz * um,
         g.z_end * um * 1e-3, g.stepsFrom(0.0), g.r_absorber * um, g.absorber_strength);
  }
  for (size_t i = 0; i < lasers.size(); ++i) {
    const auto& L = lasers[i];
    // peak power, linear polarization: P = 21.5 GW (a0 w0 / lambda)^2; circular: x2
    const double lam_um = units.lambda0_um / L.omega;
    double P_TW = 21.5e-3 * std::pow(L.a0 * L.waist * um / lam_um, 2.0);
    if (L.pol == LaserConfig::Pol::Circular) P_TW *= 2.0;
    line("laser[%zu]: a0 = %.4g, lambda = %.4g um, w0 = %.4g um (zR = %.4g mm), tau_FWHM = %.4g fs, "
         "zeta_c = %.4g um, focus at z = %.4g mm, %s, P = %.3g TW",
         i, L.a0, lam_um, L.waist * um, L.zR() * um * 1e-3, L.length_1e * std::sqrt(2.0 * std::log(2.0)) * fs,
         L.zeta_center * um, L.z_focus * um * 1e-3, L.pol == LaserConfig::Pol::Linear ? "linear" : "circular",
         P_TW);
    if (!L.profile_fn.empty())
      line("        profile_function(r) = %s (field at z = 0; w0, zR, P: estimates from its 1/e radius)",
           L.profile_fn.text().c_str());
    if (!L.envelope_fn.empty())
      line("        envelope_function(zeta) = %s (tau_FWHM: estimate)", L.envelope_fn.text().c_str());
  }
  for (const auto& A : boxAdvice())
    line("box for bunch '%s' (v < c): r_max >= %s (%s; now %s), gaps >= %s at both ends (now %s ahead, %s behind "
         "after %s slippage): window >= %s (now %s)",
         A.name.c_str(), fmtUm(A.rmax_need, units).c_str(), A.rmax_why.c_str(), fmtUm(A.rmax_now, units).c_str(),
         fmtUm(A.gap_need, units).c_str(), fmtUm(A.front_now, units).c_str(), fmtUm(A.rear_end, units).c_str(),
         fmtUm(A.slip, units).c_str(), fmtUm(A.window_need, units).c_str(), fmtUm(grid.window(), units).c_str());
  for (size_t k = 0; k < undulators.size(); ++k) {
    const auto& U = undulators[k];
    std::string res;
    for (const auto& B : bunches) {
      std::snprintf(buf, sizeof buf, "%s%s %.4g um", res.empty() ? "" : ", ", B.name.c_str(),
                    U.resonance(B.gamma) * um);
      res += buf;
    }
    line("undulator[%zu]: %s, period = %.4g mm, K = %.4g (B0 = %.4g T), z = %.4g ... %.4g mm (%.4g periods, "
         "ramps %.3g periods)%s",
         k, U.type == UndulatorConfig::Type::Planar ? "planar" : "helical", U.period * um * 1e-3, U.K,
         U.B0() * units.B0_T(), U.z_start * um * 1e-3, (U.z_start + U.length) * um * 1e-3, U.length / U.period,
         U.ramp / U.period,
         U.type == UndulatorConfig::Type::Planar ? (U.focusing ? ", natural focusing" : ", no focusing") : "");
    if (!U.K_fn.empty()) line("        taper K_function(z) = %s (K above: at the entrance)", U.K_fn.text().c_str());
    if (!res.empty()) line("        resonant wavelength (on axis): %s", res.c_str());
  }
  for (const auto& G : magnets) {
    if (G.kind == MagnetConfig::Kind::Solenoid)
      line("%s: B_z = %.4g T, z = %.4g ... %.4g mm, fringes %.4g mm", G.label.c_str(), G.field * units.B0_T(),
           G.z_start * um * 1e-3, (G.z_start + G.length) * um * 1e-3, G.ramp * um * 1e-3);
    else
      line("%s: gradient %.4g T/m (%s electrons), z = %.4g ... %.4g mm%s", G.label.c_str(),
           G.gradient * units.B0_T() / (um * 1e-6),
           G.kind == MagnetConfig::Kind::PlasmaLens ? (G.gradient > 0 ? "focuses" : "defocuses")
                                                    : (G.gradient > 0 ? "focuses in x, defocuses in y" : "focuses in y, defocuses in x"),
           G.z_start * um * 1e-3, (G.z_start + G.length) * um * 1e-3,
           G.radius < 1e299 ? (", radius " + std::to_string(G.radius * um) + " um").c_str() : "");
  }
  for (const auto& C : chicanes)
    line("chicane at z = %.4g mm: R56 = %.4g um, delay %.4g um, reference energy %.5g MeV", C.z * um * 1e-3, C.r56 * um,
         C.delay * um, (C.gamma_ref - 1.0) * MC2_MEV);
  if (plasma.present) {
    const double kp = std::sqrt(plasma.n0);
    char rbuf[64] = "open (up to the grid edge)";
    if (plasma.r_max < 1e299) std::snprintf(rbuf, sizeof rbuf, "%.4g um", plasma.r_max * um);
    line("plasma: n0 = %.4g cm^-3 (%.4g n_c), lambda_p = %.4g um, radial profile: %s, r_max = %s",
         plasma.n0 * units.nc_cm3(), plasma.n0, 2.0 * PI / kp * um,
         plasma.radial == PlasmaConfig::Radial::Uniform     ? "uniform"
         : plasma.radial == PlasmaConfig::Radial::Parabolic ? "parabolic channel"
                                                            : "table",
         rbuf);
    if (!plasma.z_points.empty()) {
      std::string s;
      for (size_t i = 0; i < plasma.z_points.size(); ++i) {
        std::snprintf(buf, sizeof buf, "%s(%.4g mm, %.3g)", i ? " " : "", plasma.z_points[i] * um * 1e-3,
                      plasma.z_values[i]);
        s += buf;
      }
      line("        longitudinal profile: %s", s.c_str());
    }
    if (!plasma.fn.empty()) line("        times density_function(x, y, z, r) = %s", plasma.fn.text().c_str());
    line("        macroparticles per cell: %d x %d, and %d x %d for r < %.4g um", plasma.ppc_r, plasma.ppc_phi,
         plasma.ppc_r_axis, plasma.ppc_phi_axis, plasma.r_axis_refine * um);
  }
  if (plasma.present) {
    if (!plasma.free_electrons) line("%s", "        no free electrons at the start (neutral gas)");
    if (plasma.smooth_length > 0.0 || plasma.smooth_length_m0 > 0.0) {
      const double am0 = plasma.smooth_length_m0 >= 0.0 ? plasma.smooth_length_m0 : plasma.smooth_length;
      line("        radial smoothing of the plasma sources over a = %.4g um (%.4g k0^-1, %.3g c/omega_p) for m >= 1, "
           "%.4g um (%.4g k0^-1) for m = 0", plasma.smooth_length * units.k0inv_um(), plasma.smooth_length,
           plasma.smooth_length * std::sqrt(plasma.n0), am0 * units.k0inv_um(), am0);
    }
    double fsum = 0.0;
    bool anyIonization = false;
    for (size_t i = 0; i < ions.size(); ++i) {
      const auto& I = ions[i];
      fsum += I.fraction();
      line("ion[%zu] '%s': Z = %d, mass = %.5g u (%.6g m_e), atom density %.4g cm^-3 (%.4g n0), %s", i,
           I.name.c_str(), I.Z, I.mass / U_OVER_ME, I.mass, I.atoms * plasma.n0 * units.nc_cm3(), I.atoms,
           I.mobile ? "mobile" : "fixed");
      if (I.own)
        line("        own profile%s%s%s (atom density above: reference; with %d electrons per atom)",
             I.prof.z_points.empty() ? "" : ", longitudinal table",
             I.prof.radial == PlasmaConfig::Radial::Uniform ? "" : ", radial profile",
             I.prof.fn.empty() ? "" : (", density_function " + I.prof.fn.text()).c_str(), I.Z);
      if (I.ionization) {
        anyIonization = true;
        std::string lv;
        for (int k = I.Z; k < I.Zmax; ++k) {
          std::snprintf(buf, sizeof buf, "%s%d->%d: %.6g", k > I.Z ? ", " : "", k, k + 1, I.ip_eV[k]);
          lv += buf;
        }
        line("        ADK field ionization up to Z = %d, ionization energies [eV] %s", I.Zmax, lv.c_str());
      }
    }
    if (plasma.free_electrons && fsum < 1.0 - 1e-9)
      line("        immobile ion background: %.3g of the positive charge", 1.0 - fsum);
    if (anyIonization) line("        ionization energies: %s", atomicDataSource());
  }
  for (size_t i = 0; i < bunches.size(); ++i) {
    const auto& B = bunches[i];
    line("bunch[%zu] '%s': %s, Q = %.4g pC, E = %.5g MeV (gamma %.5g), dE/E = %.3g, sigma_x,y = %.4g, %.4g um, "
         "eps_n = %.4g, %.4g um, %s %.4g um at zeta = %.4g um, %d macroparticles",
         i, B.name.c_str(), B.qm < 0 ? "electrons" : "positrons", B.charge_pC, (B.gamma - 1.0) * MC2_MEV, B.gamma,
         B.energy_spread, B.sigma_x * um, B.sigma_y * um, B.emittance_x * um, B.emittance_y * um,
         B.shape == BunchConfig::Shape::Gaussian || B.shape == BunchConfig::Shape::File ? "sigma_z" : "length",
         (B.shape == BunchConfig::Shape::Gaussian || B.shape == BunchConfig::Shape::File ? B.sigma_z : B.length) * um,
         (B.shape == BunchConfig::Shape::File ? B.file_center : B.zeta_center) * um,
         B.macroparticles);
    if (B.shape == BunchConfig::Shape::File)
      line("        from %s (rms length %.4g um, extent %.4g ... %.4g um)", B.file_info.c_str(), B.sigma_z * um,
           B.zeta_min * um, B.zeta_max * um);
    if (B.fawley)
      line("        shot noise: Fawley beamlets of %d macroparticles at %.4g um", B.beamlet, B.noise_wavelength * um);
    if (B.shape == BunchConfig::Shape::Function)
      line("        profile_z(zeta) = %s on %.4g ... %.4g um", B.profile_z.text().c_str(), B.zeta_min * um,
           B.zeta_max * um);
    if (!B.gamma_fn.empty())
      line("        gamma_function(zeta) = %s (gamma above: mean over the profile)", B.gamma_fn.text().c_str());
  }
  if (!constants.empty()) {
    std::string s;
    for (const auto& kv : constants) {
      std::snprintf(buf, sizeof buf, "%s%s = %.6g", s.empty() ? "" : ", ", kv.first.c_str(), kv.second);
      s += buf;
    }
    line("constants (code units): %s", s.c_str());
  }
  if (output.hdf5 || output.openpmd) {
    std::string f = output.text ? "text" : "";
    if (output.hdf5) f += f.empty() ? "hdf5" : " + hdf5";
    if (output.openpmd) f += f.empty() ? "openpmd" : " + openpmd";
    line("output format: %s (HDF5 %s precision, compression %d)", f.c_str(), output.hdf5_single ? "single" : "double",
         output.hdf5_compression);
  }
  if (flying_focus.present) {
    const auto& F = flying_focus;
    line("flying focus: luminal Bessel pulse, a0 = %.4g on the axis, cone angle %.4g mrad (core radius %.4g um, "
         "transverse wavelength %.4g um), length_1e = %.4g um at zeta = %.4g um; fed through the absorber layer r > %.4g um",
         F.a0, F.theta * 1e3, 2.405 / std::sin(F.theta) * um, 2.0 * PI / std::sin(F.theta) * um, F.length_1e * um,
         F.zeta_center * um, grid.r_absorber * um);
    if (F.background > 0) line("        incident field in a uniform plasma of %.4g n_c", F.background);
  }
  if (output.fel)
    line("FEL diagnostics: fel.txt every step, bunching at %.5g um", output.fel_wavelength * um);
  if (output.fel && output.fel_harmonics > 1)
      line("        fel_harmonics.txt: harmonics 1 ... %d, modes m = 0 ... %d, slices of %.4g um", output.fel_harmonics,
           grid.M, output.fel_slice * um);
  if (output.openpmd) {
    const double dr = output.openpmd_dr > 0 ? output.openpmd_dr : grid.dr_min;
    double dzt = output.openpmd_dzeta;
    if (!(dzt > 0)) {
      dzt = 1e300;
      for (int n = 0; n < grid.Nz; ++n) dzt = std::min(dzt, grid.zetaNode(n + 1) - grid.zetaNode(n));
    }
    const long nr = long(std::floor(grid.r_max / dr + 1e-9)) + 1, nz = long(std::floor(grid.window() / dzt + 1e-9)) + 1;
    if (output.openpmd_uniform)
      line("openPMD: %s/openpmd/openpmd_%%T.h5, fields on a uniform grid dr = %.4g um, dzeta = %.4g um (%ld x %ld "
           "points, %d modes)%s",
           output.dir.c_str(), dr * um, dzt * um, nr, nz, 2 * grid.M + 1,
           output.openpmd_native ? ", and on the native grid (records native_*)" : "");
    else
      line("openPMD: %s/openpmd/openpmd_%%T.h5, fields on the native grid (records native_*, %d modes)",
           output.dir.c_str(), 2 * grid.M + 1);
  }
  if (output.checkpoint_every > 0 || output.walltime_hours > 0)
    line("checkpoints: every %d steps (keep %d)%s", output.checkpoint_every, output.checkpoint_keep,
         output.walltime_hours > 0 ? (", stop at " + std::to_string(output.walltime_hours) + " h wall time").c_str()
                                   : "");
  line("output: dir '%s', every %d steps; numerics: %d iterations/slice, theta %.3g, trap advance %.3g dz, beam "
       "filter %d, beam pusher %s, seed %lu, backend %s",
       output.dir.c_str(), output.every, numerics.iterations, numerics.theta, numerics.trap_advance, numerics.beam_filter,
       numerics.beam_vay ? "vay" : "boris", numerics.seed,
       numerics.backend == NumericsConfig::Backend::Cuda     ? "cuda"
       : numerics.backend == NumericsConfig::Backend::OpenMP ? "openmp"
                                                             : "serial");
  return os.str();
}

}  // namespace gem
