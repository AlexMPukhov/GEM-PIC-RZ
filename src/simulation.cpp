// ============================================================================
//  simulation.cpp -- initialization, particle push, deposition, time step
// ============================================================================
#include <functional>
#include "simulation.hpp"

#include "comm.hpp"
#include "eigen.hpp"
#include "rng.hpp"

#include <algorithm>
#include <limits>
#if defined(_OPENMP)
#include <omp.h>
#endif
#include <cmath>

namespace gem {

static const cd I1(0.0, 1.0);
enum { TAG_BF_FWD = 201, TAG_BF_BWD = 202, TAG_BF_GHOST = 203 };  // MPI: bunch fields at the start

// ============================================================================
//  setup
// ============================================================================
Simulation::Simulation(const Config& c, const std::string& restartFile, bool useEngine, const StorageBlock* sb)
    : cfg(c) {
  M = cfg.grid.M;
  Nz = cfg.grid.Nz;
  blk.a = 0; blk.b = Nz; blk.sa = 0; blk.sb = Nz;
  if (sb) { blk = *sb; partial = true; }
  sLo = std::max(0, blk.sa);
  sHi = std::min(Nz, blk.sb + 1);
  {
    std::vector<double> nodes(Nz + 1);
    for (int n = 0; n <= Nz; ++n) nodes[n] = cfg.grid.zetaNode(n);
    zg.build(nodes, !cfg.grid.zetaTable, cfg.grid.dzeta);
    zv = zg.view();
  }
  Delta = cfg.grid.stepAt(0.0);
  for (const auto& U : cfg.undulators) {
    UndulatorV v{U.type == UndulatorConfig::Type::Helical ? 1 : 0, U.focusing ? 1 : 0, U.ku(), U.B0(),
                 U.z_start, U.z_start + U.length, U.ramp, 0, 0, 0.0};
    if (!U.K_fn.empty()) {  // taper: table of B0(z) = K(z) k_u, 32 points per period
      v.nk = int(std::min(4.0e6, std::ceil(32.0 * U.length / U.period))) + 1;
      v.koff = int(undK.size());
      v.dzk = U.length / (v.nk - 1);
      for (int j = 0; j < v.nk; ++j) undK.push_back(U.K_fn({U.z_start + j * v.dzk}) * U.ku());
    }
    undV.push_back(v);
  }
  for (const auto& G : cfg.magnets)
    magV.push_back(MagnetV{G.kind == MagnetConfig::Kind::Quadrupole ? 0 : G.kind == MagnetConfig::Kind::PlasmaLens ? 1 : 2,
                           G.z_start, G.z_start + G.length,
                           G.kind == MagnetConfig::Kind::Solenoid ? G.field : G.gradient, G.radius, G.ramp});
  ext = ExternalV{int(undV.size()), undV.data(), undK.data(), int(magV.size()), magV.data()};
  n_ref = cfg.plasma.present ? cfg.plasma.n0 : 1.0;
  rng.seed(cfg.numerics.seed);
  g.build(cfg.grid);
  // fields and currents: the whole window, or (MPI) only the stored block, built there
  fold.initRange(M, Nz, g.N, blk.sa, blk.sb + 1);
  fnew.initRange(M, Nz, g.N, blk.sa, blk.sb + 1);
  for (auto* C : {&c0b, &c0s, &c1b, &c1s}) C->initRange(M, Nz, g.N, blk.sa, blk.sb);
  solver.setup(g, M);
  solver.theta = cfg.numerics.theta;
  sfilter.setup(g, solver, M, cfg.plasma.present ? cfg.plasma.smooth_length : 0.0,
                cfg.plasma.present ? cfg.plasma.smooth_length_m0 : 0.0);
  bfilter.setup(g, solver, M, cfg.numerics.beam_smooth_length);
  const size_t nf = size_t(M + 1) * (g.N + 1), nc = size_t(M + 1) * g.N;
  predR.assign(nf, 0.0); predRprev.assign(nf, 0.0);
  predP.assign(nc, 0.0); predPprev.assign(nc, 0.0);
  Jr.assign(g.N + 1, 0.0); Jp.assign(g.N, 0.0);
  jr1.assign(g.N + 1, 0.0); jp1.assign(g.N, 0.0); jz1.assign(g.N, 0.0);
  // field ionization
  adk.resize(cfg.ions.size());
  exitStates.assign(cfg.ions.size(), std::vector<double>());
  freedCharge.assign(cfg.ions.size(), 0.0);
  for (size_t s = 0; s < cfg.ions.size(); ++s) {
    const auto& I = cfg.ions[s];
    exitStates[s].assign(std::max(I.Zmax, I.Z) + 1, 0.0);
    if (I.ionization) adk[s].setup(I.ip_eV, I.Z, I.Zmax, cfg.units, 2.0 * zv.node(Nz));
  }
  if (cfg.flying_focus.present) setupFlyingFocus();
  // lasers from files
  fileFields.resize(cfg.lasers.size());
  bool tref = false;
  for (size_t k = 0; k < cfg.lasers.size(); ++k) {
    const auto& L = cfg.lasers[k];
    if (!L.file) continue;
    fileFields[k] = std::make_unique<LaserFileField>(*L.file, L.fileParams(), M, g.rho[g.N]);
    if (fileFields[k]->modeLoss > 1e-4)
      std::printf("WARNING: laser[%zu]: %.3g %% of the energy of the file is in azimuthal modes beyond m_max + 1 "
                  "(not kept)\n", k, 100 * fileFields[k]->modeLoss);
    if (fileFields[k]->edgeField > 1e-2)
      std::printf("WARNING: laser[%zu]: the envelope reaches %.3g of its peak at the edge of the file (outer radius "
                  "or ends of the time axis)\n", k, fileFields[k]->edgeField);
    if (L.boundary) hasIncident = true;
    if (!tref) { lasyTref = L.t_center + L.zeta_center - L.file_z; tref = true; }
  }
  if (hasIncident) {
    const double ra = cfg.grid.r_absorber;
    iaF = 1;
    while (iaF <= g.N && g.rho[iaF] <= ra) ++iaF;
    iaC = 0;
    while (iaC < g.N && g.c[iaC] <= ra) ++iaC;
    nLF = g.N + 1 - iaF;
    nLC = g.N - iaC;
    nLo = sLo;
    nInc = sHi - sLo;
    incF.assign(size_t(M + 1) * nInc * std::max(nLF, 0), cd(0.0));
    incC.assign(size_t(M + 1) * nInc * std::max(nLC, 0), cd(0.0));
    incS = std::numeric_limits<double>::quiet_NaN();
  }
  if (restartFile.empty()) {
    initLasers();
    if (cfg.flying_focus.present)  // the incident field of level s0 in the whole box
      for (int n = sLo; n < sHi; ++n) {
        for (int i = 1; i <= g.N; ++i) fold.Rr[fold.f(1, n, i)] += ffInc(true, n, i, s0);
        for (int i = 0; i < g.N; ++i) fold.Rp[fold.c(1, n, i)] += I1 * ffInc(false, n, i, s0);
      }
    initBunches();
    initBunchFields();
  } else {
    readCheckpoint(restartFile);  // fields, currents, particles, random generator, counters
    trimTables();
  }
  beginStep();  // matrices of the first step (the engine uploads them)
  if (useEngine && cfg.numerics.backend != NumericsConfig::Backend::Serial) engine = std::make_unique<Engine>(*this);
}

Simulation::~Simulation() = default;

std::string Simulation::backendInfo() const {
  if (!engine) return "serial (reference implementation, 1 thread)";
  return std::string(cfg.numerics.backend == NumericsConfig::Backend::Cuda ? "cuda" : "openmp") + " -- " +
         Engine::describe();
}

// ----------------------------------------------------------------------------
// Lasers live in mode m = 1. With the transverse phasor P = E_x + i E_y of the
// physical field, the mode-1 characteristics of a forward wave are
//     R_r,1 = conj(P),   R_phi,1 = i conj(P),   L = 0.
// The complex paraxial Gaussian beam at distance z from focus (k = omega):
//     E = a0 omega / q * exp(-r^2 / (w0^2 q)) * exp(-(zeta-zc)^2/L^2)
//         * exp(i (omega (zeta - zc) + cep)),     q = 1 + i z / zR
// Linear polarization: P = Re(E) e^{i theta};  circular: P = E or conj(E).
// ----------------------------------------------------------------------------
void Simulation::initLasers() {
  // lasers from files: the field of level s0, interpolated at the file's plane, else propagated to s0
  for (size_t k = 0; k < cfg.lasers.size(); ++k) {
    if (!fileFields[k]) continue;
    const LaserFileField& F = *fileFields[k];
    const bool direct = F.atFilePlane(s0);
    const int nz = sHi - sLo;  // the stored slices
    std::vector<double> zeta(nz), rf(g.rho.begin() + 1, g.rho.begin() + g.N + 1), rc(g.c.begin(), g.c.begin() + g.N);
    for (int n = 0; n < nz; ++n) zeta[n] = zv.node(sLo + n);
    std::vector<cd> out;
    F.field(s0, rf, zeta, 0, out, direct);
    for (int i = 1; i <= g.N; ++i)
      for (int n = 0; n < nz; ++n)
        for (int m = 0; m <= M; ++m) fold.Rr[fold.f(m, sLo + n, i)] += out[(size_t(i - 1) * nz + n) * (M + 1) + m];
    F.field(s0, rc, zeta, 1, out, direct);
    for (int i = 0; i < g.N; ++i)
      for (int n = 0; n < nz; ++n)
        for (int m = 0; m <= M; ++m) fold.Rp[fold.c(m, sLo + n, i)] += out[(size_t(i) * nz + n) * (M + 1) + m];
    std::printf("laser[%zu]: field of the file %s at z = %.6g mm\n", k,
                direct ? "interpolated" : "propagated from its plane", s0 * cfg.units.k0inv_um() * 1e-3);
    for (const auto& l : F.info) std::printf("          %s\n", l.c_str());
  }
  for (const auto& L : cfg.lasers) {
    if (L.file) continue;
    const double zd = 0.0 - L.z_focus;  // distance from focus at the start
    const double amp = L.a0 * L.omega;  // field in units m c omega0 / e
    const bool atFocus = (zd == 0.0);
    const cd q(1.0, zd / L.zR());
    const bool linear = (L.pol == LaserConfig::Pol::Linear);
    const cd rot(std::cos(L.pol_angle), -std::sin(L.pol_angle));  // e^{-i theta}

    const bool userEnv = !L.envelope_fn.empty(), userProf = !L.profile_fn.empty();
    // returns R_r,1 at (r, zeta)
    auto Rr1 = [&](double r, double zeta) -> cd {
      const double d = zeta - L.zeta_center;
      const double env = userEnv ? L.envelope_fn({zeta}) : std::exp(-d * d / (L.length_1e * L.length_1e));
      if (userProf) {
        // field g(r) env(zeta) cos(omega d + cep) at the start (no focusing applied); the phase is
        // taken from zeta_center (the centre of the envelope) as for the Gaussian pulse
        const cd E = amp * env * L.profile_fn({r}) * std::exp(cd(0.0, L.omega * d + L.cep));
        if (linear) return E.real() * rot;
        return (L.helicity > 0) ? std::conj(E) : E;
      }
      if (linear && atFocus) {
        // real arithmetic at focus (identical to the original prototype)
        const double E = amp * std::exp(-r * r / (L.waist * L.waist)) * env * std::cos(L.omega * d + L.cep);
        return (L.pol_angle == 0.0) ? cd(E, 0.0) : E * rot;
      }
      const cd T = std::exp(-r * r / (L.waist * L.waist) / q) / q;
      const cd E = amp * env * T * std::exp(cd(0.0, L.omega * d + L.cep));
      if (linear) return E.real() * rot;                     // conj(Re(E) e^{i theta})
      return (L.helicity > 0) ? std::conj(E) : E;            // conj(P)
    };
    for (int n = sLo; n < sHi; ++n) {
      const double zeta = zv.node(n);
      for (int i = 1; i <= g.N; ++i) fold.Rr[fold.f(1, n, i)] += Rr1(g.rho[i], zeta);
      for (int i = 0; i < g.N; ++i) fold.Rp[fold.c(1, n, i)] += I1 * Rr1(g.c[i], zeta);
    }
  }
}

// ----------------------------------------------------------------------------
// Bunches: Gaussian in x, y (at the waist), Gaussian or flat-top in zeta,
// uncorrelated momenta sigma_u = eps_n / sigma. If z_waist != 0, the bunch is
// drifted ballistically from the waist back to z = 0.
// ----------------------------------------------------------------------------
void Simulation::initBunches() {
  std::normal_distribution<double> nd(0.0, 1.0);
  std::uniform_real_distribution<double> ud(0.0, 1.0);
  const double Qe = 1.602176634e-19 / (cfg.units.codeToPC(1.0) * 1e-12);  // one electron [code charge]
  for (size_t b = 0; b < cfg.bunches.size(); ++b) {
    const auto& B = cfg.bunches[b];
    // Fawley beamlets: the sampled particles are the centres of beamlets of M macroparticles
    const int M = B.fawley ? B.beamlet : 1;
    const int ncent = B.macroparticles / M;
    std::vector<Particle> cent;
    cent.reserve(size_t(ncent));
    // axis refinement: the centres at r < r_a (around the axis, where the modes m >= 1 must vanish but the
    // macroparticle noise of the small cells is largest) are sampled axis_refine times as often, each with
    // 1/axis_refine of the charge: (axis_refine - 1) x (their number) more independent centres from the same
    // distribution, accepted only at r < r_a. The density is unchanged; only its noise near the axis drops.
    auto refineAxis = [&](const std::function<Particle()>& drawExtra) {
      if (B.axis_refine <= 1) return;
      const double ra = B.axis_radius > 0 ? B.axis_radius : g.rho[std::min(4, g.N)];
      const double ra2 = ra * ra;
      size_t nin = 0;
      for (auto& c : cent)
        if (c.x * c.x + c.y * c.y < ra2) { c.q /= B.axis_refine; ++nin; }
      const size_t need = nin * size_t(B.axis_refine - 1);
      const size_t maxTries = std::max<size_t>(1000000, 1000 * need);
      size_t got = 0, tries = 0;
      while (got < need && tries < maxTries) {
        Particle p = drawExtra();
        ++tries;
        if (p.x * p.x + p.y * p.y >= ra2) continue;
        p.q /= B.axis_refine;
        cent.push_back(p);
        ++got;
      }
      if (got < need)
        throw std::runtime_error("bunch '" + B.name + "': axis refinement found only " + std::to_string(got) + " of " +
                                 std::to_string(need) + " centres at r < axis_radius (is the bunch on the axis?)");
      std::printf("bunch '%s': axis refinement x%d at r < %.4g um: %zu centres there, %zu added (%zu in all)\n",
                  B.name.c_str(), B.axis_refine, ra * cfg.units.k0inv_um(), nin, got, cent.size());
    };
    auto finish = [&](Particle& p) {
      p.z = s0;
      p.qm = B.qm;
      p.status = 0;
      p.sp = int(b);
    };
    if (B.shape == BunchConfig::Shape::File) {
      const auto& F = B.file_particles;
      if (ncent == int(F.size()) && M == 1) {  // the particles of the file as they are
        for (const auto& f : F) {
          Particle p{};
          p.x = f.x; p.y = f.y; p.zeta = f.zeta; p.ux = f.ux; p.uy = f.uy; p.uz = f.uz;
          p.q = B.qm < 0 ? -f.q : f.q;
          finish(p);
          cent.push_back(p);
        }
      } else {
        // systematic resampling with probabilities ~ charge, equal weights; optional jitter
        double Q = 0.0;
        for (const auto& f : F) Q += f.q;
        double sd[6] = {0, 0, 0, 0, 0, 0};
        if (B.resample_smoothing > 0) {
          double m1[6] = {0, 0, 0, 0, 0, 0}, m2[6] = {0, 0, 0, 0, 0, 0};
          for (const auto& f : F) {
            const double v[6] = {f.x, f.y, f.zeta, f.ux, f.uy, f.uz};
            for (int a = 0; a < 6; ++a) { m1[a] += f.q / Q * v[a]; m2[a] += f.q / Q * v[a] * v[a]; }
          }
          for (int a = 0; a < 6; ++a) sd[a] = B.resample_smoothing * std::sqrt(std::max(0.0, m2[a] - m1[a] * m1[a]));
        }
        // one particle of the file, chosen with probability ~ charge (axis refinement, below)
        std::vector<double> cumq(F.size());
        { double c = 0.0; for (size_t i = 0; i < F.size(); ++i) { c += F[i].q; cumq[i] = c; } }
        auto drawFile = [&]() {
          const double t = ud(rng) * Q;
          const size_t i = std::min(size_t(std::upper_bound(cumq.begin(), cumq.end(), t) - cumq.begin()), F.size() - 1);
          const auto& f = F[i];
          Particle p{};
          p.x = f.x; p.y = f.y; p.zeta = f.zeta; p.ux = f.ux; p.uy = f.uy; p.uz = f.uz;
          if (B.resample_smoothing > 0) {
            p.x += sd[0] * nd(rng); p.y += sd[1] * nd(rng); p.zeta += sd[2] * nd(rng);
            p.ux += sd[3] * nd(rng); p.uy += sd[4] * nd(rng); p.uz += sd[5] * nd(rng);
          }
          p.q = (B.qm < 0 ? -Q : Q) / ncent;
          finish(p);
          return p;
        };
        const double u0 = ud(rng);
        double cum = 0.0;
        size_t j = 0;
        for (int k = 0; k < ncent; ++k) {
          const double target = (k + u0) / ncent * Q;
          while (j + 1 < F.size() && cum + F[j].q < target) cum += F[j++].q;
          const auto& f = F[j];
          Particle p{};
          p.x = f.x; p.y = f.y; p.zeta = f.zeta; p.ux = f.ux; p.uy = f.uy; p.uz = f.uz;
          if (B.resample_smoothing > 0) {
            p.x += sd[0] * nd(rng); p.y += sd[1] * nd(rng); p.zeta += sd[2] * nd(rng);
            p.ux += sd[3] * nd(rng); p.uy += sd[4] * nd(rng); p.uz += sd[5] * nd(rng);
          }
          p.q = (B.qm < 0 ? -Q : Q) / ncent;
          finish(p);
          cent.push_back(p);
        }
        refineAxis(drawFile);
      }
    } else {
      const double sux = B.emittance_x / B.sigma_x, suy = B.emittance_y / B.sigma_y;
      const double uz0 = std::sqrt(B.gamma * B.gamma - 1.0);
      // Function: inverse of the cumulative distribution of profile_z on a fine table
      // (piecewise-linear density: the inverse is exact within a table cell up to a quadratic)
      std::vector<double> zt, cdf;
      if (B.shape == BunchConfig::Shape::Function) {
        const int nt = 20000;
        zt.resize(nt + 1);
        cdf.assign(nt + 1, 0.0);
        double fprev = 0.0;
        for (int j = 0; j <= nt; ++j) {
          zt[j] = B.zeta_min + (B.zeta_max - B.zeta_min) * j / nt;
          const double f = std::max(0.0, B.profile_z({zt[j]}));
          if (j > 0) cdf[j] = cdf[j - 1] + 0.5 * (f + fprev) * (zt[j] - zt[j - 1]);
          fprev = f;
        }
        for (auto& c : cdf) c /= cdf.back();
      }
      auto sampleZeta = [&](double u) {
        const size_t j = std::min(size_t(std::upper_bound(cdf.begin(), cdf.end(), u) - cdf.begin()), cdf.size() - 1);
        const size_t i = (j == 0) ? 0 : j - 1;
        const double dc = cdf[i + 1] - cdf[i];
        return zt[i] + (dc > 0 ? (u - cdf[i]) / dc : 0.5) * (zt[i + 1] - zt[i]);
      };
      auto sample = [&]() {
        Particle p{};
        auto drawZeta = [&]() {
          return (B.shape == BunchConfig::Shape::Gaussian) ? B.zeta_center + B.sigma_z * nd(rng)
                 : (B.shape == BunchConfig::Shape::Flattop) ? B.zeta_center + B.length * (ud(rng) - 0.5)
                                                            : sampleZeta(ud(rng));
        };
        if (B.profile_r.empty()) {  // same order of the random numbers as before profile_r existed
          p.x = B.sigma_x * nd(rng);
          p.y = B.sigma_y * nd(rng);
          p.zeta = drawZeta();
        } else {
          p.zeta = drawZeta();  // rejection sampling of profile_r(x, y, r, zeta) in the disk r < r_extent
          const double R = B.r_extent;
          for (long t = 0;; ++t) {
            if (t > 100000000L)
              throw std::runtime_error("bunch '" + B.name + "': profile_r: no position accepted in 1e8 tries");
            const double x = R * (2.0 * ud(rng) - 1.0), y = R * (2.0 * ud(rng) - 1.0), rr = std::hypot(x, y);
            if (rr >= R) continue;
            const double f = B.profile_r({x, y, rr, p.zeta});
            if (f > B.profile_r_max)
              throw std::runtime_error("bunch '" + B.name + "': profile_r exceeds the bound of its table (" +
                                       std::to_string(f) + " > " + std::to_string(B.profile_r_max) +
                                       "); the profile has a peak narrower than r_extent/2000");
            if (ud(rng) * B.profile_r_max < f) { p.x = x; p.y = y; break; }
          }
        }
        p.ux = sux * nd(rng);
        p.uy = suy * nd(rng);
        if (B.gamma_fn.empty()) {
          p.uz = uz0 * (1.0 + B.energy_spread * nd(rng));
        } else {  // gamma(zeta): energy chirp; energy_spread is relative to the local value
          const double gz = B.gamma_fn({p.zeta});
          p.uz = std::sqrt(gz * gz - 1.0) * (1.0 + B.energy_spread * nd(rng));
        }
        if (!B.ux_fn.empty() || !B.uy_fn.empty()) {  // mean transverse momenta; the energy is kept
          const double rr = std::hypot(p.x, p.y);
          const double fx = B.ux_fn.empty() ? 0.0 : B.ux_fn({p.x, p.y, rr, p.zeta});
          const double fy = B.uy_fn.empty() ? 0.0 : B.uy_fn({p.x, p.y, rr, p.zeta});
          const double u2 = p.uz * p.uz - fx * fx - fy * fy - 2.0 * (fx * p.ux + fy * p.uy);
          if (!(u2 > 0))
            throw std::runtime_error("bunch '" + B.name + "': ux_function, uy_function exceed the momentum of the bunch");
          p.ux += fx;
          p.uy += fy;
          p.uz = std::sqrt(u2);
        }
        if (B.z_waist != 0.0) {  // drift back from the waist
          p.x -= p.ux / p.uz * B.z_waist;
          p.y -= p.uy / p.uz * B.z_waist;
        }
        if (B.x0 != 0.0) p.x += B.x0;
        if (B.y0 != 0.0) p.y += B.y0;
        p.q = B.charge / B.macroparticles * M;
        finish(p);
        return p;
      };
      for (int k = 0; k < ncent; ++k) cent.push_back(sample());
      refineAxis(sample);
    }
    if (M == 1) {
      beam.insert(beam.end(), cent.begin(), cent.end());
      continue;
    }
    // Fawley (PRST-AB 5, 070701, 2002): M copies at phases theta_j = 2 pi j / M (zero bunching),
    // displaced by dtheta_j = sum_h a_h cos(h theta_j) + b_h sin(h theta_j), h < M/2, with Gaussian
    // a_h, b_h of variance 2 / (h^2 N_e): the bunching <|b_h|^2> = 1/N_e of N_e electrons at random
    const double k = 2.0 * PI / B.noise_wavelength;
    const int H = (M - 1) / 2;
    for (const auto& c : cent) {
      const double Ne = std::fabs(c.q) / Qe;  // electrons in the beamlet
      double a[16] = {0}, bb[16] = {0};
      for (int h = 1; h <= H && h < 16; ++h) {
        const double sg = std::sqrt(2.0 / (h * h * Ne));
        a[h] = sg * nd(rng);
        bb[h] = sg * nd(rng);
      }
      for (int j = 0; j < M; ++j) {
        const double th = 2.0 * PI * j / M - PI * (M - 1) / M;
        double dth = 0.0;
        for (int h = 1; h <= H && h < 16; ++h) dth += a[h] * std::cos(h * th) + bb[h] * std::sin(h * th);
        Particle p = c;
        p.zeta = c.zeta + (th + dth) / k;
        p.q = c.q / M;
        beam.push_back(p);
      }
    }
  }
}

// ----------------------------------------------------------------------------
// Self-fields of the bunches at the start: the field of a charge distribution
// moving with the velocity beta of its bunch, at rest in the frame of the
// bunch. With psi = gamma phi' (phi' = potential in the rest frame),
//     -(Lap_perp + gamma^-2 d_zeta^2) psi = rho,        rho = j_z / beta,
//     E_perp = -grad_perp psi,  E_z = -gamma^-2 d_zeta psi,  B = beta z x E,
// i.e.  R = (1 + beta) E_perp,  L = (1 - beta) E_perp,  B_z = 0.
// For beta -> 1 this is the purely transverse field R = 2 E_perp with
// -Lap_perp psi = j_z. For beta < 1 the field reaches gamma^-1 times the radius
// beyond the ends of the bunch, and E_z is the longitudinal space-charge field.
// (A field that is not this co-moving field does not move with a bunch that
// slips back in zeta; it stays behind and exchanges energy with the bunch.)
//
// Discretization: Lap_perp is the operator of the field solver, symmetric with
// the weights A^c: diagonalized once per mode, S = D^1/2 (-Lap) D^-1/2 = Q
// diag(lambda) Q^T. Every radial eigenmode k then gives a tridiagonal problem
// along zeta (second difference on the slices, psi = 0 beyond the window):
//     (lambda_k - gamma^-2 D_zz) psi_k = rho_k.
// Each bunch is solved with its own gamma (bunches with gamma > 1e5: the
// limit beta = 1). The bunch currents are also deposited at the start level
// s^0, where the trapezoidal update of R needs them.
// ----------------------------------------------------------------------------
void Simulation::initBunchFields() {
  if (beam.empty()) return;
  for (const auto& p : beam) depositBeamParticle(c0b, p);
  filterBeamCurrent(c0b);  // the initial field belongs to the filtered current
  addBeamFields(1.0, true);
}

// adds sign x the co-moving self-field of the beam list to fold, per group of particles
// (bunches, and with deckGamma = false also trapped particles), with the gamma of the
// bunch block (deckGamma) or the mean gamma of the group
void Simulation::addBeamFields(double sign, bool deckGamma) {
  const int N = g.N;
  // radial eigenmodes of -Lap^c_m (computed on demand, per mode)
  std::vector<std::vector<double>> lam(M + 1), Q(M + 1);
  auto eigen = [&](int m) {
    if (!lam[m].empty()) return;
    const Tridiag& T = solver.poissonMatrix(m);
    std::vector<double> d(N), e(N, 0.0);
    for (int i = 0; i < N; ++i) d[i] = T.b[i];
    for (int i = 0; i + 1 < N; ++i) e[i] = -std::sqrt(T.cc[i] * T.a[i + 1]);  // symmetrized off-diagonal
    if (!symTridiagEigen(d, e, Q[m])) throw std::runtime_error("bunch fields: eigenvalue iteration failed");
    lam[m] = d;
  };
  // the stored slices [sLo, sHi) (whole window, or an MPI block with ghosts: the solve along zeta for v < c is then
  // distributed over the ranks, with the same arithmetic as on one rank)
  Currents C;
  C.initRange(M, Nz, N, blk.sa, blk.sb);
  const int nS = sHi - sLo;
  std::vector<cd> rho(size_t(nS) * N), psi(size_t(nS) * N), col(Nz), Er(N + 1), Ep(N);
  std::vector<double> sq(N);
  for (int i = 0; i < N; ++i) sq[i] = std::sqrt(g.Ac[i]);
  std::vector<int> groups;
  for (size_t b = 0; b < cfg.bunches.size(); ++b) groups.push_back(int(b));
  if (!deckGamma) {
    for (const auto& p : beam)
      if (p.sp < 0 && std::find(groups.begin(), groups.end(), p.sp) == groups.end()) groups.push_back(p.sp);
  }
  for (const int b : groups) {
    C.zero();
    bool any = false;
    double gsum = 0.0, wsum = 0.0;
    for (const auto& p : beam)
      if (p.sp == b) {
        depositBeamParticle(C, p);
        any = true;
        gsum += std::fabs(p.q) * gammaOf(p);
        wsum += std::fabs(p.q);
      }
    if (!any) continue;
    filterBeamCurrent(C);
    const double gam = deckGamma ? cfg.bunches[size_t(b)].gamma : gsum / wsum;
    const double beta = (gam > 1e5) ? 1.0 : std::sqrt(1.0 - 1.0 / (gam * gam));
    const double ig2 = (gam > 1e5) ? 0.0 : 1.0 / (gam * gam);
    for (int m = 0; m <= M; ++m) {
      // charge density at slice n: linear interpolation of the half-cells n-1 and n
      for (int n = sLo; n < sHi; ++n)
        for (int i = 0; i < N; ++i)
          rho[size_t(n - sLo) * N + i] =
              0.5 * (zv.wLo(n) * ((n > 0 && C.has(n - 1)) ? C.jz[C.c(m, n - 1, i)] : cd(0.0)) +
                     zv.wHi(n) * (C.has(n) ? C.jz[C.c(m, n, i)] : cd(0.0))) / beta;
      if (ig2 == 0.0) {
        for (int n = 0; n < nS; ++n) solver.poisson(m, &rho[size_t(n) * N], &psi[size_t(n) * N]);
      } else {
        eigen(m);
        const std::vector<double>& Z = Q[m];
        // rho_k(n) = sum_i Q_ik sqrt(A_i) rho_i(n)
        std::vector<cd> rk(size_t(nS) * N);
        for (int n = 0; n < nS; ++n) {
          const cd* r = &rho[size_t(n) * N];
          cd* out = &rk[size_t(n) * N];
          for (int i = 0; i < N; ++i) {
            const cd v = sq[i] * r[i];
            if (v == cd(0.0)) continue;
            const double* zi = &Z[size_t(i) * N];
            for (int k = 0; k < N; ++k) out[k] += zi[k] * v;
          }
        }
        // tridiagonal solve along zeta for every eigenmode (Thomas; with MPI the forward sweep runs from rank 0
        // upwards and the back substitution downwards, passing the boundary values of the own block [a, b))
        std::vector<double> lo(Nz), di(Nz), up(Nz);
        for (int n = 0; n < Nz; ++n) {
          const double hm = (n > 0) ? zv.len(n - 1) : zv.len(0), hp = zv.len(n);
          const double f = 2.0 * ig2 / (hm + hp);
          lo[n] = -f / hm; up[n] = -f / hp; di[n] = f / hm + f / hp;
        }
        const int A = blk.a, Bk = blk.b, P = blk.P, r = blk.r;
        const bool dist = partial && P > 1;
        std::vector<double> cp(size_t(Bk - A) * N), cpIn(N, 0.0), cpOut(N);
        std::vector<cd> colIn(N, cd(0.0)), colOut(N);
        auto R = [&](int n, int k) -> cd& { return rk[size_t(n - sLo) * N + k]; };
        if (dist && r > 0) {
          Buffer B = comm::recv(r - 1, TAG_BF_FWD);
          cpIn = B.getVec<double>();
          colIn = B.getVec<cd>();
        }
        for (int k = 0; k < N; ++k) {
          const double lk = lam[m][k];
          for (int n = A; n < Bk; ++n) {
            double& c = cp[size_t(n - A) * N + k];
            if (n == 0) {
              const double den = di[0] + lk;
              c = up[0] / den;
              R(0, k) /= den;
            } else {
              const double cprev = (n == A) ? cpIn[k] : cp[size_t(n - 1 - A) * N + k];
              const cd colprev = (n == A) ? colIn[k] : R(n - 1, k);
              const double den = di[n] + lk - lo[n] * cprev;
              c = up[n] / den;
              R(n, k) = (R(n, k) - lo[n] * colprev) / den;
            }
          }
          cpOut[k] = cp[size_t(Bk - 1 - A) * N + k];
          colOut[k] = R(Bk - 1, k);
        }
        if (dist && r < P - 1) {
          Buffer B;
          B.putVec(cpOut);
          B.putVec(colOut);
          comm::send(r + 1, TAG_BF_FWD, std::move(B));
          colIn = comm::recv(r + 1, TAG_BF_BWD).getVec<cd>();  // solution at slice b
        }
        for (int k = 0; k < N; ++k)
          for (int n = Bk - 1; n >= A; --n) {
            if (n >= Nz - 1) continue;  // the last slice of the window: no update
            const cd next = (n == Bk - 1) ? colIn[k] : R(n + 1, k);
            R(n, k) -= cp[size_t(n - A) * N + k] * next;
          }
        if (dist && r > 0) {
          std::vector<cd> first(N);
          for (int k = 0; k < N; ++k) first[k] = R(A, k);
          Buffer B;
          B.putVec(first);
          comm::send(r - 1, TAG_BF_BWD, std::move(B));
        }
        // ghost slices of the block from the neighbours (their own slices)
        if (dist) {
          auto sendRows = [&](int dest, int n0, int n1) {
            Buffer B;
            B.put<int>(n0);
            B.put<int>(n1);
            B.putRaw(&rk[size_t(n0 - sLo) * N], sizeof(cd) * size_t(n1 - n0) * N);
            comm::send(dest, TAG_BF_GHOST, std::move(B));
          };
          auto recvRows = [&](int src) {
            Buffer B = comm::recv(src, TAG_BF_GHOST);
            const int n0 = B.get<int>(), n1 = B.get<int>();
            std::vector<cd> v(size_t(n1 - n0) * N);
            B.getRaw(v.data(), sizeof(cd) * v.size());
            for (int n = std::max(n0, sLo); n < std::min(n1, sHi); ++n)
              std::copy_n(&v[size_t(n - n0) * N], N, &rk[size_t(n - sLo) * N]);
          };
          if (r > 0) sendRows(r - 1, A, std::min(Bk, A + blk.g));
          if (r < P - 1) sendRows(r + 1, std::max(A, Bk - blk.g), Bk);
          if (r < P - 1) recvRows(r + 1);
          if (r > 0) recvRows(r - 1);
        }
        // psi_i(n) = sum_k Q_ik psi_k(n) / sqrt(A_i)
        for (int n = 0; n < nS; ++n) {
          const cd* in = &rk[size_t(n) * N];
          cd* out = &psi[size_t(n) * N];
          for (int i = 0; i < N; ++i) {
            const double* zi = &Z[size_t(i) * N];
            cd a = 0.0;
            for (int k = 0; k < N; ++k) a += zi[k] * in[k];
            out[i] = a / sq[i];
          }
        }
      }
      for (int n = sLo; n < sHi; ++n) {
        solver.gradient(m, &psi[size_t(n - sLo) * N], Er.data(), Ep.data());
        for (int i = 1; i <= N; ++i) {
          fold.Rr[fold.f(m, n, i)] += sign * ((1.0 + beta) * Er[i]);
          if (beta < 1.0) fold.Lr[fold.f(m, n, i)] += sign * ((1.0 - beta) * Er[i]);
        }
        for (int i = 0; i < N; ++i) {
          fold.Rp[fold.c(m, n, i)] += sign * ((1.0 + beta) * Ep[i]);
          if (beta < 1.0) fold.Lp[fold.c(m, n, i)] += sign * ((1.0 - beta) * Ep[i]);
        }
        if (ig2 > 0.0) {  // E_z = -gamma^-2 d_zeta psi (central difference, psi = 0 beyond the window)
          const double hm = (n > 0) ? zv.len(n - 1) : zv.len(0), hp = zv.len(n);
          for (int i = 0; i < N; ++i) {
            const cd up1 = (n + 1 < sHi) ? psi[size_t(n + 1 - sLo) * N + i] : cd(0.0);
            const cd dn1 = (n - 1 >= sLo) ? psi[size_t(n - 1 - sLo) * N + i] : cd(0.0);
            fold.Ez[fold.c(m, n, i)] -= sign * (ig2 * (up1 - dn1) / (hm + hp));
          }
        }
      }
      if (m == 0) addAzimuthalFields(C, sign, beta);
    }
  }
}

// Co-moving field of the azimuthal current j_phi (mode 0) of a bunch, e.g. a beam rotating on helical orbits
// (bunch.ux_function, uy_function): the field of a magnetized ring moving with beta -> 1. With the vector
// potential A_phi(zeta, r) of the current,
//     B_z = curl_z A,  dB_z/dr = -j_phi (Ampere),  E_phi = d_zeta A_phi,  B_r = -beta E_phi,
// i.e. R_phi = (1 + beta) E_phi, L_phi = (1 - beta) E_phi. Discretized with the operators of the field update so
// that the field is steady in s for a source at rest in the window: the R_phi update needs g(b) = -j_phi
// (b on the faces, b = 0 at the wall), and the B_z update (Faraday, d_zeta b = CE(E)) gives E_phi = d_zeta A with
// CE(A) = b. Without this field the beam would start with B_z = 0 and get a spurious kick in the first steps.
void Simulation::addAzimuthalFields(const Currents& C, double sign, double beta) {
  // the stored slices [sLo, sHi) (whole window, or an MPI block with ghosts; the slices are independent except for
  // d_zeta A, whose values at the outermost ghost slices are not used)
  const int N = g.N, nS = sHi - sLo;
  std::vector<cd> A(size_t(nS) * N, 0.0), b(size_t(nS) * (N + 1), 0.0);
  bool any = false;
  for (int n = sLo; n < sHi; ++n) {
    // j_phi at slice n: linear interpolation of the half-cells n-1 and n
    std::vector<cd> jp(N);
    for (int i = 0; i < N; ++i) {
      jp[i] = 0.5 * (zv.wLo(n) * ((n > 0 && C.has(n - 1)) ? C.jp[C.c(0, n - 1, i)] : cd(0.0)) +
                     zv.wHi(n) * (C.has(n) ? C.jp[C.c(0, n, i)] : cd(0.0)));
      if (jp[i] != cd(0.0)) any = true;
    }
    cd* bn = &b[size_t(n - sLo) * (N + 1)];
    bn[N] = 0.0;
    for (int i = N - 1; i >= 0; --i) bn[i] = bn[i + 1] + jp[i] * (g.rho[i + 1] - g.rho[i]);
    // A on the centres from CE(A) = b: c_0 A_0 = Af_0 b_0, c_i A_i = c_{i-1} A_{i-1} + Af_i b_i
    cd* an = &A[size_t(n - sLo) * N];
    an[0] = g.Af[0] * bn[0] / g.c[0];
    for (int i = 1; i < N; ++i) an[i] = (g.c[i - 1] * an[i - 1] + g.Af[i] * bn[i]) / g.c[i];
  }
  if (!any) return;
  for (int n = sLo; n < sHi; ++n) {
    const double hm = (n > 0) ? zv.len(n - 1) : zv.len(0), hp = zv.len(n);
    for (int i = 0; i < N; ++i) {
      const cd up1 = (n + 1 < sHi) ? A[size_t(n + 1 - sLo) * N + i] : cd(0.0);
      const cd dn1 = (n - 1 >= sLo) ? A[size_t(n - 1 - sLo) * N + i] : cd(0.0);
      const cd Ep = (up1 - dn1) / (hm + hp);
      fold.Rp[fold.c(0, n, i)] += sign * ((1.0 + beta) * Ep);
      fold.Lp[fold.c(0, n, i)] += sign * ((1.0 - beta) * Ep);
    }
    for (int i = 0; i <= N; ++i) fold.Bz[fold.f(0, n, i)] += sign * b[size_t(n - sLo) * (N + 1) + i];
  }
}

// ============================================================================
//  field interpolation
// ============================================================================
// add w * (mode-m field components at radius r, slice n) to the accumulators
void Simulation::accumulate(const Fields& F, int m, int n, double w, double r, const RadialGrid::W& wc,
                            const RadialGrid::W& wf, cd& er, cd& ep, cd& ez, cd& br, cd& bp, cd& bz) const {
  const int N = g.N;
  const size_t cb = F.c(m, n, 0), fb = F.f(m, n, 0);
  auto EpC = [&](int i) { return 0.5 * (F.Rp[cb + i] + F.Lp[cb + i]); };
  auto BrC = [&](int i) { return 0.5 * (F.Lp[cb + i] - F.Rp[cb + i]); };
  auto EzC = [&](int i) { return F.Ez[cb + i]; };
  // linear extrapolation of a centre quantity to the axis
  const double c0 = g.c[0], c1 = g.c[1];
  auto axis = [&](cd f0, cd f1) { return (c1 * f0 - c0 * f1) / (c1 - c0); };
  cd ep0 = 0.0, br0 = 0.0;
  if (m == 1) { ep0 = axis(EpC(0), EpC(1)); br0 = axis(BrC(0), BrC(1)); }

  // centres
  cd vep, vbr, vez;
  if (wc.mode == 0) {
    vep = wc.wa * EpC(wc.a) + wc.wb * EpC(wc.a + 1);
    vbr = wc.wa * BrC(wc.a) + wc.wb * BrC(wc.a + 1);
    vez = wc.wa * EzC(wc.a) + wc.wb * EzC(wc.a + 1);
  } else if (wc.mode == 2) {
    vep = EpC(N - 1); vbr = BrC(N - 1); vez = EzC(N - 1);
  } else {  // r < c_0: regular behaviour on the axis
    double x = r / c0;
    vez = EzC(0) * ((m == 0) ? 1.0 : std::pow(x, m));
    if (m == 0) { vep = EpC(0) * x; vbr = BrC(0) * x; }
    else if (m == 1) { vep = ep0 + (EpC(0) - ep0) * x; vbr = br0 + (BrC(0) - br0) * x; }
    else { double fx = std::pow(x, m - 1); vep = EpC(0) * fx; vbr = BrC(0) * fx; }
  }
  // faces (face 0 for m = 1 from regularity: E_+ = B_+ = 0 on the axis)
  auto ErF = [&](int i) -> cd {
    if (i == 0) return (m == 1) ? -I1 * ep0 : cd(0.0);
    return 0.5 * (F.Rr[fb + i] + F.Lr[fb + i]);
  };
  auto BpF = [&](int i) -> cd {
    if (i == 0) return (m == 1) ? I1 * br0 : cd(0.0);
    return 0.5 * (F.Rr[fb + i] - F.Lr[fb + i]);
  };
  auto BzF = [&](int i) { return F.Bz[fb + i]; };
  cd ver, vbp, vbz;
  if (wf.mode == 0) {
    ver = wf.wa * ErF(wf.a) + wf.wb * ErF(wf.a + 1);
    vbp = wf.wa * BpF(wf.a) + wf.wb * BpF(wf.a + 1);
    vbz = wf.wa * BzF(wf.a) + wf.wb * BzF(wf.a + 1);
  } else {
    ver = ErF(N); vbp = BpF(N); vbz = BzF(N);
  }
  er += w * ver; ep += w * vep; ez += w * vez; br += w * vbr; bp += w * vbp; bz += w * vbz;
}

// Cartesian fields at (x,y), weighted combination of slices n1, n2
EB Simulation::gather(const Fields& F, int n1, double w1, int n2, double w2, double x, double y) const {
  double r = std::sqrt(x * x + y * y);
  double cs = 1.0, sn = 0.0;
  if (r > 0.0) { cs = x / r; sn = y / r; }
  auto wc = g.locC_lin(r), wf = g.locF_lin(r);  // gather: linear in r
  const cd eip(cs, sn);
  cd em(1.0, 0.0);
  double Er = 0, Ep = 0, Ez = 0, Br = 0, Bp = 0, Bz = 0;
  for (int m = 0; m <= M; ++m) {
    cd er = 0, ep = 0, ez = 0, br = 0, bp = 0, bz = 0;
    if (w1 != 0.0) accumulate(F, m, n1, w1, r, wc, wf, er, ep, ez, br, bp, bz);
    if (w2 != 0.0) accumulate(F, m, n2, w2, r, wc, wf, er, ep, ez, br, bp, bz);
    const double fac = (m == 0) ? 1.0 : 2.0;
    Er += fac * std::real(er * em); Ep += fac * std::real(ep * em);
    Ez += fac * std::real(ez * em); Br += fac * std::real(br * em);
    Bp += fac * std::real(bp * em); Bz += fac * std::real(bz * em);
    em *= eip;
  }
  return {Er * cs - Ep * sn, Er * sn + Ep * cs, Ez, Br * cs - Bp * sn, Br * sn + Bp * cs, Bz};
}

// ============================================================================
//  deposition
// ============================================================================
// current density j_m = sum W v S e^{-i m phi} / (2 pi V h_k) into half-cell k
void Simulation::deposit(Currents& C, int k, double wk, double x, double y, double vx, double vy, double vz,
                         double W, bool ion, int kc) {
  if (k < 0 || k >= Nz || wk == 0.0 || W == 0.0) return;
  if (kc < 0 && !C.has(k)) return;  // outside the stored half-cells (MPI block)
  const int kd = (kc < 0) ? k : kc;  // index in C
  const int N = g.N;
  double r = std::sqrt(x * x + y * y);
  double cs = 1.0, sn = 0.0;
  if (r > 0.0) { cs = x / r; sn = y / r; }
  const double vr = cs * vx + sn * vy, vp = -sn * vx + cs * vy;
  auto wc = g.locC(r), wf = g.locF(r);
  const double norm = W * wk / (2.0 * PI * zv.len(k));
  const cd emi(cs, -sn);
  cd em(1.0, 0.0);
  for (int m = 0; m <= M; ++m) {
    // centres: j_phi, j_z
    auto addC = [&](int i, double sp, double sz) {
      double f = norm / g.Vc[i];
      C.jp[C.c(m, kd, i)] += f * sp * vp * em;
      C.jz[C.c(m, kd, i)] += f * sz * vz * em;
    };
    if (wc.mode == 0) { addC(wc.a, wc.wa, wc.wa); addC(wc.a + 1, wc.wb, wc.wb); }
    else if (wc.mode == 2) addC(N - 1, 1.0, 1.0);
    else addC(0, 1.0, (m == 0) ? 1.0 : std::pow(r / g.c[0], m));
    // faces: j_r (face 0 is not used by the solver)
    auto addF = [&](int i, double s) {
      if (i == 0) return;
      C.jr[C.f(m, kd, i)] += norm * s / g.Vf[i] * vr * em;
    };
    if (wf.mode == 0) { addF(wf.a, wf.wa); addF(wf.a + 1, wf.wb); }
    else addF(N, 1.0);
    em *= emi;
  }
  // m = 0 charge density (diagnostic), separately for ions
  std::vector<double>& dens = ion ? C.rhoi : C.rho;
  auto addR = [&](int i, double s) { dens[C.d(kd, i)] += norm * s / g.Vc[i]; };
  if (wc.mode == 0) { addR(wc.a, wc.wa); addR(wc.a + 1, wc.wb); }
  else if (wc.mode == 2) addR(N - 1, 1.0);
  else addR(0, 1.0);
}

// plasma.smooth_length: the currents and densities that the streaming particles (plasma
// electrons, ions, freed electrons) deposited into half-cell n are smoothed radially,
// (1 - a^2 Lap)^-1, per mode (smoothing.hpp). The beam list is not smoothed, nor are the
// fields acting on the particles.
void Simulation::smoothStream(Currents& C, int n) {
  if (!sfilter.on()) return;
  const int N = g.N;
  std::vector<cd> jr(size_t(N) + 1, cd(0.0)), jp(size_t(N), cd(0.0)), jz(size_t(N), cd(0.0));
  for (int m = 0; m <= M; ++m) {
    for (int i = 0; i <= N; ++i) jr[i] = C.jr[C.f(m, n, i)];
    for (int i = 0; i < N; ++i) { jp[i] = C.jp[C.c(m, n, i)]; jz[i] = C.jz[C.c(m, n, i)]; }
    sfilter.vector(m, jr.data(), jp.data());
    sfilter.scalar(m, jz.data());
    for (int i = 1; i <= N; ++i) C.jr[C.f(m, n, i)] = jr[i];
    for (int i = 0; i < N; ++i) { C.jp[C.c(m, n, i)] = jp[i]; C.jz[C.c(m, n, i)] = jz[i]; }
  }
  sfilter.scalar0(&C.rho[C.d(n, 0)]);
  sfilter.scalar0(&C.rhoi[C.d(n, 0)]);
}

// Numerical Cherenkov: the dispersion relation of the scheme has omega = 0 at
// the zeta-Nyquist wavenumber (kappa h = pi) for every transverse k, so this
// odd-even mode is resonant with every v ~ c source. Beam particles keep their
// zeta from step to step, so the Nyquist part of their current (shot noise)
// drives it secularly. Passes of the (1,2,1)/4 filter remove that part
// (response cos^2(kappa h/2)); zero beyond the ends of the window.
void Simulation::filterBeamCurrent(Currents& C) const {
  const int passes = cfg.numerics.beam_filter;
  if (passes <= 0) { smoothBeamRadially(C); return; }
  // over the stored half-cells [k0, k0 + nk): the whole window, or an MPI block with ghosts (the values within
  // `passes` of the ends of a block are then not used)
  const int k0 = C.k0, k1 = C.k0 + C.nk;
  auto column = [&](std::vector<cd>& v, int L) {
    for (int m = 0; m <= M; ++m)
      for (int i = 0; i < L; ++i)
        for (int p = 0; p < passes; ++p) {
          cd prev(0.0);
          for (int k = k0; k < k1; ++k) {
            cd& x = v[(size_t(m) * C.nk + k - k0) * L + i];
            const cd cur = x;
            const cd next = (k + 1 < k1) ? v[(size_t(m) * C.nk + k + 1 - k0) * L + i] : cd(0.0);
            x = 0.25 * prev + 0.5 * cur + 0.25 * next;
            prev = cur;
          }
        }
  };
  column(C.jr, g.N + 1);
  column(C.jp, g.N);
  column(C.jz, g.N);
  smoothBeamRadially(C);
}

// numerics.beam_smooth_length: the beam current of every half-cell smoothed in r as the plasma sources
// (SourceFilter, charge-conserving), so that thin, intense rings or sheets of beam current have a width of at least
// a. A beam structure of width w and length L has a near field that changes along s on the scale w^2 / L; a step
// Delta >> w^2 / L makes that response stiff (the field at the current sheet does not follow the moving beam).
void Simulation::smoothBeamRadially(Currents& C) const {
  if (!bfilter.on()) return;
  const int N = g.N;
  std::vector<cd> jr(size_t(N) + 1, cd(0.0)), jp(size_t(N), cd(0.0)), jz(size_t(N), cd(0.0));
  for (int n = C.k0; n < C.k0 + C.nk; ++n) {
    for (int m = 0; m <= M; ++m) {
      bool any = false;
      for (int i = 0; i <= N; ++i) { jr[i] = C.jr[C.f(m, n, i)]; any = any || jr[i] != cd(0.0); }
      for (int i = 0; i < N; ++i) {
        jp[i] = C.jp[C.c(m, n, i)]; jz[i] = C.jz[C.c(m, n, i)];
        any = any || jp[i] != cd(0.0) || jz[i] != cd(0.0);
      }
      if (!any) continue;
      bfilter.vector(m, jr.data(), jp.data());
      bfilter.scalar(m, jz.data());
      for (int i = 1; i <= N; ++i) C.jr[C.f(m, n, i)] = jr[i];
      for (int i = 0; i < N; ++i) { C.jp[C.c(m, n, i)] = jp[i]; C.jz[C.c(m, n, i)] = jz[i]; }
    }
    bfilter.scalar0(&C.rho[C.d(n, 0)]);
  }
}

// beam-type deposit (current at level s^1): W = q / v_z, linear in zeta
// between half-cell centres
void Simulation::depositBeamParticle(Currents& C, const Particle& p) {
  const int nsub = std::max(1, int(std::ceil(Delta / cfg.numerics.beam_dz_max)));
  double ux, uy, uz;
  velocityAtPosition(ext, p, Delta / nsub, ux, uy, uz);
  const double gm = std::sqrt(1.0 + ux * ux + uy * uy + uz * uz);
  double vx = ux / gm, vy = uy / gm, vz = uz / gm;
  if (vz <= 0.0) return;
  double f;
  const int j = zv.beamCell(p.zeta, f);
  double W = p.q / vz;
  const bool ion = isIon(p.sp);
  deposit(C, j - 1, 1.0 - f, p.x, p.y, vx, vy, vz, W, ion);
  deposit(C, j, f, p.x, p.y, vx, vy, vz, W, ion);
}

// ============================================================================
//  particle pushers
// ============================================================================
// beam list: leapfrog along s with sub-steps; fields of level s^0. A particle
// that becomes slower than vz_min stops at its current (z, zeta) and joins the
// streaming particles of the current march at its zeta.
void Simulation::pushBeam() {
  if (cfg.numerics.beam_field == BEAM_FIELD_MEAN) beamSave = beam;
  for (auto& p : beam) pushBeamParticle(p);
  for (auto& p : beam)
    if (p.status == 3) {
      Particle q = p;
      q.status = 0;
      pending.push_back(q);
    }
  // compaction in order; beamMap[i]: index of the particle i of s^0 in the new list, or -1
  beamMap.assign(beam.size(), -1);
  size_t k = 0;
  for (size_t i = 0; i < beam.size(); ++i)
    if (beam[i].status == 0) {
      beamMap[i] = int(k);
      beam[k++] = beam[i];
    }
  beam.resize(k);
}

// numerics.beam_field = mean: after the march of s^1 the beam list is pushed again from s^0, in the mean
// (F^0 + F^1)/2 of the fields of the two levels, and its current of s^1 is deposited anew. The fields of s^1
// are not solved again: the change of the current is of second order in the step. Particles that the first
// push lost or returned to streaming keep that outcome; a particle that only the second push would lose
// keeps its first-push state. Particles trapped during the march are deposited unfiltered, as in the march.
void Simulation::remeanBeam() {
  size_t n = 0;
  for (int k : beamMap) if (k >= 0) ++n;
  wNew = 0.5;
  for (size_t i = 0; i < beamSave.size(); ++i) {
    if (beamMap[i] < 0) continue;
    Particle q = beamSave[i];
    pushBeamParticle(q);
    if (q.status == 0) beam[size_t(beamMap[i])] = q;
  }
  wNew = 0.0;
  c1b.zero();
  for (size_t i = 0; i < n; ++i) depositBeamParticle(c1b, beam[i]);
  filterBeamCurrent(c1b);
  for (size_t i = n; i < beam.size(); ++i) depositBeamParticle(c1b, beam[i]);
}

void Simulation::pushBeamParticle(Particle& p) {
  const int nsub = std::max(1, int(std::ceil(Delta / cfg.numerics.beam_dz_max)));
  const double dz = Delta / nsub;
  const double rmax = g.rho[g.N];
  {
    for (int sub = 0; sub < nsub && p.status == 0; ++sub) {
      double gm = gammaOf(p);
      double vz = p.uz / gm;
      if (vz < cfg.numerics.vz_min) { p.status = 3; if (wNew == 0.0) ++n_restreamed; break; }
      double dt = dz / vz;
      double w;
      const int na = zv.cell(p.zeta, w);
      if (na < 0 || na >= Nz) { p.status = 2; break; }
      EB f = gather(fold, na, 1.0 - w, na + 1, w, p.x, p.y);
      if (wNew != 0.0) {  // numerics.beam_field = mean, second push: (F^0 + F^1)/2
        const EB f1 = gather(fnew, na, 1.0 - w, na + 1, w, p.x, p.y);
        const double b = 1.0 - wNew;
        f = {b * f.Ex + wNew * f1.Ex, b * f.Ey + wNew * f1.Ey, b * f.Ez + wNew * f1.Ez,
             b * f.Bx + wNew * f1.Bx, b * f.By + wNew * f1.By, b * f.Bz + wNew * f1.Bz};
      }
      addExternal(ext, p.x, p.y, p.z, f, dz);
      if (cfg.numerics.beam_vay) vay(p.ux, p.uy, p.uz, f, p.qm, dt);
      else boris(p.ux, p.uy, p.uz, f, p.qm, dt);
      gm = gammaOf(p);
      p.x += p.ux / gm * dt; p.y += p.uy / gm * dt;
      p.zeta += (p.uz / gm - 1.0) * dt;
      p.z += dz;
      if (p.zeta < 0.0 || std::hypot(p.x, p.y) >= rmax) p.status = 2;
    }
  }
}

// fresh plasma: rings at the front of the window, quiet start in r^2.
// At every position one free electron (if plasma.free_electrons) and one
// macroparticle per ion block are loaded; an ion macroparticle represents
// `atoms` times as many atoms as the electron represents electrons and has
// the charge Z * w, so the loaded plasma is neutral.
void Simulation::loadFresh() {
  stream.clear();
  const auto& P = cfg.plasma;
  if (!P.present) return;
  std::uniform_real_distribution<double> ud(0.0, 1.0);
  const double zmid = s0 + 0.5 * Delta;
  const unsigned long long idBase = (unsigned long long)(step_no + 1) << 32;  // key of the ionization RNG
  // ion blocks with their own profile are loaded wherever their density is positive
  double rlim = P.r_max;
  for (const auto& I : cfg.ions)
    if (I.own) rlim = std::max(rlim, I.prof.r_max);
  for (int i = 0; i < g.N; ++i) {
    if (g.rho[i] >= rlim) break;
    double r2a = g.rho[i] * g.rho[i], r2b = g.rho[i + 1] * g.rho[i + 1];
    const bool ax = (g.rho[i + 1] <= P.r_axis_refine);
    const int npr = ax ? std::max(P.ppc_r, P.ppc_r_axis) : P.ppc_r;
    const int nph = ax ? std::max(P.ppc_phi, P.ppc_phi_axis) : P.ppc_phi;
    for (int j = 0; j < npr; ++j) {
      double r = std::sqrt(r2a + (j + 0.5) / npr * (r2b - r2a));
      const double nAx = P.density(r, zmid);  // axisymmetric part n0 f(z) g(r)
      if (nAx <= 0.0) {  // nothing to load here unless an own ion profile is positive (axisymmetric part)
        bool need = false;
        for (const auto& I : cfg.ions) need = need || (I.own && I.prof.density(r, zmid) > 0.0);
        if (!need) continue;
      }
      const double Q0 = -PI * (r2b - r2a) * Delta / (npr * nph);  // free-electron macro charge per unit density
      double ph0 = 2.0 * PI * ud(rng);
      for (int k = 0; k < nph; ++k) {
        double ph = ph0 + 2.0 * PI * k / nph;
        Particle p{};
        p.x = r * std::cos(ph); p.y = r * std::sin(ph); p.z = s0;
        p.qm = -1.0; p.zeta = zv.node(Nz); p.status = 0; p.sp = SPECIES_PLASMA;
        // density_function(x, y, z, r): the weight of each particle (may break the axial symmetry)
        const double n = (nAx <= 0.0) ? 0.0 : P.fn.empty() ? nAx : P.density(p.x, p.y, zmid);
        const double Q = n * Q0;
        if (n > 0.0) {
          p.q = Q; p.w = -Q;
          if (P.free_electrons) { p.id = idBase + stream.size(); stream.push_back(p); }
        }
        for (size_t s = 0; s < cfg.ions.size(); ++s) {
          const auto& I = cfg.ions[s];
          Particle ion = p;
          if (I.own) {
            const double na = I.prof.density(p.x, p.y, zmid);  // atom density
            if (na <= 0.0) continue;
            ion.w = -na * Q0;
          } else {
            if (n <= 0.0) continue;
            ion.w = -Q * I.atoms;
          }
          if (I.own && I.Z > 0) {  // the electrons of the pre-ionized atoms (before the ions, as above)
            Particle e = p;
            e.q = -I.Z * ion.w; e.w = I.Z * ion.w;
            e.id = idBase + stream.size();
            stream.push_back(e);
          }
          ion.Zi = I.Z;
          ion.q = I.Z * ion.w;
          ion.qm = I.qm();
          ion.sp = SPECIES_ION_BASE - int(s);
          ion.id = idBase + stream.size();
          stream.push_back(ion);
        }
      }
    }
  }
}

// streaming particle from its zeta (normally zeta_{n+1}; later for particles
// that join the march inside the half-cell) to zeta_n (GT2, dt = hl/(1 - v_z))
void Simulation::pushStreamParticle(Particle& p, int n, double s1) { pushStreamParticle(p, n, s1, c1s, n, born); }

void Simulation::pushStreamSlice(int n, double s1, int threads) {
  auto one = [&](Particle& p, Currents& C, int kc, std::vector<Particle>& bo) {
    if (p.status == STATUS_WAITING) {
      if (p.zeta > zv.node(n)) p.status = 0;
      else return;
    }
    pushStreamParticle(p, n, s1, C, kc, bo);
  };
#if defined(_OPENMP)
  if (threads > 1 && stream.size() >= size_t(4 * threads)) {
    // contiguous chunks in thread order: the freed electrons keep the serial order
    const size_t np = stream.size();
    std::vector<Currents> Cb(static_cast<size_t>(threads));
    std::vector<std::vector<Particle>> bo(static_cast<size_t>(threads));
#pragma omp parallel num_threads(threads)
    {
      const int t = omp_get_thread_num(), nt = omp_get_num_threads();
      Currents& C = Cb[size_t(t)];
      C.init(M, 1, g.N);
      const size_t i0 = np * size_t(t) / size_t(nt), i1 = np * size_t(t + 1) / size_t(nt);
      for (size_t i = i0; i < i1; ++i) one(stream[i], C, 0, bo[size_t(t)]);
    }
    for (int t = 0; t < threads; ++t) {
      const Currents& C = Cb[size_t(t)];
      if (C.jr.empty()) continue;
      for (int m = 0; m <= M; ++m) {
        for (int i = 0; i <= g.N; ++i) c1s.jr[c1s.f(m, n, i)] += C.jr[C.f(m, 0, i)];
        for (int i = 0; i < g.N; ++i) {
          c1s.jp[c1s.c(m, n, i)] += C.jp[C.c(m, 0, i)];
          c1s.jz[c1s.c(m, n, i)] += C.jz[C.c(m, 0, i)];
        }
      }
      for (int i = 0; i < g.N; ++i) {
        c1s.rho[c1s.d(n, i)] += C.rho[size_t(i)];
        c1s.rhoi[c1s.d(n, i)] += C.rhoi[size_t(i)];
      }
      born.insert(born.end(), bo[size_t(t)].begin(), bo[size_t(t)].end());
    }
    return;
  }
#else
  (void)threads;
#endif
  for (auto& p : stream) one(p, c1s, n, born);
}

void Simulation::pushStreamParticle(Particle& p, int n, double s1, Currents& C, int kc,
                                    std::vector<Particle>& bornOut) {
  const bool ion = isIon(p.sp);
  const IonConfig* I = (ion && p.sp != SPECIES_TRAPPED_ION) ? &cfg.ions[ionIndex(p.sp)] : nullptr;
  // length of the path in zeta within this half-cell
  const double h = zv.len(n), zn = zv.node(n);
  const double hl = (p.zeta >= zv.node(n + 1) - 1e-9 * h) ? h : std::max(0.0, p.zeta - zn);
  const double wm = 0.5 * hl / h;  // weight of slice n+1 at the midpoint of the path
  if (I && !I->mobile) {
    // fixed ion: no motion (dt = h), fields only needed for ionization
    deposit(C, n, 1.0, p.x, p.y, 0.0, 0.0, 0.0, p.q * hl / Delta, true, kc);  // density diagnostic only
    const double zs = p.zeta;
    p.zeta = zn;
    if (I->ionization && p.Zi < I->Zmax && sliceFieldBound >= adk[ionIndex(p.sp)].minField()) {
      const double w1 = (zs - zn) / h;  // field at the start of the path
      const EB f1 = gather(fnew, n + 1, w1, n, 1.0 - w1, p.x, p.y);
      const EB f0 = gather(fnew, n, 1.0, n, 0.0, p.x, p.y);
      ionize(p, std::hypot(f1.Ex, f1.Ey, f1.Ez), std::hypot(f0.Ex, f0.Ey, f0.Ez), hl, bornOut);
    }
    return;
  }
  double g0 = gammaOf(p);
  double omv0 = (g0 - p.uz) / g0;  // 1 - v_z > 0
  double dt0 = hl / omv0;
  double xm = p.x + 0.5 * dt0 * p.ux / g0, ym = p.y + 0.5 * dt0 * p.uy / g0;
  // ionizing ions need the field at both slices (the rate is averaged, see
  // ionize); all other particles use the averaged field directly
  const bool ionizing = I && I->ionization && p.Zi < I->Zmax && sliceFieldBound >= adk[ionIndex(p.sp)].minField();
  EB f, f1{}, f0{};
  if (ionizing) {
    const double w1 = 2.0 * wm;  // start of the path
    f1 = gather(fnew, n + 1, w1, n, 1.0 - w1, xm, ym);
    f0 = gather(fnew, n, 1.0, n, 0.0, xm, ym);
    f = {0.5 * (f1.Ex + f0.Ex), 0.5 * (f1.Ey + f0.Ey), 0.5 * (f1.Ez + f0.Ez),
         0.5 * (f1.Bx + f0.Bx), 0.5 * (f1.By + f0.By), 0.5 * (f1.Bz + f0.Bz)};
  } else {
    f = gather(fnew, n + 1, wm, n, 1.0 - wm, xm, ym);
  }
  if (ext.any()) addExternal(ext, xm, ym, p.z + 0.5 * dt0 * p.uz / g0, f, dt0 * p.uz / g0);
  double ux = p.ux, uy = p.uy, uz = p.uz;
  boris(ux, uy, uz, f, p.qm, dt0);
  double g1 = std::sqrt(1.0 + ux * ux + uy * uy + uz * uz);
  double dt = 2.0 * hl / (omv0 + (g1 - uz) / g1);  // time for |dzeta| = hl
  ux = p.ux; uy = p.uy; uz = p.uz;
  boris(ux, uy, uz, f, p.qm, dt);
  g1 = std::sqrt(1.0 + ux * ux + uy * uy + uz * uz);
  double vx = 0.5 * (p.ux / g0 + ux / g1), vy = 0.5 * (p.uy / g0 + uy / g1), vz = 0.5 * (p.uz / g0 + uz / g1);
  double zeta_end = zn;
  const double sTrap = s1 + (cfg.numerics.trap_advance - 1.0) * Delta;
  // advanced trap_advance steps inside this half-cell, and fast enough (trap_gamma) -> trapped
  if (p.z + vz * dt >= sTrap && g1 >= cfg.numerics.trap_gamma) {
    double dtp = std::max(0.0, (sTrap - p.z) / std::max(vz, 1e-3));
    ux = p.ux; uy = p.uy; uz = p.uz;
    boris(ux, uy, uz, f, p.qm, dtp);
    g1 = std::sqrt(1.0 + ux * ux + uy * uy + uz * uz);
    vx = 0.5 * (p.ux / g0 + ux / g1); vy = 0.5 * (p.uy / g0 + uy / g1);
    vz = 0.5 * (p.uz / g0 + uz / g1);
    dt = dtp;
    zeta_end = p.zeta - (1.0 - vz) * dtp;
    p.status = 1;
  }
  double xo = p.x, yo = p.y;
  p.x += vx * dt; p.y += vy * dt;
  p.z = (p.status == 1) ? s1 : p.z + vz * dt;
  p.ux = ux; p.uy = uy; p.uz = uz;
  p.zeta = zeta_end;
  // streaming current: W = q dt / Delta, into half-cell n
  deposit(C, n, 1.0, 0.5 * (xo + p.x), 0.5 * (yo + p.y), vx, vy, vz, p.q * dt / Delta, ion, kc);
  if (std::hypot(p.x, p.y) >= g.rho[g.N]) p.status = 2;
  // field ionization with the mid-step field (after the deposit: the new
  // charge state acts from the next half-cell on)
  if (ionizing && p.status == 0)
    ionize(p, std::hypot(f1.Ex, f1.Ey, f1.Ez), std::hypot(f0.Ex, f0.Ey, f0.Ez), dt, bornOut);
}

// Upper bound of |E| anywhere on slice n: sum over the modes of the largest
// nodal amplitudes, times 2 for the linear extrapolation to the axis used by
// the gather. Used to skip ionization where no atom can be ionized.
double Simulation::fieldBound(const Fields& F, int n) const {
  double bound = 0.0;
  for (int m = 0; m <= M; ++m) {
    double er = 0, ep = 0, ez = 0;
    for (int i = 0; i <= g.N; ++i) er = std::max(er, 0.5 * std::abs(F.Rr[F.f(m, n, i)] + F.Lr[F.f(m, n, i)]));
    for (int i = 0; i < g.N; ++i) {
      ep = std::max(ep, 0.5 * std::abs(F.Rp[F.c(m, n, i)] + F.Lp[F.c(m, n, i)]));
      ez = std::max(ez, std::abs(F.Ez[F.c(m, n, i)]));
    }
    bound += (m == 0 ? 1.0 : 2.0) * std::sqrt(er * er + ep * ep + ez * ez);
  }
  return 2.0 * bound;
}

// ----------------------------------------------------------------------------
// Monte Carlo field ionization of one ion macroparticle during the time dt.
// E1, E0 = |E| at the two slices of the step. The rate of each level is the
// average (w(E1) + w(E0))/2, i.e. the trapezoidal rule for int w dt. (The rate
// of the averaged field, w((E1+E0)/2), would underestimate the ionization:
// w is very convex and the carrier peaks span only a few slices; for H near
// threshold at 25 points per wavelength the error is -6 %, with the averaged
// rate it is below 1e-4.) Sequential ionization k -> k+1 -> ... is sampled
// with exponential waiting times, so several levels can be ionized in one
// step; for rates constant over the step this reproduces the rate equations
// dN_k/dt = w_{k-1} N_{k-1} - w_k N_k exactly (in the mean).
// Every event raises the charge state of the whole macroparticle by one and
// creates an electron macroparticle (charge -w) at rest relative to the ion.
// ----------------------------------------------------------------------------
void Simulation::ionize(Particle& p, double E1, double E0, double dt, std::vector<Particle>& bornOut) {
  const int s = ionIndex(p.sp);
  const AdkRates& R = adk[s];
  if (E1 < R.minField() && E0 < R.minField()) return;
  const uint64_t base = ionKeyParticle(ionKey, p.id);  // counter-based (rng.hpp)
  unsigned draw = 0;
  double tleft = dt;
  const int Z0 = p.Zi;
  while (p.Zi < R.zmax()) {
    const double w = 0.5 * (R.rate(p.Zi, E1) + R.rate(p.Zi, E0));
    if (!(w > 0.0)) break;
    const double tw = -std::log(1.0 - uniform01(base + draw++)) / w;
    if (tw > tleft) break;
    tleft -= tw;
    ++p.Zi;
    Particle e = p;
    e.q = -p.w;
    e.qm = -1.0;
    e.sp = SPECIES_IONIZED_BASE - s;
    e.Zi = 0;
    bornOut.push_back(e);
  }
  if (p.Zi != Z0) {
    p.q = p.Zi * p.w;
    p.qm = p.Zi / cfg.ions[s].mass;
  }
}

// ============================================================================
//  one step s^0 -> s^1
// ============================================================================
void Simulation::beginStep() {
  Delta = cfg.grid.stepAt(s0);
  solver.prepare(Delta, zg);
}

// thin chicanes passed in the last step (z_a, z_b]: the beam particles move by
// r56 delta - delay in zeta (delta = gamma / gamma_ref - 1). Their co-moving
// self-field moves with them (subtracted at the old, added at the new
// positions; vacuum), and the beam current of the level s^0 is deposited anew.
void Simulation::applyChicanes(double za, double zb) {
  for (const auto& C : cfg.chicanes) {
    if (!(za < C.z && C.z <= zb)) continue;
    if (engine) engine->downloadFull();
    addBeamFields(-1.0, false);
    double sh = 0.0, sh2 = 0.0, w = 0.0;
    const double top = zv.node(Nz);
    for (auto& p : beam) {
      const double d = C.r56 * (gammaOf(p) / C.gamma_ref - 1.0) - C.delay;
      p.zeta += d;
      sh += std::fabs(p.q) * d; sh2 += std::fabs(p.q) * d * d; w += std::fabs(p.q);
      if (p.zeta < 0.0 || p.zeta >= top) p.status = 2;
    }
    const size_t n0 = beam.size();
    beam.erase(std::remove_if(beam.begin(), beam.end(), [](const Particle& p) { return p.status != 0; }), beam.end());
    addBeamFields(1.0, false);
    c0b.zero();
    for (const auto& p : beam) depositBeamParticle(c0b, p);
    filterBeamCurrent(c0b);
    if (w > 0) {
      const double um = cfg.units.k0inv_um(), m = sh / w;
      std::printf("chicane at z = %.6g mm: beam moved by %.4g um (rms spread %.3g um)%s\n", C.z * um * 1e-3, m * um,
                  std::sqrt(std::max(0.0, sh2 / w - m * m)) * um,
                  n0 > beam.size() ? (", " + std::to_string(n0 - beam.size()) + " particles left the window").c_str() : "");
    }
    if (engine) engine = std::make_unique<Engine>(*this);  // upload the changed state
  }
}

void Simulation::step() {
  beginStep();
  const double za = s0;
  if (engine) {  // parallel backend: same algorithm, state on the device
    engine->step();
    engine->download();
    if (!cfg.chicanes.empty()) applyChicanes(za, s0);
    return;
  }
  const double s1 = s0 + Delta;
  march(s1);
  if (cfg.numerics.beam_field == BEAM_FIELD_MEAN) remeanBeam();
  // charge states of the ions leaving the rear of the window
  for (auto& v : exitStates) std::fill(v.begin(), v.end(), 0.0);
  for (const auto& p : stream)
    if (p.status == 0 && isIon(p.sp) && p.sp != SPECIES_TRAPPED_ION) exitStates[ionIndex(p.sp)][p.Zi] += p.w;
  prepareIncident(s1, 0, Nz);
  absorber(s1);
  std::swap(fold, fnew);
  std::swap(c0b, c1b);
  std::swap(c0s, c1s);
  s0 = s1;
  ++step_no;
  if (!cfg.chicanes.empty()) applyChicanes(za, s0);
}

// one march of the level s^1: beam push and deposit, fresh plasma, slices from the front to the rear
void Simulation::march(double s1) {
  const int N = g.N;
  c1b.zero(); c1s.zero();
  std::fill(freedCharge.begin(), freedCharge.end(), 0.0);
  pending.swap(pendingNext);  // particles that reached s^0 in the last step with v_z < vz_beam
  pendingNext.clear();
  pushBeam();                 // adds beam particles that became slower than vz_min
  for (const auto& p : beam) depositBeamParticle(c1b, p);
  filterBeamCurrent(c1b);  // (particles trapped during the march are deposited unfiltered)
  loadFresh();
  // particles joining the march at their zeta wait in the streaming list
  for (const auto& p : pending) {
    Particle q = p;
    q.status = STATUS_WAITING;
    stream.push_back(q);
  }
  pending.clear();
  const uint64_t keyStep = ionKeyStep(cfg.numerics.seed, step_no);
  std::fill(predRprev.begin(), predRprev.end(), 0.0);
  std::fill(predPprev.begin(), predPprev.end(), 0.0);

  for (int n = Nz - 1; n >= 0; --n) {
    // prediction of j^1 at half-cell n-1: beam part (known) + streaming part of s^0
    for (int m = 0; m <= M; ++m) {
      for (int i = 0; i <= N; ++i)
        predR[m * (N + 1) + i] = (n > 0) ? c1b.jr[c1b.f(m, n - 1, i)] + c0s.jr[c0s.f(m, n - 1, i)] : 0.0;
      for (int i = 0; i < N; ++i)
        predP[m * N + i] = (n > 0) ? c1b.jp[c1b.c(m, n - 1, i)] + c0s.jp[c0s.c(m, n - 1, i)] : 0.0;
    }
    if (cfg.numerics.iterations > 1) backup = stream;
    for (int it = 0; it < cfg.numerics.iterations; ++it) {
      if (it == 0) fnew.copySlice(fold, n);  // predictor: slow evolution in s
      else stream = backup;
      born.clear();  // electrons freed by ionization: only the last pass counts
      if (!cfg.ions.empty()) sliceFieldBound = std::max(fieldBound(fnew, n + 1), fieldBound(fnew, n));
      c1s.zeroHalf(n);
      ionKey = ionKeySlice(keyStep, n, it);
      for (auto& p : stream) {
        if (p.status == STATUS_WAITING) {  // joins the march inside half-cell n
          if (p.zeta > zv.node(n)) p.status = 0;
          else continue;
        }
        pushStreamParticle(p, n, s1);
      }
      smoothStream(c1s, n);
      // field solve, all modes
      for (int m = 0; m <= M; ++m) {
        // four-point current sums, weighted by the interpolation to slice n
        // (wl = wh = 1 on a uniform grid); at the last slice (n = 0) the missing
        // half-cell -1 is replaced by half-cell 0, at s^1 with the current of
        // this pass (a zero there would halve the source of R(0): a growing spike)
        const int nm = (n > 0) ? n - 1 : 0;
        const double wl = zv.wLo(n), wh = zv.wHi(n);
        const double th1 = 2.0 * cfg.numerics.theta, th0 = 2.0 - th1;  // 1, 1 for the trapezoid
        for (int i = 0; i <= N; ++i) {
          const cd a = c0b.jr[c0b.f(m, nm, i)] + c0s.jr[c0s.f(m, nm, i)];
          const cd pr = (n > 0) ? predR[m * (N + 1) + i] : c1b.jr[c1b.f(m, 0, i)] + c1s.jr[c1s.f(m, 0, i)];
          Jr[i] = th0 * (wl * a + wh * c0b.jr[c0b.f(m, n, i)] + wh * c0s.jr[c0s.f(m, n, i)]) +
                  th1 * (wl * pr + wh * c1b.jr[c1b.f(m, n, i)] + wh * c1s.jr[c1s.f(m, n, i)]);
        }
        for (int i = 0; i < N; ++i) {
          const cd a = c0b.jp[c0b.c(m, nm, i)] + c0s.jp[c0s.c(m, nm, i)];
          const cd pr = (n > 0) ? predP[m * N + i] : c1b.jp[c1b.c(m, 0, i)] + c1s.jp[c1s.c(m, 0, i)];
          Jp[i] = th0 * (wl * a + wh * c0b.jp[c0b.c(m, n, i)] + wh * c0s.jp[c0s.c(m, n, i)]) +
                  th1 * (wl * pr + wh * c1b.jp[c1b.c(m, n, i)] + wh * c1s.jp[c1s.c(m, n, i)]);
        }
        // j^1 at half-cell n = beam + streaming
        for (int i = 0; i <= N; ++i) jr1[i] = c1b.jr[c1b.f(m, n, i)] + c1s.jr[c1s.f(m, n, i)];
        for (int i = 0; i < N; ++i) {
          jp1[i] = c1b.jp[c1b.c(m, n, i)] + c1s.jp[c1s.c(m, n, i)];
          jz1[i] = c1b.jz[c1b.c(m, n, i)] + c1s.jz[c1s.c(m, n, i)];
        }
        solver.solveSlice(m, n, &fold.Rr[fold.f(m, n, 0)], &fold.Rp[fold.c(m, n, 0)], &fold.Ez[fold.c(m, n, 0)],
                          &fold.Bz[fold.f(m, n, 0)], &fnew.Rr[fnew.f(m, n + 1, 0)], &fnew.Lr[fnew.f(m, n + 1, 0)],
                          &fnew.Rp[fnew.c(m, n + 1, 0)], &fnew.Lp[fnew.c(m, n + 1, 0)],
                          &fnew.Ez[fnew.c(m, n + 1, 0)], &fnew.Bz[fnew.f(m, n + 1, 0)], jr1.data(), jp1.data(),
                          jz1.data(), Jr.data(), Jp.data(), &fnew.Rr[fnew.f(m, n, 0)], &fnew.Lr[fnew.f(m, n, 0)],
                          &fnew.Rp[fnew.c(m, n, 0)], &fnew.Lp[fnew.c(m, n, 0)], &fnew.Ez[fnew.c(m, n, 0)],
                          &fnew.Bz[fnew.f(m, n, 0)]);
      }
    }
    // electrons freed in this slice start streaming from slice n
    for (const auto& e : born) freedCharge[ionizedIndex(e.sp)] += e.w;
    stream.insert(stream.end(), born.begin(), born.end());
    born.clear();
    // trapped particles -> beam list; lost particles removed
    for (auto& p : stream) {
      if (p.status == 1) {
        Particle t = p;
        t.status = 0;
        t.z = s1;
        if (t.uz / gammaOf(t) >= cfg.numerics.vz_beam) {
          // fast: continue on the beam list (streaming species get their trapped id)
          if (t.sp == SPECIES_PLASMA) t.sp = SPECIES_TRAPPED;
          else if (isIonizationElectron(t.sp)) t.sp = SPECIES_TRAPPED_IONIZED_BASE - ionizedIndex(t.sp);
          else if (isIon(t.sp)) t.sp = SPECIES_TRAPPED_ION;
          beam.push_back(t);
          depositBeamParticle(c1b, t);
          ++n_trapped_total;
        } else {
          // slow: continue streaming in the next step from its zeta
          pendingNext.push_back(t);
        }
      }
    }
    stream.erase(std::remove_if(stream.begin(), stream.end(),
                                [](const Particle& p) { return p.status != 0 && p.status != STATUS_WAITING; }),
                 stream.end());
    // half-cell n is now final: correct R at slice n+1 for its predicted j^1_{n+1/2}
    if (n + 1 <= Nz - 1) {
      const double q = 0.5 * cfg.numerics.theta * Delta * zv.wLo(n + 1);
      for (int m = 0; m <= M; ++m) {
        for (int i = 1; i <= N; ++i) {
          cd jt = c1b.jr[c1b.f(m, n, i)] + c1s.jr[c1s.f(m, n, i)];
          fnew.Rr[fnew.f(m, n + 1, i)] -= q * (jt - predRprev[m * (N + 1) + i]);
        }
        for (int i = 0; i < N; ++i) {
          cd jt = c1b.jp[c1b.c(m, n, i)] + c1s.jp[c1s.c(m, n, i)];
          fnew.Rp[fnew.c(m, n + 1, i)] -= q * (jt - predPprev[m * N + i]);
        }
      }
    }
    predRprev.swap(predR);
    predPprev.swap(predP);
  }
}

// damp R near the outer radial boundary; with a flying focus only the difference to the incident field
// (total-field / scattered-field: the layer feeds the converging wave and absorbs what the plasma adds)
void Simulation::absorber(double s) {
  const double ra = cfg.grid.r_absorber, rm = g.rho[g.N];
  if (ra >= rm) return;
  auto damp = [&](double r) {
    if (r <= ra) return 1.0;
    double x = (r - ra) / (rm - ra);
    return std::max(0.0, 1.0 - cfg.grid.absorber_strength * x * x);
  };
  const bool ff = cfg.flying_focus.present;
  for (int m = 0; m <= M; ++m)
    for (int n = 0; n < Nz; ++n) {
      for (int i = 0; i <= g.N; ++i) {
        const double d = damp(g.rho[i]);
        if (d == 1.0) continue;
        const cd inc = ((ff && m == 1) ? cd(ffInc(true, n, i, s)) : cd(0.0)) + incidentF(m, n, i);
        cd& R = fnew.Rr[fnew.f(m, n, i)];
        R = inc + d * (R - inc);
      }
      for (int i = 0; i < g.N; ++i) {
        const double d = damp(g.c[i]);
        if (d == 1.0) continue;
        const cd inc = ((ff && m == 1) ? I1 * ffInc(false, n, i, s) : cd(0.0)) + incidentC(m, n, i);
        cd& R = fnew.Rp[fnew.c(m, n, i)];
        R = inc + d * (R - inc);
      }
    }
}

// inject = boundary: incident field of the file lasers in the absorber layer at level s, slices [n0, n1)
void Simulation::prepareIncident(double s, int n0, int n1) {
  if (!hasIncident) return;
  n0 = std::max(sLo, n0);
  n1 = std::min(sHi, n1);
  if (s == incS && n0 >= incN0 && n1 <= incN1) return;
  std::fill(incF.begin(), incF.end(), cd(0.0));
  std::fill(incC.begin(), incC.end(), cd(0.0));
  std::vector<double> zeta(n1 - n0), rf(g.rho.begin() + iaF, g.rho.begin() + g.N + 1),
      rc(g.c.begin() + iaC, g.c.begin() + g.N);
  for (int n = n0; n < n1; ++n) zeta[n - n0] = zv.node(n);
  std::vector<cd> out;
  const int nz = n1 - n0;
  for (size_t k = 0; k < cfg.lasers.size(); ++k) {
    if (!fileFields[k] || !cfg.lasers[k].boundary) continue;
    fileFields[k]->field(s, rf, zeta, 0, out, false);
    for (int i = 0; i < nLF; ++i)
      for (int n = 0; n < nz; ++n)
        for (int m = 0; m <= M; ++m) incF[(size_t(m) * nInc + n0 + n - nLo) * nLF + i] += out[(size_t(i) * nz + n) * (M + 1) + m];
    fileFields[k]->field(s, rc, zeta, 1, out, false);
    for (int i = 0; i < nLC; ++i)
      for (int n = 0; n < nz; ++n)
        for (int m = 0; m <= M; ++m) incC[(size_t(m) * nInc + n0 + n - nLo) * nLC + i] += out[(size_t(i) * nz + n) * (M + 1) + m];
  }
  incS = s;
  incN0 = n0;
  incN1 = n1;
}

// incident field of the flying focus at s = 0: G(r, zeta) = sum_j c_j J0(k_j r) e^{i w_j (zeta - zeta_c) + i cep},
// c_j the Gaussian spectrum of the on-axis envelope, k_j = sqrt(2 w_j kappa - kappa^2) (k_z = w_j - kappa)
void Simulation::setupFlyingFocus() {
  const auto& F = cfg.flying_focus;
  const double kap = F.kappa(), Lz = F.length_1e, amp = F.a0;  // a0 omega0, omega0 = 1
  const int nw = 241;
  const double w0 = 1.0 - 6.0 / Lz, w1 = 1.0 + 6.0 / Lz, dw = (w1 - w0) / (nw - 1);
  std::vector<double> wj(nw), kj(nw);
  std::vector<cd> cj(nw);
  for (int j = 0; j < nw; ++j) {
    wj[j] = w0 + dw * j;
    kj[j] = std::sqrt(std::max(0.0, 2.0 * wj[j] * kap - kap * kap - F.background));
    // int dw (L / (2 sqrt(pi))) e^{-(w - 1)^2 L^2 / 4} e^{i w d} = e^{i d} e^{-d^2 / L^2}
    cj[j] = amp * dw * Lz / (2.0 * std::sqrt(PI)) * std::exp(-0.25 * (wj[j] - 1.0) * (wj[j] - 1.0) * Lz * Lz) *
            std::polar(1.0, F.cep);
  }
  auto fill = [&](std::vector<cd>& G, const std::vector<double>& rr, int nr, int i0) {
    G.assign(size_t(sHi - sLo) * nr, cd(0.0));
    std::vector<double> J(size_t(nr) * nw);
    for (int i = i0; i < nr; ++i)
      for (int j = 0; j < nw; ++j) J[size_t(i) * nw + j] = std::cyl_bessel_j(0.0, kj[j] * rr[i]);
    std::vector<cd> e(nw);
    for (int n = sLo; n < sHi; ++n) {
      const double d = zv.node(n) - F.zeta_center;
      for (int j = 0; j < nw; ++j) e[j] = cj[j] * std::polar(1.0, wj[j] * d);
      for (int i = i0; i < nr; ++i) {
        cd sum = 0.0;
        const double* Ji = &J[size_t(i) * nw];
        for (int j = 0; j < nw; ++j) sum += Ji[j] * e[j];
        G[size_t(n - sLo) * nr + i] = sum;
      }
    }
  };
  std::vector<double> rf(g.rho.begin(), g.rho.begin() + g.N + 1), rc(g.c.begin(), g.c.begin() + g.N);
  ffN0 = sLo;
  fill(ffGf, rf, g.N + 1, 1);
  fill(ffGc, rc, g.N, 0);
}

bool Simulation::finite() const {
  for (const auto& v : fold.Rr)
    if (!std::isfinite(v.real()) || !std::isfinite(v.imag())) return false;
  for (const auto& v : fold.Ez)
    if (!std::isfinite(v.real()) || !std::isfinite(v.imag())) return false;
  return true;
}

}  // namespace gem
