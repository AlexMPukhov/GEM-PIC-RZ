// ============================================================================
//  simulation.cpp -- initialization (lasers, bunches, their fields), see also simulation_*.cpp
// ============================================================================
#include "simulation.hpp"

#include "eigen.hpp"
#include "rng.hpp"

#include <algorithm>
#if defined(_OPENMP)
#include <omp.h>
#endif
#include <cmath>

namespace gem {

static const cd I1(0.0, 1.0);

// ============================================================================
//  setup
// ============================================================================
Simulation::Simulation(const Config& c, const std::string& restartFile, bool useEngine) : cfg(c) {
  M = cfg.grid.M;
  Nz = cfg.grid.Nz;
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
  fold.init(M, Nz, g.N);
  fnew.init(M, Nz, g.N);
  for (auto* C : {&c0b, &c0s, &c1b, &c1s}) C->init(M, Nz, g.N);
  solver.setup(g, M);
  solver.theta = cfg.numerics.theta;
  sfilter.setup(g, solver, M, cfg.plasma.present ? cfg.plasma.smooth_length : 0.0,
                cfg.plasma.present ? cfg.plasma.smooth_length_m0 : 0.0);
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
  if (restartFile.empty()) {
    initLasers();
    if (cfg.flying_focus.present)  // the incident field of level s0 in the whole box
      for (int n = 0; n < Nz; ++n) {
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
  for (const auto& L : cfg.lasers) {
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
    for (int n = 0; n < Nz; ++n) {
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
      for (int k = 0; k < ncent; ++k) {
        Particle p{};
        p.x = B.sigma_x * nd(rng);
        p.y = B.sigma_y * nd(rng);
        p.zeta = (B.shape == BunchConfig::Shape::Gaussian) ? B.zeta_center + B.sigma_z * nd(rng)
                 : (B.shape == BunchConfig::Shape::Flattop) ? B.zeta_center + B.length * (ud(rng) - 0.5)
                                                            : sampleZeta(ud(rng));
        p.ux = sux * nd(rng);
        p.uy = suy * nd(rng);
        if (B.gamma_fn.empty()) {
          p.uz = uz0 * (1.0 + B.energy_spread * nd(rng));
        } else {  // gamma(zeta): energy chirp; energy_spread is relative to the local value
          const double gz = B.gamma_fn({p.zeta});
          p.uz = std::sqrt(gz * gz - 1.0) * (1.0 + B.energy_spread * nd(rng));
        }
        if (B.z_waist != 0.0) {  // drift back from the waist
          p.x -= p.ux / p.uz * B.z_waist;
          p.y -= p.uy / p.uz * B.z_waist;
        }
        if (B.x0 != 0.0) p.x += B.x0;
        if (B.y0 != 0.0) p.y += B.y0;
        p.q = B.charge / B.macroparticles * M;
        finish(p);
        cent.push_back(p);
      }
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

}  // namespace gem
