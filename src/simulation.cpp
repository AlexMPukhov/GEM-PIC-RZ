// ============================================================================
//  simulation.cpp -- initialization, particle push, deposition, time step
// ============================================================================
#include "simulation.hpp"

#include "eigen.hpp"
#include "rng.hpp"

#include <algorithm>
#include <cmath>

namespace gem {

static const cd I1(0.0, 1.0);

// ============================================================================
//  setup
// ============================================================================
Simulation::Simulation(const Config& c, const std::string& restartFile) : cfg(c) {
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
  if (restartFile.empty()) {
    initLasers();
    initBunches();
    initBunchFields();
  } else {
    readCheckpoint(restartFile);  // fields, currents, particles, random generator, counters
    trimTables();
  }
  beginStep();  // matrices of the first step (the engine uploads them)
  if (cfg.numerics.backend != NumericsConfig::Backend::Serial) engine = std::make_unique<Engine>(*this);
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
  Currents C;
  C.init(M, Nz, N);
  std::vector<cd> rho(size_t(Nz) * N), psi(size_t(Nz) * N), col(Nz), Er(N + 1), Ep(N);
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
      for (int n = 0; n < Nz; ++n)
        for (int i = 0; i < N; ++i)
          rho[size_t(n) * N + i] = 0.5 * (zv.wLo(n) * ((n > 0) ? C.jz[C.c(m, n - 1, i)] : cd(0.0)) +
                                          zv.wHi(n) * C.jz[C.c(m, n, i)]) / beta;
      if (ig2 == 0.0) {
        for (int n = 0; n < Nz; ++n) solver.poisson(m, &rho[size_t(n) * N], &psi[size_t(n) * N]);
      } else {
        eigen(m);
        const std::vector<double>& Z = Q[m];
        // rho_k(n) = sum_i Q_ik sqrt(A_i) rho_i(n)
        std::vector<cd> rk(size_t(Nz) * N);
        for (int n = 0; n < Nz; ++n) {
          const cd* r = &rho[size_t(n) * N];
          cd* out = &rk[size_t(n) * N];
          for (int i = 0; i < N; ++i) {
            const cd v = sq[i] * r[i];
            if (v == cd(0.0)) continue;
            const double* zi = &Z[size_t(i) * N];
            for (int k = 0; k < N; ++k) out[k] += zi[k] * v;
          }
        }
        // tridiagonal solve along zeta for every eigenmode
        std::vector<double> lo(Nz), di(Nz), up(Nz), cp(Nz);
        for (int n = 0; n < Nz; ++n) {
          const double hm = (n > 0) ? zv.len(n - 1) : zv.len(0), hp = zv.len(n);
          const double f = 2.0 * ig2 / (hm + hp);
          lo[n] = -f / hm; up[n] = -f / hp; di[n] = f / hm + f / hp;
        }
        for (int k = 0; k < N; ++k) {
          const double lk = lam[m][k];
          for (int n = 0; n < Nz; ++n) col[n] = rk[size_t(n) * N + k];
          double den = di[0] + lk;
          cp[0] = up[0] / den;
          col[0] /= den;
          for (int n = 1; n < Nz; ++n) {
            den = di[n] + lk - lo[n] * cp[n - 1];
            cp[n] = up[n] / den;
            col[n] = (col[n] - lo[n] * col[n - 1]) / den;
          }
          for (int n = Nz - 2; n >= 0; --n) col[n] -= cp[n] * col[n + 1];
          for (int n = 0; n < Nz; ++n) rk[size_t(n) * N + k] = col[n];
        }
        // psi_i(n) = sum_k Q_ik psi_k(n) / sqrt(A_i)
        for (int n = 0; n < Nz; ++n) {
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
      for (int n = 0; n < Nz; ++n) {
        solver.gradient(m, &psi[size_t(n) * N], Er.data(), Ep.data());
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
            const cd up1 = (n + 1 < Nz) ? psi[size_t(n + 1) * N + i] : cd(0.0);
            const cd dn1 = (n > 0) ? psi[size_t(n - 1) * N + i] : cd(0.0);
            fold.Ez[fold.c(m, n, i)] -= sign * (ig2 * (up1 - dn1) / (hm + hp));
          }
        }
      }
    }
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
                         double W, bool ion) {
  if (k < 0 || k >= Nz || wk == 0.0 || W == 0.0) return;
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
      C.jp[C.c(m, k, i)] += f * sp * vp * em;
      C.jz[C.c(m, k, i)] += f * sz * vz * em;
    };
    if (wc.mode == 0) { addC(wc.a, wc.wa, wc.wa); addC(wc.a + 1, wc.wb, wc.wb); }
    else if (wc.mode == 2) addC(N - 1, 1.0, 1.0);
    else addC(0, 1.0, (m == 0) ? 1.0 : std::pow(r / g.c[0], m));
    // faces: j_r (face 0 is not used by the solver)
    auto addF = [&](int i, double s) {
      if (i == 0) return;
      C.jr[C.f(m, k, i)] += norm * s / g.Vf[i] * vr * em;
    };
    if (wf.mode == 0) { addF(wf.a, wf.wa); addF(wf.a + 1, wf.wb); }
    else addF(N, 1.0);
    em *= emi;
  }
  // m = 0 charge density (diagnostic), separately for ions
  std::vector<double>& dens = ion ? C.rhoi : C.rho;
  auto addR = [&](int i, double s) { dens[size_t(k) * N + i] += norm * s / g.Vc[i]; };
  if (wc.mode == 0) { addR(wc.a, wc.wa); addR(wc.a + 1, wc.wb); }
  else if (wc.mode == 2) addR(N - 1, 1.0);
  else addR(0, 1.0);
}

// Numerical Cherenkov: the dispersion relation of the scheme has omega = 0 at
// the zeta-Nyquist wavenumber (kappa h = pi) for every transverse k, so this
// odd-even mode is resonant with every v ~ c source. Beam particles keep their
// zeta from step to step, so the Nyquist part of their current (shot noise)
// drives it secularly. Passes of the (1,2,1)/4 filter remove that part
// (response cos^2(kappa h/2)); zero beyond the ends of the window.
void Simulation::filterBeamCurrent(Currents& C) const {
  const int passes = cfg.numerics.beam_filter;
  if (passes <= 0) return;
  auto column = [&](std::vector<cd>& v, int L) {
    for (int m = 0; m <= M; ++m)
      for (int i = 0; i < L; ++i)
        for (int p = 0; p < passes; ++p) {
          cd prev(0.0);
          for (int k = 0; k < Nz; ++k) {
            cd& x = v[(size_t(m) * Nz + k) * L + i];
            const cd cur = x;
            const cd next = (k + 1 < Nz) ? v[(size_t(m) * Nz + k + 1) * L + i] : cd(0.0);
            x = 0.25 * prev + 0.5 * cur + 0.25 * next;
            prev = cur;
          }
        }
  };
  column(C.jr, g.N + 1);
  column(C.jp, g.N);
  column(C.jz, g.N);
}

// beam-type deposit (current at level s^1): W = q / v_z, linear in zeta
// between half-cell centres
void Simulation::depositBeamParticle(Currents& C, const Particle& p) {
  double gm = gammaOf(p);
  double vx = p.ux / gm, vy = p.uy / gm, vz = p.uz / gm;
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
  const int nsub = std::max(1, int(std::ceil(Delta / cfg.numerics.beam_dz_max)));
  const double dz = Delta / nsub;
  const double rmax = g.rho[g.N];
  for (auto& p : beam) {
    for (int sub = 0; sub < nsub && p.status == 0; ++sub) {
      double gm = gammaOf(p);
      double vz = p.uz / gm;
      if (vz < cfg.numerics.vz_min) { p.status = 3; ++n_restreamed; break; }
      double dt = dz / vz;
      double w;
      const int na = zv.cell(p.zeta, w);
      if (na < 0 || na >= Nz) { p.status = 2; break; }
      EB f = gather(fold, na, 1.0 - w, na + 1, w, p.x, p.y);
      addExternal(ext, p.x, p.y, p.z, f);
      if (cfg.numerics.beam_vay) vay(p.ux, p.uy, p.uz, f, p.qm, dt);
      else boris(p.ux, p.uy, p.uz, f, p.qm, dt);
      gm = gammaOf(p);
      p.x += p.ux / gm * dt; p.y += p.uy / gm * dt;
      p.zeta += (p.uz / gm - 1.0) * dt;
      p.z += dz;
      if (p.zeta < 0.0 || std::hypot(p.x, p.y) >= rmax) p.status = 2;
    }
  }
  for (auto& p : beam)
    if (p.status == 3) {
      Particle q = p;
      q.status = 0;
      pending.push_back(q);
    }
  beam.erase(std::remove_if(beam.begin(), beam.end(), [](const Particle& p) { return p.status != 0; }), beam.end());
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
void Simulation::pushStreamParticle(Particle& p, int n, double s1) {
  const bool ion = isIon(p.sp);
  const IonConfig* I = (ion && p.sp != SPECIES_TRAPPED_ION) ? &cfg.ions[ionIndex(p.sp)] : nullptr;
  // length of the path in zeta within this half-cell
  const double h = zv.len(n), zn = zv.node(n);
  const double hl = (p.zeta >= zv.node(n + 1) - 1e-9 * h) ? h : std::max(0.0, p.zeta - zn);
  const double wm = 0.5 * hl / h;  // weight of slice n+1 at the midpoint of the path
  if (I && !I->mobile) {
    // fixed ion: no motion (dt = h), fields only needed for ionization
    deposit(c1s, n, 1.0, p.x, p.y, 0.0, 0.0, 0.0, p.q * hl / Delta, true);  // density diagnostic only
    const double zs = p.zeta;
    p.zeta = zn;
    if (I->ionization && p.Zi < I->Zmax && sliceFieldBound >= adk[ionIndex(p.sp)].minField()) {
      const double w1 = (zs - zn) / h;  // field at the start of the path
      const EB f1 = gather(fnew, n + 1, w1, n, 1.0 - w1, p.x, p.y);
      const EB f0 = gather(fnew, n, 1.0, n, 0.0, p.x, p.y);
      ionize(p, std::hypot(f1.Ex, f1.Ey, f1.Ez), std::hypot(f0.Ex, f0.Ey, f0.Ez), hl);
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
  if (ext.any()) addExternal(ext, xm, ym, p.z + 0.5 * dt0 * p.uz / g0, f);
  double ux = p.ux, uy = p.uy, uz = p.uz;
  boris(ux, uy, uz, f, p.qm, dt0);
  double g1 = std::sqrt(1.0 + ux * ux + uy * uy + uz * uz);
  double dt = 2.0 * hl / (omv0 + (g1 - uz) / g1);  // time for |dzeta| = hl
  ux = p.ux; uy = p.uy; uz = p.uz;
  boris(ux, uy, uz, f, p.qm, dt);
  g1 = std::sqrt(1.0 + ux * ux + uy * uy + uz * uz);
  double vx = 0.5 * (p.ux / g0 + ux / g1), vy = 0.5 * (p.uy / g0 + uy / g1), vz = 0.5 * (p.uz / g0 + uz / g1);
  double zeta_end = zn;
  if (p.z + vz * dt >= s1) {  // reaches s^1 inside this half-cell -> trapped
    double dtp = std::max(0.0, (s1 - p.z) / std::max(vz, 1e-3));
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
  deposit(c1s, n, 1.0, 0.5 * (xo + p.x), 0.5 * (yo + p.y), vx, vy, vz, p.q * dt / Delta, ion);
  if (std::hypot(p.x, p.y) >= g.rho[g.N]) p.status = 2;
  // field ionization with the mid-step field (after the deposit: the new
  // charge state acts from the next half-cell on)
  if (ionizing && p.status == 0)
    ionize(p, std::hypot(f1.Ex, f1.Ey, f1.Ez), std::hypot(f0.Ex, f0.Ey, f0.Ez), dt);
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
void Simulation::ionize(Particle& p, double E1, double E0, double dt) {
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
    born.push_back(e);
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
  const int N = g.N;
  const double s1 = s0 + Delta;
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
      // field solve, all modes
      for (int m = 0; m <= M; ++m) {
        // four-point current sums, weighted by the interpolation to slice n
        // (wl = wh = 1 on a uniform grid); at the last slice (n = 0) the missing
        // half-cell -1 is replaced by half-cell 0, at s^1 with the current of
        // this pass (a zero there would halve the source of R(0): a growing spike)
        const int nm = (n > 0) ? n - 1 : 0;
        const double wl = zv.wLo(n), wh = zv.wHi(n);
        for (int i = 0; i <= N; ++i) {
          const cd a = c0b.jr[c0b.f(m, nm, i)] + c0s.jr[c0s.f(m, nm, i)];
          const cd pr = (n > 0) ? predR[m * (N + 1) + i] : c1b.jr[c1b.f(m, 0, i)] + c1s.jr[c1s.f(m, 0, i)];
          Jr[i] = wl * a + wh * c0b.jr[c0b.f(m, n, i)] + wh * c0s.jr[c0s.f(m, n, i)] + wl * pr +
                  wh * c1b.jr[c1b.f(m, n, i)] + wh * c1s.jr[c1s.f(m, n, i)];
        }
        for (int i = 0; i < N; ++i) {
          const cd a = c0b.jp[c0b.c(m, nm, i)] + c0s.jp[c0s.c(m, nm, i)];
          const cd pr = (n > 0) ? predP[m * N + i] : c1b.jp[c1b.c(m, 0, i)] + c1s.jp[c1s.c(m, 0, i)];
          Jp[i] = wl * a + wh * c0b.jp[c0b.c(m, n, i)] + wh * c0s.jp[c0s.c(m, n, i)] + wl * pr +
                  wh * c1b.jp[c1b.c(m, n, i)] + wh * c1s.jp[c1s.c(m, n, i)];
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
      const double q = 0.25 * Delta * zv.wLo(n + 1);
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
  // charge states of the ions leaving the rear of the window
  for (auto& v : exitStates) std::fill(v.begin(), v.end(), 0.0);
  for (const auto& p : stream)
    if (p.status == 0 && isIon(p.sp) && p.sp != SPECIES_TRAPPED_ION) exitStates[ionIndex(p.sp)][p.Zi] += p.w;
  absorber();
  std::swap(fold, fnew);
  std::swap(c0b, c1b);
  std::swap(c0s, c1s);
  s0 = s1;
  ++step_no;
  if (!cfg.chicanes.empty()) applyChicanes(za, s0);
}

// damp R near the outer radial boundary
void Simulation::absorber() {
  const double ra = cfg.grid.r_absorber, rm = g.rho[g.N];
  if (ra >= rm) return;
  auto damp = [&](double r) {
    if (r <= ra) return 1.0;
    double x = (r - ra) / (rm - ra);
    return std::max(0.0, 1.0 - cfg.grid.absorber_strength * x * x);
  };
  for (int m = 0; m <= M; ++m)
    for (int n = 0; n < Nz; ++n) {
      for (int i = 0; i <= g.N; ++i) fnew.Rr[fnew.f(m, n, i)] *= damp(g.rho[i]);
      for (int i = 0; i < g.N; ++i) fnew.Rp[fnew.c(m, n, i)] *= damp(g.c[i]);
    }
}

bool Simulation::finite() const {
  for (const auto& v : fold.Rr)
    if (!std::isfinite(v.real()) || !std::isfinite(v.imag())) return false;
  for (const auto& v : fold.Ez)
    if (!std::isfinite(v.real()) || !std::isfinite(v.imag())) return false;
  return true;
}

}  // namespace gem
