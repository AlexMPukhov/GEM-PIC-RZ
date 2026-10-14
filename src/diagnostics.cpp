// ============================================================================
//  diagnostics.cpp -- console output, ASCII field/particle files, bunch table
//
//  Files in output.dir (every output.every steps):
//    axis_NNNNN.txt   zeta[um]  E_x  E_z  n_e/n0  Zn_i/n0            (on axis)
//    map_NNNNN.txt    zeta[um]  x[um]  E_x  E_z  n_e/n0  Zn_i/n0      (x-z plane, gnuplot pm3d)
//    beam_NNNNN.txt   x y zeta[um]  ux uy uz  q[pC]  species
//    data_NNNNN.h5    (output.format = hdf5 or both) all of the above and more:
//                     every azimuthal mode on the native grid (output_hdf5.cpp)
//  Every step:
//    bunches.txt      statistics of every bunch, of the self-trapped plasma
//                     electrons ('trapped') and of the trapped electrons freed
//                     by field ionization ('inj_<ion name>')
//    ionization.txt   per ionizing species: charge-state distribution of the
//                     ions leaving the rear of the window, freed charge
//  Fields are in units of m c omega0 / e; n_e/n0 = -rho_e / (e n0) (m = 0);
//  Zn_i/n0 = charge density of the MOBILE ions / (e n0), 0 without ion blocks.
// ============================================================================
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <string>
#include <vector>

#include "comm.hpp"
#include "simulation.hpp"

namespace gem {

namespace {

struct Stats {
  size_t n = 0;
  double Q = 0;  // |charge| [code]
  double gmean = 1, dgrel = 0, sx = 0, sy = 0, ex = 0, ey = 0, zm = 0, sz = 0;
};

// statistics from the (combined) moments of one species
Stats finish(const Moments& M) {
  Stats s;
  s.n = size_t(M.n);
  s.Q = M.Q;
  if (s.Q <= 0) return s;
  const double Q = s.Q;
  auto var = [&](double a1, double a2) { return std::max(0.0, a2 / Q - (a1 / Q) * (a1 / Q)); };
  auto cov = [&](double a1, double b1, double ab) { return ab / Q - (a1 / Q) * (b1 / Q); };
  s.gmean = M.g1 / Q;
  s.dgrel = std::sqrt(var(M.g1, M.g2)) / s.gmean;
  s.sx = std::sqrt(var(M.x1, M.x2));
  s.sy = std::sqrt(var(M.y1, M.y2));
  s.ex = std::sqrt(std::max(0.0, var(M.x1, M.x2) * var(M.u1, M.u2) - std::pow(cov(M.x1, M.u1, M.xu), 2)));
  s.ey = std::sqrt(std::max(0.0, var(M.y1, M.y2) * var(M.v1, M.v2) - std::pow(cov(M.y1, M.v1, M.yv), 2)));
  s.zm = M.z1 / Q;
  s.sz = std::sqrt(var(M.z1, M.z2));
  return s;
}

Stats statsOf(const DiagSums& D, int species) {
  const auto it = D.moments.find(species);
  return it == D.moments.end() ? Stats{} : finish(it->second);
}

}  // namespace

void Moments::add(const Particle& p) {
  const double w = std::fabs(p.q), gm = gammaOf(p);
  ++n;
  Q += w;
  g1 += w * gm; g2 += w * gm * gm;
  x1 += w * p.x; x2 += w * p.x * p.x; u1 += w * p.ux; u2 += w * p.ux * p.ux; xu += w * p.x * p.ux;
  y1 += w * p.y; y2 += w * p.y * p.y; v1 += w * p.uy; v2 += w * p.uy * p.uy; yv += w * p.y * p.uy;
  z1 += w * p.zeta; z2 += w * p.zeta * p.zeta;
}
void Moments::add(const Moments& o) {
  n += o.n; Q += o.Q; g1 += o.g1; g2 += o.g2; x1 += o.x1; x2 += o.x2; u1 += o.u1; u2 += o.u2; xu += o.xu;
  y1 += o.y1; y2 += o.y2; v1 += o.v1; v2 += o.v2; yv += o.yv; z1 += o.z1; z2 += o.z2;
  br += o.br; bi += o.bi;
}
void DiagSums::add(const DiagSums& o) {
  maxEx = std::max(maxEx, o.maxEx);
  minEz = std::min(minEz, o.minEz);
  maxEz = std::max(maxEz, o.maxEz);
  nimax = std::max(nimax, o.nimax);
  Wrad += o.Wrad;
  W0 += o.W0;
  if (o.Pmax > Pmax) { Pmax = o.Pmax; zPmax = o.zPmax; }
  for (const auto& kv : o.slices) {
    auto& v = slices[kv.first];
    if (v.size() < kv.second.size()) v.resize(kv.second.size(), 0.0);
    for (size_t j = 0; j < kv.second.size(); ++j) v[j] += kv.second[j];
  }
  for (const auto& kv : o.moments) moments[kv.first].add(kv.second);
  for (size_t s = 0; s < o.freed.size(); ++s) {
    if (freed.size() <= s) freed.resize(s + 1, 0.0);
    freed[s] += o.freed[s];
  }
}

// on-axis field extrema over the slices [n0, n1), moments of the beam list
DiagSums Simulation::diagSums(int n0, int n1) const {
  DiagSums D;
  for (int n = n0; n < n1; ++n) {
    EB f = gather(fold, n, 1.0, n, 0.0, 1e-9, 0.0);
    D.maxEx = std::max(D.maxEx, std::fabs(f.Ex));
    D.minEz = std::min(D.minEz, f.Ez);
    D.maxEz = std::max(D.maxEz, f.Ez);
  }
  for (int k = n0; k < n1; ++k)
    D.nimax = std::max(D.nimax, (c0s.rhoi[size_t(k) * g.N] + c0b.rhoi[size_t(k) * g.N]) / n_ref);
  for (const auto& p : beam) D.moments[p.sp].add(p);
  D.freed = freedCharge;
  if (cfg.output.fel) {
    // bunching at the FEL wavelength
    const double k = 2.0 * PI / cfg.output.fel_wavelength;
    const double ls = cfg.output.fel_slice;
    const size_t ns = size_t(std::ceil(zv.node(Nz) / ls)) + 1;
    for (const auto& p : beam) {
      Moments& m = D.moments[p.sp];
      const double w = std::fabs(p.q), c = w * std::cos(k * p.zeta), s = w * std::sin(k * p.zeta);
      m.br += c;
      m.bi += s;
      if (p.sp < 0) continue;
      auto& v = D.slices[p.sp];
      if (v.empty()) v.assign(3 * ns, 0.0);
      const size_t j = std::min(ns - 1, size_t(std::max(0.0, p.zeta / ls)));
      v[3 * j] += w; v[3 * j + 1] += c; v[3 * j + 2] += s;
    }
    // energy and power: |E_r|^2 + |B_phi|^2 = (|R_r|^2 + |L_r|^2) / 2, E_r B_phi* (real part) = (|R_r|^2 - |L_r|^2) / 4
    const int N = g.N;
    for (int n = n0; n < std::min(n1, Nz); ++n) {
      const double hn = zv.len(n);
      double P = 0.0;
      for (int m = 0; m <= M; ++m) {
        const double fac = (m == 0) ? 1.0 : 2.0;
        double e = 0.0, s = 0.0;
        for (int i = 1; i <= N; ++i) {
          const size_t j = fold.f(m, n, i);
          const double R = std::norm(fold.Rr[j]), L = std::norm(fold.Lr[j]);
          e += g.Af[i] * (0.5 * (R + L) + std::norm(fold.Bz[j]));
          s += g.Af[i] * (R - L);
        }
        for (int i = 0; i < N; ++i) {
          const size_t j = fold.c(m, n, i);
          const double R = std::norm(fold.Rp[j]), L = std::norm(fold.Lp[j]);
          e += g.Ac[i] * (0.5 * (R + L) + std::norm(fold.Ez[j]));
          s += g.Ac[i] * (R - L);
        }
        if (m == 0) D.W0 += PI * fac * hn * e;
        else { D.Wrad += PI * fac * hn * e; P += PI * s; }
      }
      if (P > D.Pmax) { D.Pmax = P; D.zPmax = zv.node(n); }
    }
  }
  return D;
}

void Simulation::diagnostics(bool files) { report(diagSums(0, Nz), files); }

// console line, bunches.txt, ionization.txt and (files) the field/particle files
void Simulation::report(const DiagSums& D, bool files) {
  const double um = cfg.units.k0inv_um();
  const Stats tr = statsOf(D, SPECIES_TRAPPED);
  std::printf("step %5d  z = %8.4f mm  |Ex|max = %.3f  Ez = [%.4f, %.4f]  trapped %.3f pC (N = %zu)", step_no,
              s0 * um * 1e-3, D.maxEx, D.minEz, D.maxEz, cfg.units.codeToPC(tr.Q), tr.n);
  if (!cfg.ions.empty()) std::printf("  Zn_i/n0(axis)max = %.3f", D.nimax);
  std::printf("\n");
  for (size_t s = 0; s < cfg.ions.size(); ++s) {
    if (!cfg.ions[s].ionization) continue;
    const auto& v = exitStates[s];
    double tot = 0, zsum = 0;
    for (size_t z = 0; z < v.size(); ++z) { tot += v[z]; zsum += z * v[z]; }
    const Stats inj = statsOf(D, SPECIES_TRAPPED_IONIZED_BASE - int(s));
    const double fr = s < D.freed.size() ? D.freed[s] : 0.0;
    std::printf("      %-10s <Z> leaving the window = %.4f   freed %.4g pC   trapped (ionization injection) %.4f pC"
                " (N = %zu)\n",
                cfg.ions[s].name.c_str(), tot > 0 ? zsum / tot : 0.0, cfg.units.codeToPC(fr),
                cfg.units.codeToPC(inj.Q), inj.n);
  }
  for (size_t b = 0; b < cfg.bunches.size(); ++b) {
    const Stats s = statsOf(D, int(b));
    std::printf("      %-10s N = %6zu  Q = %7.3f pC  E = %9.2f MeV  dE/E = %.2e  sig_x = %.3f um  eps_nx = %.3f um"
                "  sig_z = %.2f um\n",
                cfg.bunches[b].name.c_str(), s.n, cfg.units.codeToPC(s.Q), (s.gmean - 1.0) * MC2_MEV, s.dgrel,
                s.sx * um, s.ex * um, s.sz * um);
  }
  std::fflush(stdout);
  writeBunchTable(D);
  writeIonizationTable(D);
  if (cfg.output.fel) writeFelTable(D);
  if (files && cfg.output.text) writeFiles();
  if (files && cfg.output.hdf5) writeHdf5();
  if (files && cfg.output.openpmd) writeOpenPMD();
}

// After a restart from step N the tables may hold lines of steps > N from the
// interrupted run: keep the comments and the lines of steps <= N.
void Simulation::trimTables() {
  if (comm::rank() != 0) return;  // only rank 0 writes the output
  for (const char* t : {"/bunches.txt", "/ionization.txt", "/fel.txt"}) {
    const std::string name = cfg.output.dir + t;
    FILE* fp = std::fopen(name.c_str(), "r");
    if (!fp) continue;
    std::string keep;
    char line[4096];
    while (std::fgets(line, sizeof line, fp)) {
      int st = 0;
      if (line[0] == '#' || (std::sscanf(line, "%d", &st) == 1 && st <= step_no)) keep += line;
    }
    std::fclose(fp);
    if (FILE* out = std::fopen(name.c_str(), "w")) {
      std::fputs(keep.c_str(), out);
      std::fclose(out);
    }
  }
}

// fel.txt: radiation energy and peak power (modes m >= 1), bunching factor of every bunch
void Simulation::writeFelTable(const DiagSums& D) {
  const std::string name = cfg.output.dir + "/fel.txt";
  const bool fresh = (step_no == 0) || !std::filesystem::exists(name);
  FILE* fp = std::fopen(name.c_str(), fresh ? "w" : "a");
  if (!fp) return;
  const Units& u = cfg.units;
  const double L = u.k0inv_um() * 1e-6, E0 = MC2_MEV * 1e6 / L, eps0 = 8.8541878128e-12;
  const double WJ = eps0 * E0 * E0 * L * L * L, PW = eps0 * C_LIGHT * E0 * E0 * L * L;
  if (fresh) {
    std::fprintf(fp, "# FEL diagnostics at lambda = %.6g um (bunching b = |<exp(i 2 pi zeta / lambda)>|, charge-weighted)\n",
                 cfg.output.fel_wavelength * u.k0inv_um());
    std::fprintf(fp, "# b_max: largest bunching over slices of %.4g um holding >= half the charge of the fullest slice\n",
                 cfg.output.fel_slice * u.k0inv_um());
    std::fprintf(fp, "# step  z[mm]  W_rad[J](m>=1)  P_max[W](m>=1)  zeta(P_max)[um]  W_m0[J]");
    for (const auto& B : cfg.bunches) std::fprintf(fp, "  b_%s  b_max_%s", B.name.c_str(), B.name.c_str());
    std::fprintf(fp, "\n");
  }
  std::fprintf(fp, "%6d %10.5f %13.5e %13.5e %10.4f %13.5e", step_no, s0 * u.k0inv_um() * 1e-3, D.Wrad * WJ,
               D.Pmax * PW, D.zPmax * u.k0inv_um(), D.W0 * WJ);
  for (size_t b = 0; b < cfg.bunches.size(); ++b) {
    const auto it = D.moments.find(int(b));
    const double bf = (it == D.moments.end() || it->second.Q <= 0)
                          ? 0.0 : std::hypot(it->second.br, it->second.bi) / it->second.Q;
    double bmax = 0.0;
    const auto sl = D.slices.find(int(b));
    if (sl != D.slices.end()) {
      const auto& v = sl->second;
      double qmax = 0.0;
      for (size_t j = 0; j < v.size(); j += 3) qmax = std::max(qmax, v[j]);
      for (size_t j = 0; j < v.size(); j += 3)
        if (v[j] > 0 && v[j] >= 0.5 * qmax) bmax = std::max(bmax, std::hypot(v[j + 1], v[j + 2]) / v[j]);
    }
    std::fprintf(fp, " %12.5e %12.5e", bf, bmax);
  }
  std::fprintf(fp, "\n");
  std::fclose(fp);
}

void Simulation::writeBunchTable(const DiagSums& D) {
  namespace fs = std::filesystem;
  fs::create_directories(cfg.output.dir);
  const std::string name = cfg.output.dir + "/bunches.txt";
  const bool fresh = (step_no == 0) || !std::filesystem::exists(name);
  FILE* fp = std::fopen(name.c_str(), fresh ? "w" : "a");
  if (!fp) return;
  const double um = cfg.units.k0inv_um();
  if (fresh)
    std::fprintf(fp,
                 "# step  z[mm]  bunch  N  Q[pC]  E[MeV]  dE/E(rms)  sig_x[um]  sig_y[um]  eps_nx[um]  eps_ny[um]"
                 "  <zeta>[um]  sig_z[um]    (bunch 'trapped' = self-injected plasma electrons)\n");
  auto row = [&](std::string nm, const Stats& s) {
    std::replace(nm.begin(), nm.end(), ' ', '_');  // keep one column per quantity
    std::fprintf(fp, "%6d %10.5f %-12s %8zu %10.4f %12.4f %10.3e %9.4f %9.4f %9.4f %9.4f %10.3f %9.4f\n", step_no,
                 s0 * um * 1e-3, nm.c_str(), s.n, cfg.units.codeToPC(s.Q), (s.gmean - 1.0) * MC2_MEV, s.dgrel, s.sx * um,
                 s.sy * um, s.ex * um, s.ey * um, s.zm * um, s.sz * um);
  };
  for (size_t b = 0; b < cfg.bunches.size(); ++b) row(cfg.bunches[b].name, statsOf(D, int(b)));
  row("trapped", statsOf(D, SPECIES_TRAPPED));
  for (size_t s = 0; s < cfg.ions.size(); ++s)
    if (cfg.ions[s].ionization)
      row("inj_" + cfg.ions[s].name, statsOf(D, SPECIES_TRAPPED_IONIZED_BASE - int(s)));
  std::fclose(fp);
}

void Simulation::writeIonizationTable(const DiagSums& D) {
  bool any = false;
  for (const auto& I : cfg.ions) any = any || I.ionization;
  if (!any) return;
  const std::string name = cfg.output.dir + "/ionization.txt";
  const bool fresh = (step_no == 0) || !std::filesystem::exists(name);
  FILE* fp = std::fopen(name.c_str(), fresh ? "w" : "a");
  if (!fp) return;
  const double um = cfg.units.k0inv_um();
  if (fresh)
    std::fprintf(fp, "# step  z[mm]  species  <Z>  freed[pC]  f(Z=0) f(Z=1) ... f(Z=Zmax)\n"
                     "#   f(Z) = fraction of the atoms with charge state Z among the ions leaving the rear of\n"
                     "#   the window in this step (all radii); freed = charge of the electrons freed in this step\n");
  for (size_t s = 0; s < cfg.ions.size(); ++s) {
    const auto& I = cfg.ions[s];
    if (!I.ionization) continue;
    std::string nm = I.name;
    std::replace(nm.begin(), nm.end(), ' ', '_');
    const auto& v = exitStates[s];
    double tot = 0, zsum = 0;
    for (size_t z = 0; z < v.size(); ++z) { tot += v[z]; zsum += z * v[z]; }
    std::fprintf(fp, "%6d %10.5f %-10s %8.5f %11.4e", step_no, s0 * um * 1e-3, nm.c_str(), tot > 0 ? zsum / tot : 0.0,
                 cfg.units.codeToPC(s < D.freed.size() ? D.freed[s] : 0.0));
    for (size_t z = 0; z < v.size(); ++z) std::fprintf(fp, " %9.6f", tot > 0 ? v[z] / tot : 0.0);
    std::fprintf(fp, "\n");
  }
  std::fclose(fp);
}

void Simulation::writeFiles() {
  namespace fs = std::filesystem;
  const auto& O = cfg.output;
  fs::create_directories(O.dir);
  const double um = cfg.units.k0inv_um();
  char name[1024];
  // densities in half-cell k (the last slice has no half-cell: use k = Nz-1)
  auto ne = [&](int k, int i) {
    k = std::min(k, Nz - 1);
    return -(c0s.rho[size_t(k) * g.N + i] + c0b.rho[size_t(k) * g.N + i]) / n_ref;
  };
  auto ni = [&](int k, int i) {
    k = std::min(k, Nz - 1);
    return (c0s.rhoi[size_t(k) * g.N + i] + c0b.rhoi[size_t(k) * g.N + i]) / n_ref;
  };

  if (O.axis) {
    std::snprintf(name, sizeof name, "%s/axis_%05d.txt", O.dir.c_str(), step_no);
    if (FILE* fp = std::fopen(name, "w")) {
      std::fprintf(fp, "# zeta[um]  Ex  Ez  ne/n0  Zni/n0     (z = %.6g mm, fields in m c omega0/e)\n",
                   s0 * um * 1e-3);
      for (int n = 0; n < Nz; ++n) {
        EB f = gather(fold, n, 1.0, n, 0.0, 1e-9, 0.0);
        std::fprintf(fp, "%.6e %.6e %.6e %.6e %.6e\n", zv.node(n) * um, f.Ex, f.Ez, ne(n, 0), ni(n, 0));
      }
      std::fclose(fp);
    }
  }
  if (O.maps) {
    std::snprintf(name, sizeof name, "%s/map_%05d.txt", O.dir.c_str(), step_no);
    if (FILE* fp = std::fopen(name, "w")) {
      std::fprintf(fp, "# zeta[um]  x[um]  Ex  Ez  ne/n0(m=0)  Zni/n0(m=0)     (z = %.6g mm, plane phi = 0 / pi)\n",
                   s0 * um * 1e-3);
      for (int n = 0; n < Nz; n += O.map_stride) {
        for (int sgn = -1; sgn <= 1; sgn += 2)
          for (int ii = 0; ii < g.N; ++ii) {
            int i = (sgn < 0) ? g.N - 1 - ii : ii;
            double x = sgn * g.c[i];
            EB f = gather(fold, n, 1.0, n, 0.0, x, 0.0);
            std::fprintf(fp, "%.5e %.5e %.5e %.5e %.5e %.5e\n", zv.node(n) * um, x * um, f.Ex, f.Ez, ne(n, i), ni(n, i));
          }
        std::fprintf(fp, "\n");
      }
      std::fclose(fp);
    }
  }
  if (O.particles) {
    std::snprintf(name, sizeof name, "%s/beam_%05d.txt", O.dir.c_str(), step_no);
    if (FILE* fp = std::fopen(name, "w")) {
      std::fprintf(fp, "# x[um] y[um] zeta[um] ux uy uz q[pC] species  (species: bunch index, -1 self-trapped plasma "
                       "electrons, -30-s trapped electrons from ionization of ion block s)\n");
      for (const auto& p : beam)
        std::fprintf(fp, "%.6e %.6e %.6e %.6e %.6e %.6e %.6e %d\n", p.x * um, p.y * um, p.zeta * um, p.ux, p.uy,
                     p.uz, cfg.units.codeToPC(p.q), p.sp);
      std::fclose(fp);
    }
  }
}

}  // namespace gem
