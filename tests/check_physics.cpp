// ============================================================================
//  check_physics.cpp -- checks of the physics tests against analytic results
//
//  check_physics bunchfield  <map_00001.txt>   (tests/bunch_vacuum.deck)
//  check_physics ioncollapse <axis_00001.txt>  (tests/ion_collapse.deck)
//  check_physics ionH        <ionization.txt>  (tests/ionization_H.deck)
//  check_physics ionN        <ionization.txt>  (tests/ionization_N.deck)
//  check_physics ffbeat      <map> ...         (tests/flying_focus_beat.sh)
//  The ionization references integrate the ADK rate (adkRateSI, NIST energies)
//  over the analytic laser field E = a0 exp(-r^2/w0^2) exp(-zeta^2/L^2) cos(zeta)
//  (single level: 1 - exp(-int w dt); two levels: exact rate equations per
//  step) and average over the loaded disk r < 50 k0^-1.
//  Exit code 0 on success.
// ============================================================================
#include <cmath>
#include <complex>
#include <cstdlib>
#include <map>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "atomic_data.hpp"
#include "ionization.hpp"

static std::vector<std::vector<double>> read(const char* path) {
  std::vector<std::vector<double>> rows;
  std::ifstream in(path);
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#') continue;
    std::istringstream ss(line);
    std::vector<double> r;
    double x;
    while (ss >> x) r.push_back(x);
    if (!r.empty()) rows.push_back(r);
  }
  return rows;
}

// charge-state fractions after the laser pulse, averaged over the disk r < R
static std::vector<double> ionReference(const char* sym, int Z0, int nlev, double a0) {
  const gem::Element* el = gem::findElement(sym);
  const double lam = 0.8e-6, E0 = 2 * M_PI * 0.51099895e6 / lam, w0 = 2 * M_PI * 299792458.0 / lam;
  const double L = 25.0, W = 30.0, R = 50.0, h = 0.01;
  const int nr = 500;
  std::vector<double> avg(nlev + 1, 0.0);
  double norm = 0.0;
  for (int ir = 0; ir <= nr; ++ir) {
    const double r = R * ir / nr, wr = (ir == 0 || ir == nr ? 0.5 : 1.0) * r;  // trapezoid in r, weight 2 r dr
    std::vector<double> N(nlev + 1, 0.0);
    N[0] = 1.0;
    for (double z = -90.0; z < 90.0; z += h) {
      const double E = a0 * std::exp(-r * r / (W * W)) * std::exp(-z * z / (L * L)) * std::fabs(std::cos(z)) * E0;
      if (nlev == 1) {
        N[0] *= std::exp(-gem::adkRateSI(el->ip_eV[Z0], Z0 + 1, E) * h / w0);
        N[1] = 1.0 - N[0];
      } else {  // two-level chain, exact for constant rates over the step
        const double a = gem::adkRateSI(el->ip_eV[Z0], Z0 + 1, E) * h / w0;
        const double b = gem::adkRateSI(el->ip_eV[Z0 + 1], Z0 + 2, E) * h / w0;
        if (a == 0 && b == 0) continue;
        const double e5 = std::exp(-a), e6 = std::exp(-b);
        const double n6 = N[1] * e6 + (std::fabs(b - a) > 1e-12 ? N[0] * a * (e5 - e6) / (b - a) : N[0] * a * e5);
        const double n5 = N[0] * e5;
        N[2] += (N[0] - n5) + (N[1] - n6);
        N[0] = n5;
        N[1] = n6;
      }
    }
    for (int k = 0; k <= nlev; ++k) avg[k] += wr * N[k];
    norm += wr;
  }
  for (auto& v : avg) v /= norm;
  return avg;
}

int main(int argc, char** argv) {
  if (argc < 3) { std::fprintf(stderr, "usage: check_physics bunchfield|ioncollapse <file>\n"); return 2; }
  const std::string mode = argv[1];
  const double k0inv_um = 0.8 / (2.0 * M_PI);
  auto rows = read(argv[2]);
  if (rows.empty() && mode.rfind("ion", 0) != 0) { std::fprintf(stderr, "cannot read %s\n", argv[2]); return 2; }
  int fails = 0;

  if (mode == "bunchfield") {
    // map columns: zeta[um] x[um] Ex Ez ne ni ; bunch interior 1200 < zeta < 4400 k0^-1
    const double nb = 1e-3, sig = 20.0;
    const double radii[] = {56.0, 145.0, 294.0};
    for (double R : radii) {
      double sum = 0, th = 0, bestdx = 1e9, xr = 0;
      int cnt = 0;
      for (const auto& r : rows) {  // nearest radius on the grid
        const double x = r[1] / k0inv_um;
        if (x > 0 && std::fabs(x - R) < bestdx) { bestdx = std::fabs(x - R); xr = x; }
      }
      th = -(nb * sig * sig / xr) * (1.0 - std::exp(-xr * xr / (2 * sig * sig)));
      double ez2 = 0;
      for (const auto& r : rows) {
        const double z = r[0] / k0inv_um, x = r[1] / k0inv_um;
        if (z > 1200 && z < 4400 && std::fabs(x - xr) < 1e-9 * xr + 1e-6) { sum += r[2]; ez2 += r[3] * r[3]; ++cnt; }
      }
      const double ratio = sum / cnt / th, ezrms = std::sqrt(ez2 / cnt);
      const bool ok = std::fabs(ratio - 1.0) < 0.005 && ezrms < 0.05 * std::fabs(th);
      std::printf("%s bunch field at r = %.1f k0^-1: <E_r>/Gauss = %.4f (tolerance 0.5%%), rms E_z/E_r = %.3f\n",
                  ok ? "OK  " : "FAIL", xr, ratio, ezrms / std::fabs(th));
      if (!ok) ++fails;
    }
  } else if (mode == "ioncollapse") {
    // axis columns: zeta[um] Ex Ez ne ni ; head of the bunch at zeta = 4800 k0^-1.
    // test-ion model values (analytic field, full blowout, no ion self-field)
    const double dist[] = {1800, 2000, 2200, 2400, 2600};
    const double model[] = {2.88, 3.89, 5.76, 9.79, 20.59};
    const double tol[] = {0.05, 0.05, 0.05, 0.07, 0.10};
    for (int j = 0; j < 5; ++j) {
      double best = 1e9, ni = 0;
      for (const auto& r : rows) {
        const double behind = 4800.0 - r[0] / k0inv_um;
        if (std::fabs(behind - dist[j]) < best) { best = std::fabs(behind - dist[j]); ni = r[4]; }
      }
      const bool ok = std::fabs(ni / model[j] - 1.0) < tol[j];
      std::printf("%s ion density on axis %4.0f k0^-1 behind the head: %.2f (model %.2f, tolerance %.0f%%)\n",
                  ok ? "OK  " : "FAIL", dist[j], ni, model[j], 100 * tol[j]);
      if (!ok) ++fails;
    }
  } else if (mode == "ionH" || mode == "ionN") {
    // last row of ionization.txt: step z name <Z> freed f(0) ... f(Zmax)
    std::ifstream in(argv[2]);
    std::string line, last;
    while (std::getline(in, line))
      if (!line.empty() && line[0] != '#') last = line;
    std::istringstream ss(last);
    std::vector<std::string> tok;
    std::string t;
    while (ss >> t) tok.push_back(t);
    auto f = [&](int Z) { return std::stod(tok.at(5 + Z)); };
    if (mode == "ionH") {
      const auto ref = ionReference("H", 0, 1, 0.01);
      const double sim = f(1), tol = 0.05;
      const bool ok = std::fabs(sim / ref[1] - 1.0) < tol;
      std::printf("%s H ionized fraction (disk r < 50): %.5f, ADK reference %.5f (tolerance %.0f%%)\n",
                  ok ? "OK  " : "FAIL", sim, ref[1], 100 * tol);
      if (!ok) ++fails;
    } else {
      const auto ref = ionReference("N", 5, 2, 2.2);
      const double tols[] = {0.005, 0.10, 0.35};
      for (int k = 0; k < 3; ++k) {
        const double sim = f(5 + k);
        const bool ok = std::fabs(sim / ref[k] - 1.0) < tols[k];
        std::printf("%s N%d+ fraction (disk r < 50): %.5f, rate equations %.5f (tolerance %.1f%%)\n",
                    ok ? "OK  " : "FAIL", 5 + k, sim, ref[k], 100 * tols[k]);
        if (!ok) ++fails;
      }
    }
  } else if (mode == "ffbeat") {
    // flying focus of several luminal Bessel components (tests/flying_focus_beat.sh), map at y = 0:
    //   ffbeat <map> <tol_all> <tol_odd> <tol_peak[um]> <ds[um]> <a0> <L_1e[um]> <zeta_c[um]>
    //          then per component <theta[rad]> <l> <amplitude> <phase>
    // E_x(x, 0, zeta, s) = Re sum_c U_c e^{i l_c phi} (linear and circular polarization alike; phi = 0, pi), with
    // U_c = a0 A_c e^{i psi_c} int dw a(w) J_|l|(k_c(w) r) e^{i w (zeta - zeta_c) - i kappa_c s},
    // a(w) = L/(2 sqrt(pi)) e^{-(w-1)^2 L^2/4}, kappa_c = 1 - cos theta_c, k_c^2 = 2 w kappa_c - kappa_c^2 (vacuum),
    // integrated here on 401 frequencies over +-8/L, at the level s = z + ds. Compared within |zeta - zeta_c| < 2 L,
    // r < 40 um: the whole field (rms difference / a0 < tol_all) and its part odd in x (the components of odd l:
    // rms difference / its rms < tol_odd); and the x of the intensity maximum (E_x^2 summed over |zeta - zeta_c| < L/2)
    // near the axis, which the beat displaces by x_b cos(...) (difference < tol_peak).
    if (argc < 15 || (argc - 11) % 4 != 0) { std::fprintf(stderr, "ffbeat: wrong arguments\n"); return 2; }
    const double tol = std::atof(argv[3]), tolOdd = std::atof(argv[4]), tolPeak = std::atof(argv[5]),
                 ds = std::atof(argv[6]) / k0inv_um, a0 = std::atof(argv[7]), L = std::atof(argv[8]) / k0inv_um,
                 zc = std::atof(argv[9]) / k0inv_um;
    struct C { double th; int l; double A, psi; };
    std::vector<C> cs;
    for (int a = 10; a + 3 < argc; a += 4)
      cs.push_back({std::atof(argv[a]), std::atoi(argv[a + 1]), std::atof(argv[a + 2]), std::atof(argv[a + 3])});
    double s = 0;
    {
      std::ifstream in(argv[2]);
      std::string head;
      std::getline(in, head);
      const auto p = head.find("z = ");
      if (p == std::string::npos) { std::fprintf(stderr, "no z in the map header\n"); return 2; }
      s = std::atof(head.c_str() + p + 4) * 1e3 / k0inv_um + ds;
    }
    const int nw = 401;
    auto ref = [&](double x, double zeta) {
      const double r = std::fabs(x);
      double ex = 0.0;
      for (const auto& c : cs) {
        const double kap = 1.0 - std::cos(c.th), dw = 16.0 / L / (nw - 1);
        std::complex<double> U = 0.0;
        for (int j = 0; j < nw; ++j) {
          const double w = 1.0 - 8.0 / L + j * dw;
          const double a = L / (2.0 * std::sqrt(M_PI)) * std::exp(-0.25 * (w - 1) * (w - 1) * L * L) * dw;
          const double kr = std::sqrt(std::max(0.0, 2.0 * w * kap - kap * kap));
          U += a * std::cyl_bessel_j(double(std::abs(c.l)), kr * r) * std::polar(1.0, w * (zeta - zc) - kap * s);
        }
        U *= a0 * c.A * std::polar(1.0, c.psi);
        ex += (U * std::polar(1.0, x < 0 ? c.l * M_PI : 0.0)).real();
      }
      return ex;
    };
    // rows: zeta x Ex ...; pair x with -x (the map has both signs on the same centres)
    double d2 = 0, n2 = 0, o2 = 0, od2 = 0;
    long cnt = 0;
    std::map<std::pair<long, long>, double> E;  // (zeta, x) in 1e-4 um -> Ex
    for (const auto& r : rows) E[{std::lround(r[0] * 1e4), std::lround(r[1] * 1e4)}] = r[2];
    std::map<long, double> Ic, Ir;  // intensity vs x (code, reference)
    for (const auto& r : rows) {
      const double zeta = r[0] / k0inv_um, x = r[1] / k0inv_um;
      if (std::fabs(zeta - zc) > 2 * L || std::fabs(r[1]) > 40.0) continue;
      const double e = ref(x, zeta);
      d2 += (r[2] - e) * (r[2] - e);
      n2 += e * e;
      ++cnt;
      if (std::fabs(zeta - zc) < 0.5 * L && std::fabs(r[1]) < 3.0) {
        Ic[std::lround(r[1] * 1e4)] += r[2] * r[2];
        Ir[std::lround(r[1] * 1e4)] += e * e;
      }
      if (r[1] > 0) {  // odd part: (E(x) - E(-x)) / 2
        const auto it = E.find({std::lround(r[0] * 1e4), std::lround(-r[1] * 1e4)});
        if (it == E.end()) continue;
        const double oc = 0.5 * (r[2] - it->second), orf = 0.5 * (e - ref(-x, zeta));
        o2 += orf * orf;
        od2 += (oc - orf) * (oc - orf);
      }
    }
    if (cnt == 0) { std::fprintf(stderr, "no map points in the window\n"); return 2; }
    // x of the maximum: parabola through the largest node and its neighbours
    auto peak = [](const std::map<long, double>& I) {
      auto best = I.begin();
      for (auto it = I.begin(); it != I.end(); ++it)
        if (it->second > best->second) best = it;
      if (best == I.begin() || std::next(best) == I.end()) return best->first * 1e-4;
      const auto lo = std::prev(best), hi = std::next(best);
      const double x0 = lo->first * 1e-4, x1 = best->first * 1e-4, x2 = hi->first * 1e-4;
      const double y0 = lo->second, y1 = best->second, y2 = hi->second;
      const double den = (x0 - x1) * (x0 - x2) * (x1 - x2);
      const double A = (x2 * (y1 - y0) + x1 * (y0 - y2) + x0 * (y2 - y1)) / den;
      const double B = (x2 * x2 * (y0 - y1) + x1 * x1 * (y2 - y0) + x0 * x0 * (y1 - y2)) / den;
      return -B / (2 * A);
    };
    const double rel = std::sqrt(d2 / cnt) / a0, relodd = std::sqrt(od2 / std::max(o2, 1e-300));
    const double pc = peak(Ic), pr = peak(Ir);
    const bool ok = rel < tol && relodd < tolOdd && std::fabs(pc - pr) < tolPeak;
    std::printf("%s z = %.4f mm (reference at s = %.4f mm): rms(E_x - reference)/a0 = %.2e (tolerance %.1e), "
                "rms E_x/a0 = %.3f; odd part (beat components): rms difference / rms = %.2e (tolerance %.2g), "
                "rms/a0 = %.3e; intensity maximum at x = %.4f um (reference %.4f um, tolerance %.2g um)\n",
                ok ? "OK  " : "FAIL", (s - ds) * k0inv_um * 1e-3, s * k0inv_um * 1e-3, rel, tol,
                std::sqrt(n2 / cnt) / a0, relodd, tolOdd, std::sqrt(o2 / std::max(1L, cnt / 2)) / a0, pc, pr, tolPeak);
    if (!ok) ++fails;
  } else {
    std::fprintf(stderr, "unknown mode %s\n", mode.c_str());
    return 2;
  }
  return fails ? 1 : 0;
}
