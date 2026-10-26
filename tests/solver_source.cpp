// ============================================================================
//  solver_source.cpp -- reduced test of the field solver with a prescribed source
//
//  FieldSolver::solveSlice is driven directly (no particles) by a current
//      j(r, zeta, s) = f(r) exp(i (k zeta + K s))
//  in one azimuthal mode m on a uniform radial grid, marched over s with the
//  four-point current sums of Simulation::march (trapezoid, theta = 1/2). The
//  profile E_r(r) of the middle slice is compared with independent references:
//
//   1. grid convergence: a smooth m = 2 source (j_z, or a transverse current
//      j_x = J x^ of an offset Gaussian beam) at dr = 1 um against dr = 0.25 um;
//   2. deposit: the same beam current from the deposition weights of
//      Simulation::deposit (locC / locF, Vc / Vf, the r < c_0 rule) against the
//      exact azimuthal projection on the nodes (the field must not depend on it);
//   3. no spurious radial structure: the fraction of the field in the radial
//      band 1.0 < |k_r| < 1.6 / um (the 5 um structure of the FEL 2nd harmonic)
//      stays at the level of the smooth reference.
//
//  Units: k0 = 2 pi / lambda0 with lambda0 = 10.625 nm (the FEL benchmark);
//  dzeta = lambda0 / 32, dz = 0.625 mm, undulator period 1 cm.
//
//  usage: solver_source                 run the checks (ctest solver_source)
//         solver_source <dr_um> <src> <K/k_u> <steps> <out.txt>   one profile
//                       src = z (j_z), a (beam, exact projection), b (beam, deposited)
// ============================================================================
#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include "grid.hpp"
#include "solver.hpp"
#include "zgrid.hpp"

using namespace gem;

namespace {

const double PI = 3.14159265358979323846;
const double lambda0 = 10.625e-9;          // m
const double um = 2.0 * PI / lambda0 * 1e-6; // 1 um in code units
const double ku = lambda0 / 1e-2;           // k_u / k0 (period 1 cm)

struct Profile {
  std::vector<double> r;  // faces [um]
  std::vector<cd> Er;     // demodulated E_r on the faces
};

// m = 2 profile E_r(r) at the middle slice after `steps` steps of dz
Profile run(double dr_um, char src, double Kku, int steps, double rmax_um = 200.0, int Nz = 40) {
  const int m = 2;
  GridConfig gc;
  gc.radial = GridConfig::Radial::Uniform;
  gc.stretched = false;
  gc.dr_min = gc.dr_max = dr_um * um;
  gc.Nr = int(rmax_um / dr_um + 0.5);
  RadialGrid g;
  g.build(gc);
  const int N = g.N;
  const double h = 2.0 * PI / 32.0, k = 2.0, D = 0.625e-3 * 2.0 * PI / lambda0, K = Kku * ku;
  std::vector<double> nodes(Nz + 1);
  for (int n = 0; n <= Nz; ++n) nodes[n] = n * h;
  ZetaGrid zg;
  zg.build(nodes, true, h);
  FieldSolver S;
  S.setup(g, m);
  S.prepare(D, zg);

  // radial profiles of the source (faces: j_r; centres: j_phi, j_z)
  const double sig = 10.0 * um, off = -2.4 * um;  // beam: sigma 10 um, offset -x_w
  std::vector<cd> PR(N + 1, 0.0), PP(N, 0.0), PZ(N, 0.0);
  auto J = [&](double x, double y) { return std::exp(-0.5 * ((x - off) * (x - off) + y * y) / (sig * sig)); };
  if (src == 'z') {
    for (int i = 0; i < N; ++i) PZ[i] = std::pow(g.c[i] / sig, m) * std::exp(-0.5 * g.c[i] * g.c[i] / (sig * sig));
  } else if (src == 'a') {  // exact m-coefficient (1/2pi) int j e^{-i m phi} dphi of j = J x^
    const int nph = 1024;
    auto proj = [&](double r, int comp) {
      cd s = 0;
      for (int q = 0; q < nph; ++q) {
        const double ph = 2 * PI * (q + 0.5) / nph, cs = std::cos(ph), sn = std::sin(ph);
        s += (comp == 0 ? cs : -sn) * J(r * cs, r * sn) * std::exp(cd(0.0, -m * ph));
      }
      return s / double(nph);
    };
    for (int i = 1; i <= N; ++i) PR[i] = proj(g.rho[i], 0);
    for (int i = 0; i < N; ++i) PP[i] = proj(g.c[i], 1);
  } else {  // deposited with the weights of Simulation::deposit from a fine quadrature of the beam
    const double dq = 0.025 * um, L = 6 * sig + std::fabs(off);
    const int nq = int(2 * L / dq);
    for (int ix = 0; ix < nq; ++ix)
      for (int iy = 0; iy < nq; ++iy) {
        const double x = -L + (ix + 0.5) * dq, y = -L + (iy + 0.5) * dq;
        const double W = J(x, y) * dq * dq / (2 * PI), r = std::hypot(x, y);
        const double cs = r > 0 ? x / r : 1.0, sn = r > 0 ? y / r : 0.0;
        const cd em = std::exp(cd(0.0, -m * std::atan2(y, x)));
        auto wc = g.locC(r), wf = g.locF(r);
        auto addC = [&](int i, double sp) { PP[i] += W / g.Vc[i] * sp * (-sn) * em; };
        if (wc.mode == 0) { addC(wc.a, wc.wa); addC(wc.a + 1, wc.wb); }
        else if (wc.mode == 2) addC(N - 1, 1.0);
        else addC(0, 1.0);  // r < c_0: j_phi with weight 1 (as in Simulation::deposit)
        auto addF = [&](int i, double s) { if (i > 0) PR[i] += W * s / g.Vf[i] * cs * em; };
        if (wf.mode == 0) { addF(wf.a, wf.wa); addF(wf.a + 1, wf.wb); }
        else addF(N, 1.0);
      }
  }

  struct F { std::vector<cd> Rr, Lr, Rp, Lp, Ez, Bz; };
  std::vector<F> fold(Nz + 1), fnew;
  for (auto& x : fold) {
    x.Rr.assign(N + 1, 0.0); x.Lr.assign(N + 1, 0.0); x.Bz.assign(N + 1, 0.0);
    x.Rp.assign(N, 0.0); x.Lp.assign(N, 0.0); x.Ez.assign(N, 0.0);
  }
  fnew = fold;
  const int nramp = 16;  // front ramp of the source in zeta (half-cells), one undulator period ramp in s
  auto amp = [&](int n, double s) {
    const double wz = (n > Nz - 1 - nramp) ? std::pow(std::sin(0.5 * PI * (Nz - 1 - n) / nramp), 2) : 1.0;
    const double ws = std::min(1.0, s * ku / (2 * PI));
    return wz * ws * std::exp(cd(0.0, k * (n + 0.5) * h + K * s));
  };
  struct Jc { std::vector<cd> r, p, z; };
  auto fill = [&](int n, double s, Jc& c) {
    const cd a = amp(std::max(n, 0), s);
    c.r.resize(N + 1); c.p.resize(N); c.z.resize(N);
    for (int i = 0; i <= N; ++i) c.r[i] = a * PR[i];
    for (int i = 0; i < N; ++i) { c.p[i] = a * PP[i]; c.z[i] = a * PZ[i]; }
  };
  Jc a0, b0, a1, b1;
  std::vector<cd> Jr(N + 1), Jp(N);
  double s0 = 0.0;
  for (int st = 0; st < steps; ++st) {
    const double s1 = s0 + D;
    for (int n = Nz - 1; n >= 0; --n) {
      fnew[n] = fold[n];
      fill(n - 1, s0, a0); fill(n, s0, b0); fill(n - 1, s1, a1); fill(n, s1, b1);
      for (int i = 0; i <= N; ++i) Jr[i] = a0.r[i] + b0.r[i] + a1.r[i] + b1.r[i];
      for (int i = 0; i < N; ++i) Jp[i] = a0.p[i] + b0.p[i] + a1.p[i] + b1.p[i];
      F &o = fold[n], &u = fnew[n + 1], &w = fnew[n];
      S.solveSlice(m, n, o.Rr.data(), o.Rp.data(), o.Ez.data(), o.Bz.data(), u.Rr.data(), u.Lr.data(), u.Rp.data(),
                   u.Lp.data(), u.Ez.data(), u.Bz.data(), b1.r.data(), b1.p.data(), b1.z.data(), Jr.data(), Jp.data(),
                   w.Rr.data(), w.Lr.data(), w.Rp.data(), w.Lp.data(), w.Ez.data(), w.Bz.data());
    }
    std::swap(fold, fnew);
    s0 = s1;
  }
  const int nm = Nz / 2;
  const cd ph = std::exp(cd(0.0, -(k * nm * h + K * s0)));
  Profile P;
  for (int i = 0; i <= N; ++i) {
    P.r.push_back(g.rho[i] / um);
    P.Er.push_back(0.5 * (fold[nm].Rr[i] + fold[nm].Lr[i]) * ph);
  }
  return P;
}

// max |E - E_ref| / max |E_ref| over r0 < r < r1 [um] (E_ref linearly interpolated)
double deviation(const Profile& p, const Profile& ref, double r0 = 2.0, double r1 = 40.0) {
  double mx = 0, peak = 0;
  for (const auto& e : ref.Er) peak = std::max(peak, std::abs(e));
  for (size_t i = 0; i < p.r.size(); ++i) {
    if (p.r[i] < r0 || p.r[i] > r1) continue;
    size_t j = 1;
    while (j + 1 < ref.r.size() && ref.r[j] < p.r[i]) ++j;
    const double t = (p.r[i] - ref.r[j - 1]) / (ref.r[j] - ref.r[j - 1]);
    const cd e = ref.Er[j - 1] + t * (ref.Er[j] - ref.Er[j - 1]);
    mx = std::max(mx, std::abs(p.Er[i] - e) / peak);
  }
  return mx;
}

// rms fraction of the field (3 < r < 32 um, Hann window) in the radial band 1.0 < |k_r| < 1.6 / um
double band(const Profile& p) {
  std::vector<double> r;
  std::vector<cd> x;
  for (size_t i = 0; i < p.r.size(); ++i)
    if (p.r[i] > 3 && p.r[i] < 32) { r.push_back(p.r[i]); x.push_back(p.Er[i]); }
  const size_t n = x.size();
  for (size_t i = 0; i < n; ++i) x[i] *= 0.5 - 0.5 * std::cos(2 * PI * (i + 0.5) / n);
  const double dr = r[1] - r[0], L = n * dr;
  double tot = 0, bnd = 0;
  for (const auto& v : x) tot += std::norm(v);
  // Parseval over the discrete wavenumbers 2 pi j / L
  for (int j = -int(n) / 2; j < int(n) / 2; ++j) {
    const double kk = 2 * PI * j / L;
    if (std::fabs(kk) <= 1.0 || std::fabs(kk) >= 1.6) continue;
    cd s = 0;
    for (size_t i = 0; i < n; ++i) s += x[i] * std::exp(cd(0.0, -kk * (r[i] - r[0])));
    bnd += std::norm(s) / n;
  }
  return std::sqrt(bnd / tot);
}

}  // namespace

int main(int argc, char** argv) {
  if (argc == 6) {
    const Profile p = run(std::atof(argv[1]), argv[2][0], std::atof(argv[3]), std::atoi(argv[4]));
    FILE* f = std::fopen(argv[5], "w");
    std::fprintf(f, "# r[um] |Er| Re Im   (m = 2, demodulated)  band fraction %.4e\n", band(p));
    for (size_t i = 0; i < p.r.size(); ++i)
      std::fprintf(f, "%g %.6e %.6e %.6e\n", p.r[i], std::abs(p.Er[i]), p.Er[i].real(), p.Er[i].imag());
    std::fclose(f);
    return 0;
  }
  int fail = 0;
  auto check = [&](const char* what, double v, double lim) {
    const bool ok = v <= lim;
    std::printf("%-62s %.3e  (limit %.1e)  %s\n", what, v, lim, ok ? "ok" : "FAIL");
    if (!ok) ++fail;
  };
  const int steps = 160;  // 10 undulator periods: the response to the K = 2 k_u source is stationary
  const Profile z1 = run(1.0, 'z', 2.0, steps), z4 = run(0.25, 'z', 2.0, steps);
  check("j_z, K = 2 k_u: dr = 1 um vs 0.25 um, max |dE|/max|E|", deviation(z1, z4), 0.03);
  const Profile a1 = run(1.0, 'a', -1.0, steps), a4 = run(0.25, 'a', -1.0, steps);
  check("beam j_x, K = -k_u: dr = 1 um vs 0.25 um", deviation(a1, a4), 0.03);
  // the deposit is less accurate in its first cells near the axis: 3.7 % of the peak field at dr = 1 um, 0.9 % at 0.5 um
  // (second order in dr); beyond 6 um
  // the field equals the one of the exact projection
  const Profile b1 = run(1.0, 'b', -1.0, steps), b2 = run(0.5, 'b', -1.0, steps), a2 = run(0.5, 'a', -1.0, steps);
  check("beam j_x: deposited vs exact projection, dr = 1 um, r < 5 um", deviation(b1, a1, 0.5, 5.0), 0.08);
  check("beam j_x: deposited vs exact projection, dr = 0.5 um, r < 5 um", deviation(b2, a2, 0.25, 5.0), 0.04);
  check("beam j_x: deposited vs exact projection, dr = 1 um, 6 < r < 40 um", deviation(b1, a1, 6.0, 40.0), 0.01);
  check("band 1.0-1.6/um fraction, exact projection", band(a1), 0.01);
  check("band 1.0-1.6/um fraction, deposited", band(b1), 0.01);
  check("band 1.0-1.6/um fraction, j_z", band(z1), 0.01);
  std::printf("%s\n", fail ? "solver_source: FAILED" : "solver_source: all checks passed");
  return fail ? 1 : 0;
}
