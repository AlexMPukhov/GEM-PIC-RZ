// ============================================================================
//  laser_file.cpp -- laser pulses from LASY / openPMD envelope files (see laser_file.hpp)
// ============================================================================
#include "laser_file.hpp"

#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>

#if defined(GEMPIC_HAVE_HDF5)
#include "h5read.hpp"
#endif

namespace gem {

namespace {
const cd I1(0.0, 1.0);
// weights of a trapezoidal rule on [0, (n-1) h] (units of h) with Gregory's corrections at the left end, up to
// the 4th differences: the integrands of the Hankel transforms, r f(r) and kappa F(kappa), are odd at 0, where
// the plain trapezoid has an error O(h^2) (Euler-Maclaurin); the right end is where the data have died out
double gregoryLeft(int p, int n) {
  double w = (p == 0 || p == n - 1) ? 0.5 : 1.0;
  if (n < 10) return w;
  const double a1 = 1.0 / 12, a2 = -1.0 / 24, a3 = 19.0 / 720, a4 = -3.0 / 160;
  const double c[5] = {-a1 + a2 - a3 + a4, a1 - 2 * a2 + 3 * a3 - 4 * a4, a2 - 3 * a3 + 6 * a4, a3 - 4 * a4, a4};
  if (p < 5) w += c[p];
  return w;
}

// Weights W_j of the forward Hankel transform F(kappa) = int f(r) J_nu(kappa r) r dr = sum_j W_j f_j of the data
// f_j on the uniform nodes r_j = r0 + j dr. The integral is exact for the interpolant f = r^nu g(r^2), g linear in
// u = r^2 between the nodes (the regular form of an azimuthal mode nu near the axis); g at r = 0 is extrapolated
// linearly in u from the first two nodes (nu > 0), and g is constant on [0, r0] if r0 > 0. Exactness matters: a
// quadrature rule of the integrand r f J (odd in r) leaves an error ~ dr^2 f(0) at every kappa, which the inverse
// transform (measure kappa dkappa) turns into a large error on the axis. Uses
//   int r^{nu+1} J_nu(k r) dr = r^{nu+1} J_{nu+1} / k,   int r^{nu+3} J_nu(k r) dr = r^{nu+3} J_{nu+1} / k - 2 r^{nu+2} J_{nu+2} / k^2.
std::vector<double> hankelWeights(int nu, double kap, double r0, double dr, int nr) {
  std::vector<double> W(nr, 0.0);
  // nodes: an extra node at r = 0 if the data start at r0 > 0 (index -1)
  auto rr = [&](int j) { return j < 0 ? 0.0 : r0 + j * dr; };
  auto B1 = [&](double r) {
    if (kap == 0.0) return nu == 0 ? 0.5 * r * r : 0.0;
    return std::pow(r, nu + 1) * std::cyl_bessel_j(nu + 1.0, kap * r) / kap;
  };
  auto B3 = [&](double r) {
    if (kap == 0.0) return nu == 0 ? 0.25 * r * r * r * r : 0.0;
    return std::pow(r, nu + 3) * std::cyl_bessel_j(nu + 1.0, kap * r) / kap -
           2.0 * std::pow(r, nu + 2) * std::cyl_bessel_j(nu + 2.0, kap * r) / (kap * kap);
  };
  // coefficient of g at node j -> weights of f: g_j = f_j / r_j^nu; g(0) from extrapolation
  auto addG = [&](int j, double c) {
    const double r = rr(j);
    if (j >= 0 && (r > 0.0 || nu == 0)) { W[j] += c / std::pow(r, nu); return; }
    if (j < 0) { W[0] += c / std::pow(r0, nu); return; }  // g constant = g(r0) on [0, r0]
    // node r = 0 of a mode nu > 0: g(0) = g_a + (g_a - g_b) u_a / (u_b - u_a) from the nodes 1 and 2
    const int ja = 1, jb = 2;
    if (jb >= nr) { W[ja] += c / std::pow(rr(ja), nu); return; }
    const double ua = rr(ja) * rr(ja), ub = rr(jb) * rr(jb), x = ua / (ub - ua);
    W[ja] += c * (1.0 + x) / std::pow(rr(ja), nu);
    W[jb] -= c * x / std::pow(rr(jb), nu);
  };
  const int jfirst = r0 > 0.0 ? -1 : 0;
  double b1a = B1(rr(jfirst)), b3a = B3(rr(jfirst));
  for (int j = jfirst; j < nr - 1; ++j) {
    const double b1b = B1(rr(j + 1)), b3b = B3(rr(j + 1));
    const double I1 = b1b - b1a;
    if (j < 0) {
      addG(-1, I1);  // [0, r0]: g constant = g(r0)
    } else {
      const double ua = rr(j) * rr(j), du = rr(j + 1) * rr(j + 1) - ua;
      const double Iu = (b3b - b3a) - ua * I1;  // int r^{nu+1} (u - u_j) J dr
      addG(j, I1 - Iu / du);
      addG(j + 1, Iu / du);
    }
    b1a = b1b;
    b3a = b3b;
  }
  return W;
}
std::string fmt(const char* f, double a, double b = 0, double c = 0, double d = 0) {
  char buf[256];
  std::snprintf(buf, sizeof buf, f, a, b, c, d);
  return buf;
}
}  // namespace

// ---------------------------------------------------------------- FFT
void fft(std::vector<cd>& a, int sign) {
  const size_t n = a.size();
  if (n < 2) return;
  if (n & (n - 1)) throw std::runtime_error("fft: length is not a power of two");
  for (size_t i = 1, j = 0; i < n; ++i) {
    size_t bit = n >> 1;
    for (; j & bit; bit >>= 1) j ^= bit;
    j ^= bit;
    if (i < j) std::swap(a[i], a[j]);
  }
  for (size_t len = 2; len <= n; len <<= 1) {
    const double ang = sign * 2.0 * PI / double(len);
    const cd wl(std::cos(ang), std::sin(ang));
    for (size_t i = 0; i < n; i += len) {
      cd w(1.0, 0.0);
      for (size_t k = 0; k < len / 2; ++k) {
        const cd u = a[i + k], v = a[i + k + len / 2] * w;
        a[i + k] = u + v;
        a[i + k + len / 2] = u - v;
        if ((k & 63) == 63) w = std::polar(1.0, ang * double(k + 1));  // limit the drift of the recurrence
        else w *= wl;
      }
    }
  }
}

// ---------------------------------------------------------------- envelope -> R modes
// E_x + i E_y = [A env + conj(B env)] / 2 with A = (px + i py) phi, B = (px - i py) phi, phi = e^{-i omegaL tau};
// W = E_r + i E_phi = e^{-i theta} (E_x + i E_y) = sum_n w_n e^{i n theta},
// w_n = [A a_{-n-1} + conj(B a_{n+1})] / 2;  E_r,m = (w_m + conj(w_-m)) / 2, E_phi,m = (w_m - conj(w_-m)) / (2 i)
void envelopeToR(const cd* a, int K, cd px, cd py, double omegaL, double tau, int M, cd* Rr, cd* Rp) {
  const cd phi = std::polar(1.0, -omegaL * tau);
  const cd A = (px + I1 * py) * phi, B = (px - I1 * py) * phi;
  auto ak = [&](int k) { return (k < -K || k > K) ? cd(0.0) : a[k + K]; };
  auto w = [&](int n) { return 0.5 * (A * ak(-n - 1) + std::conj(B * ak(n + 1))); };
  for (int m = 0; m <= M; ++m) {
    const cd wp = w(m), wm = std::conj(w(-m));
    Rr[m] = wp + wm;           // R = 2 E
    Rp[m] = -I1 * (wp - wm);
  }
}

// ---------------------------------------------------------------- reading
LaserFileData readLaserFile(const std::string& path, const std::string& mesh, long iteration, const Units& u) {
#if !defined(GEMPIC_HAVE_HDF5)
  (void)mesh; (void)iteration; (void)u;
  throw std::runtime_error(path + ": reading laser files needs HDF5 support (install libhdf5-dev and rebuild)");
#else
  H5Eset_auto2(H5E_DEFAULT, nullptr, nullptr);
  h5r::File F(path);
  LaserFileData D;
  D.path = path;
  D.mesh = mesh.empty() ? "laserEnvelope" : mesh;
  if (!F.hasAttr("/", "openPMD")) F.fail("not an openPMD file (no attribute 'openPMD')");
  std::string itRoot;
  D.iteration = F.iteration(iteration, itRoot);
  std::string mp = F.attrS("/", "meshesPath");
  if (mp.empty()) mp = "meshes/";
  if (mp.back() != '/') mp += '/';
  const std::string mroot = itRoot + std::to_string(D.iteration) + "/" + mp;
  const std::string p = mroot + D.mesh;
  if (!F.exists(p)) {
    std::string l;
    for (const auto& c : F.children(mroot)) l += (l.empty() ? "" : ", ") + c;
    F.fail("no mesh '" + D.mesh + "' in iteration " + std::to_string(D.iteration) + " (meshes: " +
           (l.empty() ? "none" : l) + "); a LASY envelope is 'laserEnvelope', else give file_mesh");
  }
  if (F.isGroup(p)) F.fail("mesh '" + D.mesh + "' is not a scalar record (an envelope is one complex dataset)");
  D.software = F.attrS("/", "software");
  const std::string geo = F.attrS(p, "geometry");
  if (geo == "cartesian") D.cartesian = true;
  else if (geo != "thetaMode") F.fail("mesh '" + D.mesh + "': geometry '" + geo + "' (need thetaMode or cartesian)");
  std::vector<std::string> lab = F.attrSV(p, "axisLabels");
  std::vector<double> sp = F.attrD(p, "gridSpacing"), off = F.attrD(p, "gridGlobalOffset");
  std::vector<double> pos(lab.size(), 0.0);
  if (F.hasAttr(p, "position")) pos = F.attrD(p, "position");
  if (F.attrS(p, "dataOrder") == "F") {
    std::reverse(lab.begin(), lab.end());
    std::reverse(sp.begin(), sp.end());
    std::reverse(off.begin(), off.end());
    std::reverse(pos.begin(), pos.end());
  }
  const size_t na = D.cartesian ? 3 : 2;
  if (lab.size() != na || sp.size() != na || off.size() != na || pos.size() != na)
    F.fail("mesh '" + D.mesh + "': expected " + std::to_string(na) + " axes (axisLabels, gridSpacing, gridGlobalOffset)");
  auto axis = [&](const char* name) {
    for (size_t k = 0; k < lab.size(); ++k)
      if (lab[k] == name) return int(k);
    if (std::string(name) == "t")
      F.fail("mesh '" + D.mesh + "': no time axis 't' (axisLabels " + lab[0] + ", " + lab[1] +
             (na > 2 ? ", " + lab[2] : std::string()) + "); a field given in z must be converted to t first");
    F.fail(std::string("mesh '") + D.mesh + "': no axis '" + name + "'");
  };
  const double gu = F.attr1(p, "gridUnitSI", 1.0);
  const double w0 = u.omega0_per_fs() * 1e15, Lm = u.k0inv_um() * 1e-6;
  const double E0 = MC2_MEV * 1e6 / Lm;  // m c omega0 / e [V/m]
  // envelope metadata (openPMD LaserEnvelope extension)
  if (!F.hasAttr(p, "angularFrequency")) F.fail("mesh '" + D.mesh + "' has no attribute angularFrequency");
  D.omegaL = F.attrD(p, "angularFrequency")[0] / w0;
  if (!(D.omegaL > 0)) F.fail("angularFrequency must be positive");
  if (F.hasAttr(p, "polarization")) {
    const auto v = F.attrC(p, "polarization");
    if (v.size() != 2) F.fail("polarization must have 2 entries");
    D.px = v[0];
    D.py = v[1];
  }
  const double pn = std::sqrt(std::norm(D.px) + std::norm(D.py));
  if (!(pn > 0)) F.fail("polarization vector is zero");
  D.px /= pn;
  D.py /= pn;
  D.field = F.hasAttr(p, "envelopeField") ? F.attrS(p, "envelopeField") : "normalized_vector_potential";
  const bool vecpot = (D.field == "normalized_vector_potential");
  if (!vecpot && D.field != "electric_field")
    F.fail("envelopeField '" + D.field + "' (need electric_field or normalized_vector_potential)");
  const double uSI = F.attr1(p, "unitSI", 1.0);

  // data, reordered to [comp][t][r] or [t][y][x]
  const std::vector<hsize_t> dims = F.dims(p);
  const std::vector<cd> raw = F.dataC(p);
  const int it = axis("t");
  if (D.cartesian) {
    if (dims.size() != 3) F.fail("cartesian envelope: dataset must have 3 dimensions");
    const int iy = axis("y"), ix = axis("x");
    D.nt = int(dims[it]); D.ny = int(dims[iy]); D.nx = int(dims[ix]);
    D.t0 = (off[it] + pos[it] * sp[it]) * gu * w0; D.dt = sp[it] * gu * w0;
    D.y0 = (off[iy] + pos[iy] * sp[iy]) * gu / Lm; D.dy = sp[iy] * gu / Lm;
    D.x0 = (off[ix] + pos[ix] * sp[ix]) * gu / Lm; D.dx = sp[ix] * gu / Lm;
    size_t st[3] = {dims[1] * dims[2], dims[2], 1};
    D.env.resize(size_t(D.nt) * D.ny * D.nx);
    for (int l = 0; l < D.nt; ++l)
      for (int y = 0; y < D.ny; ++y)
        for (int x = 0; x < D.nx; ++x)
          D.env[(size_t(l) * D.ny + y) * D.nx + x] = raw[l * st[it] + y * st[iy] + x * st[ix]];
  } else {
    if (dims.size() != 3) F.fail("thetaMode envelope: dataset must have 3 dimensions (components, 2 axes)");
    const int ir = axis("r");
    D.ncomp = int(dims[0]);
    if (D.ncomp % 2 != 1) F.fail("thetaMode envelope: the number of components must be odd (2 K + 1)");
    D.nt = int(dims[1 + it]); D.nr = int(dims[1 + ir]);
    D.t0 = (off[it] + pos[it] * sp[it]) * gu * w0; D.dt = sp[it] * gu * w0;
    D.r0 = (off[ir] + pos[ir] * sp[ir]) * gu / Lm; D.dr = sp[ir] * gu / Lm;
    if (D.r0 < -1e-9 * D.dr) F.fail("thetaMode envelope: the radial axis starts below r = 0");
    D.r0 = std::max(0.0, D.r0);
    const size_t st[2] = {dims[2], 1};
    D.env.resize(size_t(D.ncomp) * D.nt * D.nr);
    for (int c = 0; c < D.ncomp; ++c)
      for (int l = 0; l < D.nt; ++l)
        for (int j = 0; j < D.nr; ++j)
          D.env[(size_t(c) * D.nt + l) * D.nr + j] = raw[c * dims[1] * dims[2] + l * st[it] + j * st[ir]];
  }
  if (D.nt < 4 || !(D.dt > 0)) F.fail("the time axis needs at least 4 points and a positive spacing");
  // units: E [V/m] -> code; normalized vector potential -> E = i omegaL a - d a / dt
  if (!vecpot) {
    for (auto& v : D.env) v *= uSI / E0;
  } else {
    const size_t inner = D.cartesian ? size_t(D.ny) * D.nx : size_t(D.nr);
    const int outer = D.cartesian ? 1 : D.ncomp;
    std::vector<cd> col(D.nt);
    for (int c = 0; c < outer; ++c)
      for (size_t q = 0; q < inner; ++q) {
        const size_t base = size_t(c) * D.nt * inner + q;
        for (int l = 0; l < D.nt; ++l) col[l] = D.env[base + l * inner] * uSI;
        for (int l = 0; l < D.nt; ++l) {
          const int la = std::max(0, l - 1), lb = std::min(D.nt - 1, l + 1);
          const cd dadt = (col[lb] - col[la]) / ((lb - la) * D.dt);
          D.env[base + l * inner] = I1 * D.omegaL * col[l] - dadt;
        }
      }
  }
  // peak |env| over theta
  if (D.cartesian) {
    for (const auto& v : D.env) D.peak = std::max(D.peak, std::abs(v));
  } else {
    const int Kf = (D.ncomp - 1) / 2, nth = 4 * Kf + 4;
    for (int l = 0; l < D.nt; ++l)
      for (int j = 0; j < D.nr; ++j)
        for (int q = 0; q < nth; ++q) {
          const double th = 2.0 * PI * q / nth;
          cd e = D.env[(size_t(0) * D.nt + l) * D.nr + j];
          for (int m = 1; m <= Kf; ++m)
            e += D.env[(size_t(2 * m - 1) * D.nt + l) * D.nr + j] * std::cos(m * th) +
                 D.env[(size_t(2 * m) * D.nt + l) * D.nr + j] * std::sin(m * th);
          D.peak = std::max(D.peak, std::abs(e));
        }
  }
  if (!(D.peak > 0)) F.fail("the envelope is zero everywhere");
  return D;
#endif
}

// ---------------------------------------------------------------- field of a file laser
LaserFileField::LaserFileField(const LaserFileData& d, const Params& p, int M_, double rNeed_)
    : P(p), M(M_), K(M_ + 1), omegaL(d.omegaL), px(d.px), py(d.py), rNeed(rNeed_) {
  nt = d.nt;
  t0 = d.t0;
  dt = d.dt;
  const int NK = 2 * K + 1;
  double pAll = 0.0, pKept = 0.0;
  if (!d.cartesian) {
    nr = d.nr;
    r0 = d.r0;
    dr = d.dr;
    A.assign(size_t(NK) * nt * nr, cd(0.0));
    const int Kf = (d.ncomp - 1) / 2;
    for (int l = 0; l < nt; ++l)
      for (int j = 0; j < nr; ++j) {
        const double wr = std::max(r0 + j * dr, 0.25 * dr);
        auto put = [&](int k, cd v) {
          const double e = std::norm(v) * wr;
          pAll += e;
          if (k >= -K && k <= K) { A[(size_t(k + K) * nt + l) * nr + j] = v; pKept += e; }
        };
        put(0, d.env[(size_t(0) * nt + l) * nr + j]);
        for (int m = 1; m <= Kf; ++m) {
          const cd C = d.env[(size_t(2 * m - 1) * nt + l) * nr + j], S = d.env[(size_t(2 * m) * nt + l) * nr + j];
          put(m, 0.5 * (C + I1 * S));
          put(-m, 0.5 * (C - I1 * S));
        }
      }
  } else {
    // polar resampling around the axis x = y = 0: largest disk inside the file's grid
    dr = std::min(d.dx, d.dy);
    const double R = std::min({d.x0 + (d.nx - 1) * d.dx, -d.x0, d.y0 + (d.ny - 1) * d.dy, -d.y0});
    if (!(R > dr)) throw std::runtime_error(d.path + ": the cartesian grid must contain the axis x = y = 0 well inside");
    r0 = 0.0;
    nr = int(std::floor(R / dr + 1e-9)) + 1;
    const int nth = p.n_theta > 0 ? p.n_theta : std::max(32, 4 * K + 4);
    A.assign(size_t(NK) * nt * nr, cd(0.0));
    std::vector<double> cs(nth), sn(nth);
    for (int q = 0; q < nth; ++q) { cs[q] = std::cos(2 * PI * q / nth); sn[q] = std::sin(2 * PI * q / nth); }
    for (int l = 0; l < nt; ++l)
      for (int j = 0; j < nr; ++j) {
        const double r = j * dr, wr = std::max(r, 0.25 * dr);
        std::vector<cd> e(nth);
        for (int q = 0; q < nth; ++q) {
          const double fx = (r * cs[q] - d.x0) / d.dx, fy = (r * sn[q] - d.y0) / d.dy;
          const int ix = std::min(std::max(int(std::floor(fx)), 0), d.nx - 2);
          const int iy = std::min(std::max(int(std::floor(fy)), 0), d.ny - 2);
          const double wx = std::min(1.0, std::max(0.0, fx - ix)), wy = std::min(1.0, std::max(0.0, fy - iy));
          auto at = [&](int y, int x) { return d.env[(size_t(l) * d.ny + y) * d.nx + x]; };
          e[q] = (1 - wy) * ((1 - wx) * at(iy, ix) + wx * at(iy, ix + 1)) +
                 wy * ((1 - wx) * at(iy + 1, ix) + wx * at(iy + 1, ix + 1));
          pAll += std::norm(e[q]) * wr / nth;
        }
        for (int k = -K; k <= K; ++k) {  // env = sum a_k e^{-i k theta}
          cd s = 0.0;
          for (int q = 0; q < nth; ++q) s += e[q] * std::polar(1.0, k * 2 * PI * q / nth);
          s /= double(nth);
          A[(size_t(k + K) * nt + l) * nr + j] = s;
          pKept += std::norm(s) * wr;
        }
      }
    info.push_back("polar resampling of the cartesian grid: " + std::to_string(nr) + " radii up to the largest disk "
                  "around the axis, " + std::to_string(nth) + " angles");
  }
  modeLoss = pAll > 0 ? std::max(0.0, 1.0 - pKept / pAll) : 0.0;
  for (auto& v : A) v *= P.scale;
  double amax = 0.0, aedge = 0.0;
  for (int a = 0; a < NK; ++a) {
    double m = 0.0;
    for (int l = 0; l < nt; ++l)
      for (int j = 0; j < nr; ++j) {
        const double v = std::abs(A[(size_t(a) * nt + l) * nr + j]);
        m = std::max(m, v);
        if (l == 0 || l == nt - 1 || j == nr - 1) aedge = std::max(aedge, v);
      }
    amax = std::max(amax, m);
    if (m > 0) active.push_back(a);
  }
  edgeField = amax > 0 ? aedge / amax : 0.0;
}

void LaserFileField::interp(double r, double t, cd* a) const {
  const int NK = 2 * K + 1;
  for (int k = 0; k < NK; ++k) a[k] = 0.0;
  const double x = (t - t0) / dt;
  if (!(x >= 0.0) || x > nt - 1) return;
  const int l = std::min(int(x), nt - 2);
  const double wl = x - l;
  const double y = (r - r0) / dr;
  if (y > nr - 1) return;
  int j = 0;
  double wj = 0.0, f = 1.0;
  if (y < 0.0) {  // inside the first node (cell-centred files): regular towards the axis
    f = r0 > 0 ? r / r0 : 1.0;
  } else {
    j = std::min(int(y), nr - 2);
    wj = y - j;
  }
  for (int ka : active) {
    const cd* B = &A[(size_t(ka) * nt) * nr];
    cd v = (1 - wl) * ((1 - wj) * B[size_t(l) * nr + j] + wj * B[size_t(l) * nr + j + 1]) +
           wl * ((1 - wj) * B[size_t(l + 1) * nr + j] + wj * B[size_t(l + 1) * nr + j + 1]);
    if (y < 0.0 && ka != K) v *= f;  // modes k != 0 vanish on the axis
    a[ka] = v;
  }
}

void LaserFileField::prepare() const {
  if (prepared) return;
  prepared = true;
  const int na = int(active.size());
  nt2 = 1;
  while (nt2 < std::max(1, P.padding) * nt) nt2 <<= 1;
  const int off = (nt2 - nt) / 2;
  Tstart = t0 - off * dt;
  dOm = 2.0 * PI / (nt2 * dt);
  // FFT in t of every active mode and radius: F[(a * nr + j) * nt2 + q]
  std::vector<cd> F(size_t(na) * nr * nt2);
#pragma omp parallel for schedule(dynamic)
  for (int aj = 0; aj < na * nr; ++aj) {
    const int a = active[aj / nr], j = aj % nr;
    std::vector<cd> s(nt2, cd(0.0));
    for (int l = 0; l < nt; ++l) s[off + l] = A[(size_t(a) * nt + l) * nr + j];
    fft(s, +1);  // F_q = sum_l f_l e^{+i Omega_q (t_l - Tstart)}
    std::copy(s.begin(), s.end(), F.begin() + size_t(aj) * nt2);
  }
  // spectral band: frequencies with power above 1e-12 of the largest that propagate (omega > omega_p)
  std::vector<double> pw(nt2, 0.0);
  for (int a = 0; a < na; ++a)
    for (int j = 0; j < nr; ++j) {
      const double wr = std::max(r0 + j * dr, 0.25 * dr);
      for (int q = 0; q < nt2; ++q) pw[q] += std::norm(F[(size_t(a) * nr + j) * nt2 + q]) * wr;
    }
  double pmax = 0.0;
  for (double v : pw) pmax = std::max(pmax, v);
  int smin = nt2, smax = -nt2;
  const double wp = std::sqrt(std::max(0.0, P.background));
  for (int q = 0; q < nt2; ++q) {
    const int qs = q < nt2 / 2 ? q : q - nt2;
    if (pw[q] > 1e-12 * pmax && omegaL + qs * dOm > wp) { smin = std::min(smin, qs); smax = std::max(smax, qs); }
  }
  if (smax < smin) throw std::runtime_error("laser file: no propagating frequency in the spectrum");
  qlo = smin;
  nq = smax - smin + 1;
  // Hankel grid: kappa_p = p dk; the sampling is periodic in r with 2 pi / dk: keep the images beyond the
  // radii in use (file radius R_f + the larger of R_f and the largest radius evaluated), with a margin of 2
  const double Rf = r0 + (nr - 1) * dr, Rn = std::max(Rf, rNeed);
  dk = 2.0 * PI / (2.0 * (Rf + Rn));
  const double wmax = omegaL + (qlo + nq - 1) * dOm;
  const double kmax = std::min(PI / dr, std::sqrt(std::max(0.0, wmax * wmax - P.background)));
  nk = std::max(2, int(std::floor(kmax / dk)) + 1);
  G.assign(size_t(na) * nq * nk, cd(0.0));
  for (int a = 0; a < na; ++a) {
    const int order = std::abs(active[a] - K);
#pragma omp parallel for schedule(dynamic)
    for (int pk = 0; pk < nk; ++pk) {
      const double kap = pk * dk, wk = gregoryLeft(pk, nk);
      const std::vector<double> W = hankelWeights(order, kap, r0, dr, nr);
      for (int q = 0; q < nq; ++q) {
        const int qi = ((qlo + q) % nt2 + nt2) % nt2;
        cd s = 0.0;
        for (int j = 0; j < nr; ++j) s += F[(size_t(a) * nr + j) * nt2 + qi] * W[j];
        G[(size_t(a) * nq + q) * nk + pk] = s * kap * dk * wk;  // inverse weights folded in
      }
    }
  }
  beta.assign(size_t(nq) * nk, std::numeric_limits<double>::quiet_NaN());
  for (int q = 0; q < nq; ++q) {
    const double w = omegaL + (qlo + q) * dOm;
    for (int pk = 0; pk < nk; ++pk) {
      const double kz2 = w * w - P.background - pk * dk * pk * dk;
      if (w > 0 && kz2 > 0) beta[size_t(q) * nk + pk] = std::sqrt(kz2) - w;
    }
  }
  const double lam = 2.0 * PI / omegaL;
  info.push_back("propagator: FFT " + std::to_string(nt2) + " points (band " + std::to_string(nq) + " frequencies, " +
                 fmt("omega/omega_L = %.3f .. %.3f", (omegaL + qlo * dOm) / omegaL, wmax / omegaL) + "), " +
                 std::to_string(nk) + " transverse wavenumbers" + fmt(" (k_perp <= %.3g k_L, dk = %.3g / lambda_L)",
                 kmax / omegaL, dk * lam));
}

void LaserFileField::field(double s, const std::vector<double>& radii, const std::vector<double>& zetas, int comp,
                           std::vector<cd>& out, bool direct) const {
  const int NK = 2 * K + 1, nz = int(zetas.size()), nrad = int(radii.size()), M1 = M + 1;
  out.assign(size_t(nrad) * nz * M1, cd(0.0));
  if (direct) {
#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < nrad; ++i) {
      std::vector<cd> a(NK), Rr(M1), Rp(M1);
      for (int n = 0; n < nz; ++n) {
        const double t = tau(zetas[n]);
        interp(radii[i], t, a.data());
        envelopeToR(a.data(), K, px, py, omegaL, t, M, Rr.data(), Rp.data());
        for (int m = 0; m <= M; ++m) out[(size_t(i) * nz + n) * M1 + m] = comp == 0 ? Rr[m] : Rp[m];
      }
    }
    return;
  }
  prepare();
  const int na = int(active.size());
  const double ds = s - P.file_z;
  std::vector<cd> ph(size_t(nq) * nk);
  for (size_t k = 0; k < ph.size(); ++k) ph[k] = std::isnan(beta[k]) ? cd(0.0) : std::polar(1.0, beta[k] * ds);
  const double Twin = nt2 * dt, norm = 1.0 / nt2;
#pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < nrad; ++i) {
    // H[a][q] = sum_p G J(kappa_p r) e^{i beta (s - z_f)}
    std::vector<cd> H(size_t(na) * nq, cd(0.0)), a(NK), Rr(M1), Rp(M1);
    std::vector<double> J(nk);
    for (int ai = 0; ai < na; ++ai) {
      const int order = std::abs(active[ai] - K);
      for (int pk = 0; pk < nk; ++pk) J[pk] = std::cyl_bessel_j(double(order), pk * dk * radii[i]);
      for (int q = 0; q < nq; ++q) {
        const cd* g = &G[(size_t(ai) * nq + q) * nk];
        const cd* e = &ph[size_t(q) * nk];
        cd sum = 0.0;
        for (int pk = 0; pk < nk; ++pk) sum += g[pk] * e[pk] * J[pk];
        H[size_t(ai) * nq + q] = sum;
      }
    }
    for (int n = 0; n < nz; ++n) {
      const double t = tau(zetas[n]), x = t - Tstart;
      for (int k = 0; k < NK; ++k) a[k] = 0.0;
      if (x >= 0.0 && x < Twin) {
        const cd step = std::polar(1.0, -dOm * x);
        for (int ai = 0; ai < na; ++ai) {
          cd e = std::polar(1.0, -qlo * dOm * x), sum = 0.0;
          const cd* h = &H[size_t(ai) * nq];
          for (int q = 0; q < nq; ++q) {
            sum += h[q] * e;
            e *= step;
          }
          a[active[ai]] = sum * norm;
        }
      }
      envelopeToR(a.data(), K, px, py, omegaL, t, M, Rr.data(), Rp.data());
      for (int m = 0; m <= M; ++m) out[(size_t(i) * nz + n) * M1 + m] = comp == 0 ? Rr[m] : Rp[m];
    }
  }
}

}  // namespace gem
