// ============================================================================
//  h5_extract.cpp -- read a GEM-PIC-RZ HDF5 file (data_NNNNN.h5) and print
//  ASCII tables; also an example of how to read the files from C++.
//
//  h5_extract <file> info                  attributes and datasets
//  h5_extract <file> axis                  zeta[um] Ex Ez ne/n0 Zni/n0 on the axis
//  h5_extract <file> map [stride]          zeta[um] x[um] Ex Ez ne/n0 Zni/n0 in the
//                                          plane phi = 0 / pi (gnuplot pm3d blocks)
//  h5_extract <file> lineout <r_um> <phi_deg>
//                                          zeta[um] Ex Ey Ez Bx By Bz at (r, phi)
//  h5_extract <file> beam                  x y zeta[um] ux uy uz q[pC] species
//  h5_extract <file> spectrum zlo zhi [k1 k2 nk]  spectral energy density dW/dk of the modes m >= 1
//                                          (radiation) over zlo <= zeta < zhi [um], k in units of k0
//  h5_extract <file> energy [zlo zhi]      field energy (|E|^2+|B|^2)/2 dV of the window (or of the
//                                          slices zlo <= zeta < zhi [um]), per mode group, and its
//                                          zeta-Nyquist part (grid-scale odd-even noise, e.g.
//                                          numerical Cherenkov radiation), in code units
//
//  The tables of `axis`, `map` and `beam` have the format of the text output
//  (axis_/map_/beam_NNNNN.txt). Fields are in m c omega0 / e (multiply by
//  the file attribute E0_V_per_m for V/m, B0_T for T).
//
//  Field reconstruction (function fieldAt) follows the particle gather of the
//  code: linear interpolation in r between the grid points of every mode,
//  regular behaviour ~ r^m inside the first centre, the m = 1 transverse
//  field on the axis extrapolated linearly, and
//      F(r, phi) = F_0 + 2 Re sum_m F_m exp(i m phi).
//
//  Build: part of the CMake / Makefile build when HDF5 is found, or
//      g++ -O2 -std=c++17 h5_extract.cpp $(pkg-config --cflags --libs hdf5)
// ============================================================================
#include <hdf5.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <stdexcept>
#include <string>
#include <vector>

using cd = std::complex<double>;

// ---------------------------------------------------------------- reading
struct File {
  hid_t id;
  explicit File(const char* path) : id(H5Fopen(path, H5F_ACC_RDONLY, H5P_DEFAULT)) {
    if (id < 0) throw std::runtime_error(std::string("cannot open ") + path);
  }
  ~File() { H5Fclose(id); }

  double attrD(const char* name, const char* obj = "/") const {
    double v = 0;
    hid_t a = H5Aopen_by_name(id, obj, name, H5P_DEFAULT, H5P_DEFAULT);
    if (a < 0 || H5Aread(a, H5T_NATIVE_DOUBLE, &v) < 0) throw std::runtime_error(std::string("attribute ") + name);
    H5Aclose(a);
    return v;
  }
  long long attrI(const char* name) const {
    long long v = 0;
    hid_t a = H5Aopen_by_name(id, "/", name, H5P_DEFAULT, H5P_DEFAULT);
    if (a < 0 || H5Aread(a, H5T_NATIVE_LLONG, &v) < 0) throw std::runtime_error(std::string("attribute ") + name);
    H5Aclose(a);
    return v;
  }
  std::string attrS(const char* name, const char* obj = "/") const {
    hid_t a = H5Aopen_by_name(id, obj, name, H5P_DEFAULT, H5P_DEFAULT);
    if (a < 0) return "";
    hid_t t = H5Aget_type(a);
    std::string s(H5Tget_size(t), '\0');
    H5Aread(a, t, &s[0]);
    H5Tclose(t);
    H5Aclose(a);
    return s.c_str();
  }
  bool exists(const char* path) const { return H5Lexists(id, path, H5P_DEFAULT) > 0; }

  // dataset -> vector (converted to double / complex<double>), dims returned
  template <class T>
  std::vector<T> read(const char* path, std::vector<hsize_t>* dims = nullptr) const {
    hid_t ds = H5Dopen2(id, path, H5P_DEFAULT);
    if (ds < 0) throw std::runtime_error(std::string("no dataset ") + path);
    hid_t sp = H5Dget_space(ds);
    std::vector<hsize_t> d(H5Sget_simple_extent_ndims(sp));
    H5Sget_simple_extent_dims(sp, d.data(), nullptr);
    size_t n = 1;
    for (hsize_t x : d) n *= x;
    std::vector<T> v(n);
    hid_t mt = memType<T>();
    if (n > 0 && H5Dread(ds, mt, H5S_ALL, H5S_ALL, H5P_DEFAULT, v.data()) < 0)
      throw std::runtime_error(std::string("cannot read ") + path);
    if (H5Tget_class(mt) == H5T_COMPOUND) H5Tclose(mt);
    H5Sclose(sp);
    H5Dclose(ds);
    if (dims) *dims = d;
    return v;
  }
  template <class T> static hid_t memType();
};
template <> hid_t File::memType<double>() { return H5T_NATIVE_DOUBLE; }
template <> hid_t File::memType<int>() { return H5T_NATIVE_INT; }
template <> hid_t File::memType<cd>() {  // {r, i} compound = layout of std::complex<double>
  hid_t t = H5Tcreate(H5T_COMPOUND, sizeof(cd));
  H5Tinsert(t, "r", 0, H5T_NATIVE_DOUBLE);
  H5Tinsert(t, "i", sizeof(double), H5T_NATIVE_DOUBLE);
  return t;
}

// ---------------------------------------------------------------- field reconstruction
struct Data {
  int M, N, Nz;
  double um;  // um: length unit in micrometres
  std::vector<double> rho, c, zs;  // zs: zeta of the slices
  std::vector<cd> Er, Bphi, Bz, Ephi, Br, Ez;  // [m][n][i]
  std::vector<double> ne, ni;                  // [k][i]
  cd face(const std::vector<cd>& F, int m, int n, int i) const { return F[(size_t(m) * (Nz + 1) + n) * (N + 1) + i]; }
  cd cent(const std::vector<cd>& F, int m, int n, int i) const { return F[(size_t(m) * (Nz + 1) + n) * N + i]; }

  explicit Data(const File& f) {
    M = int(f.attrI("modes")) - 1;
    N = int(f.attrI("Nr"));
    Nz = int(f.attrI("Nzeta"));
    um = f.attrD("k0inv_m") * 1e6;
    rho = f.read<double>("/grid/r_faces");
    c = f.read<double>("/grid/r_centres");
    zs = f.read<double>("/grid/zeta_slices");
    Er = f.read<cd>("/fields/Er"); Bphi = f.read<cd>("/fields/Bphi"); Bz = f.read<cd>("/fields/Bz");
    Ephi = f.read<cd>("/fields/Ephi"); Br = f.read<cd>("/fields/Br"); Ez = f.read<cd>("/fields/Ez");
    ne = f.read<double>("/densities/ne_over_n0");
    ni = f.read<double>("/densities/Zni_over_n0");
  }

  // Cartesian E, B at (x, y) on slice n
  void fieldAt(int n, double x, double y, double out[6]) const {
    const double r = std::sqrt(x * x + y * y);
    double cs = 1.0, sn = 0.0;
    if (r > 0.0) { cs = x / r; sn = y / r; }
    // linear weights: centres (mode 1: r < c_0, mode 2: beyond the last centre), faces
    int ac = 0, modeC = 0, af = 0, modeF = 0;
    double wc = 0, wf = 0;
    if (r < c[0]) modeC = 1;
    else if (r >= c[N - 1]) { modeC = 2; ac = N - 1; }
    else { ac = int(std::upper_bound(c.begin(), c.begin() + N, r) - c.begin()) - 1; wc = (r - c[ac]) / (c[ac + 1] - c[ac]); }
    if (r >= rho[N]) { modeF = 2; af = N; }
    else { af = int(std::upper_bound(rho.begin(), rho.end(), r) - rho.begin()) - 1; wf = (r - rho[af]) / (rho[af + 1] - rho[af]); }
    const cd I(0.0, 1.0), eip(cs, sn);
    cd em(1.0, 0.0);
    double Er_ = 0, Ep_ = 0, Ez_ = 0, Br_ = 0, Bp_ = 0, Bz_ = 0;
    for (int m = 0; m <= M; ++m) {
      auto lin = [&](cd f0, cd f1) { return (c[1] * f0 - c[0] * f1) / (c[1] - c[0]); };  // extrapolation to r = 0
      const cd ep0 = (m == 1) ? lin(cent(Ephi, m, n, 0), cent(Ephi, m, n, 1)) : cd(0.0);
      const cd br0 = (m == 1) ? lin(cent(Br, m, n, 0), cent(Br, m, n, 1)) : cd(0.0);
      cd vep, vbr, vez;
      if (modeC == 0) {
        vep = (1 - wc) * cent(Ephi, m, n, ac) + wc * cent(Ephi, m, n, ac + 1);
        vbr = (1 - wc) * cent(Br, m, n, ac) + wc * cent(Br, m, n, ac + 1);
        vez = (1 - wc) * cent(Ez, m, n, ac) + wc * cent(Ez, m, n, ac + 1);
      } else if (modeC == 2) {
        vep = cent(Ephi, m, n, ac); vbr = cent(Br, m, n, ac); vez = cent(Ez, m, n, ac);
      } else {  // inside the first centre: regular behaviour on the axis
        const double xr = r / c[0];
        vez = cent(Ez, m, n, 0) * ((m == 0) ? 1.0 : std::pow(xr, m));
        if (m == 0) { vep = cent(Ephi, m, n, 0) * xr; vbr = cent(Br, m, n, 0) * xr; }
        else if (m == 1) { vep = ep0 + (cent(Ephi, m, n, 0) - ep0) * xr; vbr = br0 + (cent(Br, m, n, 0) - br0) * xr; }
        else { const double fx = std::pow(xr, m - 1); vep = cent(Ephi, m, n, 0) * fx; vbr = cent(Br, m, n, 0) * fx; }
      }
      // faces; on the axis (face 0) E_r = -i E_phi, B_phi = i B_r for m = 1, else 0
      auto ErF = [&](int i) { return i == 0 ? ((m == 1) ? -I * ep0 : cd(0.0)) : face(Er, m, n, i); };
      auto BpF = [&](int i) { return i == 0 ? ((m == 1) ? I * br0 : cd(0.0)) : face(Bphi, m, n, i); };
      cd ver, vbp, vbz;
      if (modeF == 0) {
        ver = (1 - wf) * ErF(af) + wf * ErF(af + 1);
        vbp = (1 - wf) * BpF(af) + wf * BpF(af + 1);
        vbz = (1 - wf) * face(Bz, m, n, af) + wf * face(Bz, m, n, af + 1);
      } else {
        ver = ErF(N); vbp = BpF(N); vbz = face(Bz, m, n, N);
      }
      const double fac = (m == 0) ? 1.0 : 2.0;
      Er_ += fac * std::real(ver * em); Ep_ += fac * std::real(vep * em); Ez_ += fac * std::real(vez * em);
      Br_ += fac * std::real(vbr * em); Bp_ += fac * std::real(vbp * em); Bz_ += fac * std::real(vbz * em);
      em *= eip;
    }
    out[0] = Er_ * cs - Ep_ * sn; out[1] = Er_ * sn + Ep_ * cs; out[2] = Ez_;
    out[3] = Br_ * cs - Bp_ * sn; out[4] = Br_ * sn + Bp_ * cs; out[5] = Bz_;
  }
  // field energy (code units) and its Nyquist part hp_n = (2F_n - F_{n-1} - F_{n+1})/4
  // (slices with zlo <= zeta_n < zhi; per mode group: m = 0 and m >= 1)
  void energy(double& total, double& nyquist, double& w0, double& w1, double zlo = -1e300,
              double zhi = 1e300) const {
    std::vector<double> Af(N + 1), Ac(N);
    Af[0] = 0.5 * c[0] * c[0];
    for (int i = 1; i <= N; ++i) Af[i] = rho[i] * ((i < N ? c[i] : 2 * rho[N] - c[N - 1]) - c[i - 1]);
    for (int i = 0; i < N; ++i) Ac[i] = c[i] * (rho[i + 1] - rho[i]);
    total = nyquist = w0 = w1 = 0;
    for (int m = 0; m <= M; ++m) {
      const double fac = (m == 0) ? 1.0 : 2.0;  // <|F|^2>_phi = |F_0|^2 + 2 sum |F_m|^2
      double& wm = (m == 0) ? w0 : w1;
      auto acc = [&](const std::vector<cd>& F, int L, const std::vector<double>& A) {
        for (int n = 0; n < Nz; ++n) {
          if (zs[n] < zlo || zs[n] >= zhi) continue;
          for (int i = (L == N + 1 ? 1 : 0); i < L; ++i) {
            const size_t j = (size_t(m) * (Nz + 1) + n) * L + i;
            const double hn = zs[n + 1] - zs[n];  // length of half-cell n
            const double e = fac * hn * A[i] * std::norm(F[j]);
            total += e;
            wm += e;
            if (n > 0) nyquist += fac * hn * A[i] * std::norm(0.25 * (2.0 * F[j] - F[j - L] - F[j + L]));
          }
        }
      };
      acc(Er, N + 1, Af); acc(Bphi, N + 1, Af); acc(Bz, N + 1, Af);
      acc(Ephi, N, Ac); acc(Br, N, Ac); acc(Ez, N, Ac);
    }
    total *= M_PI;  // (1/2) * 2 pi
    nyquist *= M_PI;
    w0 *= M_PI;
    w1 *= M_PI;
  }
  // spectral energy density dW/dk of the modes m >= 1 (radiation of a centred
  // beam) over the slices zlo <= zeta < zhi: W(m>=1) = int_0^inf dW/dk dk
  double spectrum(double k, double zlo, double zhi) const {
    std::vector<double> Af(N + 1), Ac(N);
    Af[0] = 0.5 * c[0] * c[0];
    for (int i = 1; i <= N; ++i) Af[i] = rho[i] * ((i < N ? c[i] : 2 * rho[N] - c[N - 1]) - c[i - 1]);
    for (int i = 0; i < N; ++i) Ac[i] = c[i] * (rho[i + 1] - rho[i]);
    std::vector<cd> ep(Nz), em(Nz);
    for (int n = 0; n < Nz; ++n) {
      const double hn = zs[n + 1] - zs[n];
      ep[n] = std::polar(hn, -k * zs[n]);
      em[n] = std::polar(hn, k * zs[n]);
    }
    double S = 0;
    for (int m = 1; m <= M; ++m) {
      auto acc = [&](const std::vector<cd>& F, int L, const std::vector<double>& A) {
        for (int i = (L == N + 1 ? 1 : 0); i < L; ++i) {
          cd a = 0, b = 0;
          for (int n = 0; n < Nz; ++n) {
            if (zs[n] < zlo || zs[n] >= zhi) continue;
            const cd f = F[(size_t(m) * (Nz + 1) + n) * L + i];
            a += f * ep[n];
            b += f * em[n];
          }
          S += 2.0 * A[i] * (std::norm(a) + std::norm(b));
        }
      };
      acc(Er, N + 1, Af); acc(Bphi, N + 1, Af); acc(Bz, N + 1, Af);
      acc(Ephi, N, Ac); acc(Br, N, Ac); acc(Ez, N, Ac);
    }
    return 0.5 * S;  // pi / (2 pi)
  }
  double nE(int k, int i) const { return ne[size_t(std::min(k, Nz - 1)) * N + i]; }
  double nI(int k, int i) const { return ni[size_t(std::min(k, Nz - 1)) * N + i]; }
};

// ---------------------------------------------------------------- info
herr_t listObj(hid_t, const char* name, const H5O_info_t* info, void*) {
  if (info->type == H5O_TYPE_DATASET) std::printf("  dataset  /%s\n", name);
  else if (info->type == H5O_TYPE_GROUP) std::printf("  group    /%s\n", name[0] == '.' ? "" : name);
  return 0;
}

int main(int argc, char** argv) {
  if (argc < 3) {
    std::fprintf(stderr, "usage: h5_extract <file> info | axis | map [stride] | lineout <r_um> <phi_deg> | beam | energy [zeta_lo_um zeta_hi_um]\n"
                         "       | spectrum <zeta_lo_um> <zeta_hi_um> [k_min k_max nk]\n");
    return 2;
  }
  try {
    File f(argv[1]);
    const std::string what = argv[2];
    if (what == "info") {
      std::printf("%s: %s, step %lld, z = %.6g mm, %lld modes, Nr = %lld, Nzeta = %lld\n", argv[1],
                  f.attrS("code").c_str(), f.attrI("step"), f.attrD("z_m") * 1e3, f.attrI("modes"), f.attrI("Nr"),
                  f.attrI("Nzeta"));
      std::printf("  lambda0 = %.6g um, n0 = %.6g cm^-3, E0 = %.6g V/m, backend: %s\n", f.attrD("lambda0_m") * 1e6,
                  f.attrD("n0_m3") * 1e-6, f.attrD("E0_V_per_m"), f.attrS("backend").c_str());
#if H5_VERSION_GE(1, 12, 0)
      H5Ovisit3(f.id, H5_INDEX_NAME, H5_ITER_NATIVE, listObj, nullptr, H5O_INFO_BASIC);
#else
      H5Ovisit(f.id, H5_INDEX_NAME, H5_ITER_NATIVE, listObj, nullptr);
#endif
      return 0;
    }
    if (what == "beam") {
      if (!f.exists("/particles/beam")) throw std::runtime_error("no particles in this file (output.particles = false)");
      const double um = f.attrD("k0inv_m") * 1e6, pC = f.attrD("charge_unit_C") * 1e12;
      auto x = f.read<double>("/particles/beam/x"), y = f.read<double>("/particles/beam/y");
      auto z = f.read<double>("/particles/beam/zeta"), ux = f.read<double>("/particles/beam/ux");
      auto uy = f.read<double>("/particles/beam/uy"), uz = f.read<double>("/particles/beam/uz");
      auto q = f.read<double>("/particles/beam/charge");
      auto sp = f.read<int>("/particles/beam/species");
      std::printf("# x[um] y[um] zeta[um] ux uy uz q[pC] species\n");
      for (size_t k = 0; k < x.size(); ++k)
        std::printf("%.6e %.6e %.6e %.6e %.6e %.6e %.6e %d\n", x[k] * um, y[k] * um, z[k] * um, ux[k], uy[k], uz[k],
                    q[k] * pC, sp[k]);
      return 0;
    }
    const Data D(f);
    double e[6];
    if (what == "axis") {
      std::printf("# zeta[um]  Ex  Ez  ne/n0  Zni/n0\n");
      for (int n = 0; n < D.Nz; ++n) {
        D.fieldAt(n, 1e-9, 0.0, e);
        std::printf("%.6e %.6e %.6e %.6e %.6e\n", D.zs[n] * D.um, e[0], e[2], D.nE(n, 0), D.nI(n, 0));
      }
    } else if (what == "map") {
      const int stride = (argc > 3) ? std::atoi(argv[3]) : 4;
      std::printf("# zeta[um]  x[um]  Ex  Ez  ne/n0(m=0)  Zni/n0(m=0)\n");
      for (int n = 0; n < D.Nz; n += stride) {
        for (int sgn = -1; sgn <= 1; sgn += 2)
          for (int ii = 0; ii < D.N; ++ii) {
            const int i = (sgn < 0) ? D.N - 1 - ii : ii;
            const double x = sgn * D.c[i];
            D.fieldAt(n, x, 0.0, e);
            std::printf("%.5e %.5e %.5e %.5e %.5e %.5e\n", D.zs[n] * D.um, x * D.um, e[0], e[2], D.nE(n, i), D.nI(n, i));
          }
        std::printf("\n");
      }
    } else if (what == "energy") {
      double wt, wn, w0, w1;
      double zlo = -1e300, zhi = 1e300;  // optional zeta range [um]
      if (argc >= 5) { zlo = std::atof(argv[3]) / D.um; zhi = std::atof(argv[4]) / D.um; }
      D.energy(wt, wn, w0, w1, zlo, zhi);
      // SI: W [J] = eps0 E0^2 / k0^3 * W (code)
      const double E0 = f.attrD("E0_V_per_m"), L = f.attrD("k0inv_m");
      const double J = 8.8541878128e-12 * E0 * E0 * L * L * L;
      std::printf("# step  z[mm]  W_total  W_nyquist  W(m=0)  W(m>=1)  (code units; 1 = %.6e J)\n"
                  "%lld %.6g %.6e %.6e %.6e %.6e\n", J, f.attrI("step"), f.attrD("z_m") * 1e3, wt, wn, w0, w1);
    } else if (what == "spectrum") {
      if (argc < 5) throw std::runtime_error("spectrum needs <zeta_lo_um> <zeta_hi_um> [k_min k_max nk]");
      const double zlo = std::atof(argv[3]) / D.um, zhi = std::atof(argv[4]) / D.um;
      const double k1 = argc > 5 ? std::atof(argv[5]) : 0.5, k2 = argc > 6 ? std::atof(argv[6]) : 1.5;
      const int nk = argc > 7 ? std::atoi(argv[7]) : 201;
      std::printf("# k[k0]  lambda[um]  dW/dk (modes m >= 1, slices %s ... %s um; code units)\n", argv[3], argv[4]);
      for (int j = 0; j < nk; ++j) {
        const double k = k1 + (k2 - k1) * j / std::max(1, nk - 1);
        std::printf("%.6f %.6f %.6e\n", k, 2.0 * M_PI / k * D.um, D.spectrum(k, zlo, zhi));
      }
    } else if (what == "lineout") {
      if (argc < 5) throw std::runtime_error("lineout needs <r_um> <phi_deg>");
      const double r = std::atof(argv[3]) / D.um, ph = std::atof(argv[4]) * M_PI / 180.0;
      std::printf("# zeta[um]  Ex Ey Ez Bx By Bz  at r = %s um, phi = %s deg\n", argv[3], argv[4]);
      for (int n = 0; n < D.Nz; ++n) {
        D.fieldAt(n, r * std::cos(ph), r * std::sin(ph), e);
        std::printf("%.6e %.6e %.6e %.6e %.6e %.6e %.6e\n", D.zs[n] * D.um, e[0], e[1], e[2], e[3], e[4], e[5]);
      }
    } else {
      throw std::runtime_error("unknown command '" + what + "'");
    }
  } catch (const std::exception& ex) {
    std::fprintf(stderr, "h5_extract: %s\n", ex.what());
    return 1;
  }
  return 0;
}
