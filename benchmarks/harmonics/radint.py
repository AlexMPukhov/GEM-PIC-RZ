#!/usr/bin/env python3
"""Reference for harmonic emission in a planar undulator: coherent emission of a short Gaussian bunch from the
far-field radiation integral (Jackson 14.65), as tools/undulator_radint.py, for the trajectory in the field of
GEM-PIC-RZ (B_y = B0 g(z) cos(k_u z), linear ramps of whole periods, no natural focusing), times N_e^2 and the
form factor exp(-k^2 sigma_z^2 - k^2 theta^2 sigma_r^2).

The far-field amplitude is split by its parity under phi -> phi + pi: the part with A(phi + pi) = A(phi) is the
radiation of the odd azimuthal modes (m = 1, 3, ...: fundamental and odd harmonics), the part with
A(phi + pi) = -A(phi) that of the even modes (m = 0, 2, ...: the even harmonics, zero on the axis).
Writes radint_<tag>.txt: k/k0, dW/dk (all, odd m, even m) in code units, and prints the energies.

usage: python3 radint.py gamma K period_m periods sigma_z_m sigma_r_m charge_C lambda0_m tag"""
import sys
import numpy as np

gamma, K, period, periods = float(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4])
sigma_z, sigma_r, charge, lam0, tag = float(sys.argv[5]), float(sys.argv[6]), float(sys.argv[7]), float(sys.argv[8]), sys.argv[9]
ramp_periods = 1
c, e, eps0, me = 299792458.0, 1.602176634e-19, 8.8541878128e-12, 9.1093837015e-31
k0 = 2 * np.pi / lam0
lu = period * k0; ku = 2 * np.pi / lu; L = periods * lu; ramp = ramp_periods * lu; B0 = K * ku
sz, sr, Ne = sigma_z * k0, sigma_r * k0, charge / e
nu_r = 2 * gamma**2 / (1 + K * K / 2) * ku          # k_r / k0
nu = np.linspace(0.3, 3.6, 1321) * nu_r
theta = np.linspace(0.0, 2.5 / gamma, 126)
nphi = 32
phi = np.arange(nphi) * 2 * np.pi / nphi

z = np.linspace(0.0, L, periods * 96 + 1)
g = np.minimum(1.0, np.minimum(z / ramp, (L - z) / ramp))
By = B0 * g * np.cos(ku * z)
cum = lambda f: np.concatenate([[0.0], np.cumsum(0.5 * (f[1:] + f[:-1]) * np.diff(z))])
ux = cum(By); uz = np.sqrt(gamma**2 - 1 - ux**2)
bx, bz, dtdz = ux / gamma, uz / gamma, gamma / uz
t, x = cum(dtdz), cum(bx * dtdz)
dbx, dbz = np.gradient(bx, z) / dtdz, np.gradient(bz, z) / dtdz
wt = dtdz * np.gradient(z)
s_all, s_odd, s_even = (np.zeros(len(nu)) for _ in range(3))
dth = theta[1] - theta[0]
per_th = []
per_m = np.zeros((nphi // 2 + 1, len(nu)))
for th in theta:
    amp = []
    for p in phi:
        n = np.array([np.sin(th) * np.cos(p), np.sin(th) * np.sin(p), np.cos(th)])
        a0 = n[1] * dbz; a1 = (n[2] - bz) * dbx - (n[0] - bx) * dbz; a2 = -n[1] * dbx
        A = [n[1] * a2 - n[2] * a1, n[2] * a0 - n[0] * a2, n[0] * a1 - n[1] * a0]
        w = wt / (1 - n[0] * bx - n[2] * bz) ** 2
        E = np.exp(1j * np.outer(nu, t - (n[0] * x + n[2] * z)))
        amp.append(np.array([E @ (Ai * w) for Ai in A]))       # (3, nu): Cartesian vector amplitude
    amp = np.array(amp)                                         # (phi, 3, nu)
    F2 = np.exp(-nu**2 * (sz**2 + (np.sin(th) * sr) ** 2))
    wgt = F2 * np.sin(th) * dth * (2 * np.pi / nphi)
    sh = np.roll(amp, -nphi // 2, axis=0)                       # A(phi + pi)
    s_all += wgt * (np.abs(amp) ** 2).sum(axis=(0, 1))
    s_odd += wgt * (np.abs(0.5 * (amp + sh)) ** 2).sum(axis=(0, 1))
    s_even += wgt * (np.abs(0.5 * (amp - sh)) ** 2).sum(axis=(0, 1))
    cp, sp = np.cos(phi)[:, None], np.sin(phi)[:, None]
    comps = [amp[:, 0] * cp + amp[:, 1] * sp, -amp[:, 0] * sp + amp[:, 1] * cp, amp[:, 2]]
    for f in comps:                                             # azimuthal Fourier modes (Parseval)
        cm = np.abs(np.fft.fft(f, axis=0) / nphi) ** 2 * nphi    # sum over phi of |f|^2 = nphi sum_m |c_m|^2
        for m in range(nphi // 2 + 1):
            per_m[m] += wgt * (cm[m] + (cm[-m] if 0 < m < nphi // 2 else 0))
    per_th.append([th, (wgt * (np.abs(0.5 * (amp + sh)) ** 2).sum(axis=(0, 1))).sum(), (wgt * (np.abs(0.5 * (amp - sh)) ** 2).sum(axis=(0, 1))).sum()])
unit = e**2 * (c * k0) / (16 * np.pi**3 * eps0 * c)        # J per (d nu dOmega)
E0 = me * c**2 * k0 / e; Wcode = eps0 * E0**2 / k0**3       # J per code unit of energy
dnu = nu[1] - nu[0]
W = [Ne**2 * unit * s.sum() * dnu / Wcode for s in (s_all, s_odd, s_even)]
print(f"{tag}: k_r = {nu_r:.5f} k0; coherent energy (code units) all {W[0]:.4e}  odd m {W[1]:.4e}  even m {W[2]:.4e}"
      f"  (even/odd {W[2] / W[1]:.4e})")
for h in (1, 2, 3):
    b = (nu > (h - 0.5) * nu_r) & (nu < (h + 0.5) * nu_r)
    print(f"  band {h - 0.5:.1f}..{h + 0.5:.1f} k_r: odd {Ne**2 * unit * s_odd[b].sum() * dnu / Wcode:.4e}"
          f"  even {Ne**2 * unit * s_even[b].sum() * dnu / Wcode:.4e}")
np.savetxt(f"radint_{tag}.txt", np.c_[nu, Ne**2 * unit * np.c_[s_all, s_odd, s_even] / Wcode],
           header="k/k0  dW/dk all, odd m, even m (code units)")
pt = np.array(per_th); pt[:, 1:] *= Ne**2 * unit * dnu / Wcode
pm = Ne**2 * unit * per_m.sum(axis=1) * dnu / Wcode
print("  per azimuthal mode |m|: " + " ".join(f"{m}:{v:.3e}" for m, v in enumerate(pm[:8])))
np.savetxt(f"radint_{tag}_modes.txt", np.c_[nu, Ne**2 * unit * per_m[:8].T / Wcode], header="k/k0  dW/dk of |m| = 0..7")
np.savetxt(f"radint_{tag}_theta.txt", pt, header="theta [rad]  energy in this theta bin: odd m, even m (code units)")
