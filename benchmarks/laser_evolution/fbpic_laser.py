"""Laser evolution benchmark: FBPIC (full PIC, spectral, quasi-cylindrical, moving window) for the case of gem.deck.
Laser: a0, w0, field 1/e half length L (tau = L / c), lambda = 0.8 um, linear polarization (x), starting at focus.
Plasma: n0 behind a linear ramp of length z_ramp, immobile ions.
Window: the laser centre starts at z = 0; the window is [-zeta_c, L_w - zeta_c] as in GEM-PIC-RZ (zeta_c from the
rear of the window). GEM-PIC-RZ loads the plasma at the front of the window with n(s) (s = lab position of the
level), so in FBPIC the plasma profile is n(z - (L_w - zeta_c)).
Diagnostics every D_out (laser_fbpic.txt, one row per output):
  ct[um]  W(m=1) [J]  max|E_x| on axis [m c omega0/e]  <|k|>/k0 of the m = 1 fields  rms width of <|k|>
and the on-axis line-outs fbpic_axis_<ct>.txt: zeta[um] Ez Ex [m c omega0/e] n_e/n0.
usage: python3 fbpic_laser.py  [z_end_um] [dz_um]"""
import sys
import numpy as np
from scipy.constants import c, e, m_e, epsilon_0
from fbpic.main import Simulation
from fbpic.lpa_utils.laser import add_laser_pulse, GaussianLaser

lam = 0.8e-6
k0 = 2 * np.pi / lam
E0 = m_e * c * c * k0 / e                  # m c omega0 / e [V/m]
nc = epsilon_0 * m_e * (c * k0) ** 2 / e ** 2
# case (see gem.deck)
a0, w0, L1e = 1.5, 28.0 / k0, 16.0 / k0
n0 = 0.01 * nc
z_ramp = 10.0 * 1e-6
Lw, zeta_c = 196.0 / k0, 146.0 / k0
r_box = 15.0 * 1e-6
z_end = float(sys.argv[1]) * 1e-6 if len(sys.argv) > 1 else 300.0 * 1e-6
dz = float(sys.argv[2]) * 1e-6 if len(sys.argv) > 2 else lam / 24
D_out = 10.0 * 1e-6

zmin, zmax = -zeta_c, Lw - zeta_c
Nz = int(round((zmax - zmin) / dz))
dr = 0.125 * 1e-6
Nr = int(round(r_box / dr))
dt = (zmax - zmin) / Nz / c
front = zmax  # lab position of the window front at t = 0


def dens_func(z, r):
    return np.clip((z - front) / z_ramp, 0.0, 1.0) if z_ramp > 0 else np.where(z >= front, 1.0, 0.0)


sim = Simulation(Nz, zmax, Nr, r_box, 2, dt, zmin=zmin, n_order=-1, use_cuda=False,
                 boundaries={"z": "open", "r": "open"}, verbose_level=1)
sim.add_new_species(q=-e, m=m_e, n=n0, dens_func=dens_func, p_nz=2, p_nr=1, p_nt=4, p_zmin=zmin,
                    p_rmax=15.0 * 1e-6)
add_laser_pulse(sim, GaussianLaser(a0, w0, L1e / c, z0=0.0, zf=0.0, lambda0=lam))
sim.set_moving_window(v=c)


def diag():
    f0, f1 = sim.fld.interp[0], sim.fld.interp[1]
    r = f1.r
    dzg = f1.z[1] - f1.z[0]
    w = 2 * np.pi * r * (r[1] - r[0]) * dzg
    # FBPIC: F = F_0 + 2 Re sum_m F_m e^{-i m theta} -> <F^2>_theta of mode 1 = 2 |F_1|^2
    E2 = sum(np.abs(getattr(f1, a)) ** 2 for a in ("Er", "Et", "Ez"))
    B2 = sum(np.abs(getattr(f1, a)) ** 2 for a in ("Br", "Bt", "Bz"))
    W1 = 0.5 * epsilon_0 * np.sum(2 * (E2 + c * c * B2) * w[None, :])
    Ex_axis = 2.0 * f1.Er[:, 0].real / E0
    # spectrum of the m = 1 transverse fields along z, weighted with r dr
    k = 2 * np.pi * np.fft.fftfreq(f1.z.size, dzg)
    S = np.zeros_like(k)
    for a in ("Er", "Et", "Br", "Bt"):
        F = np.fft.fft(getattr(f1, a), axis=0)
        fac = 1.0 if a[0] == "E" else c * c
        S += fac * np.sum(np.abs(F) ** 2 * (r * (r[1] - r[0]))[None, :], axis=1)
    ka = np.abs(k) / k0
    m = ka > 0.3
    kbar = np.sum(ka[m] * S[m]) / np.sum(S[m])
    ksig = np.sqrt(np.sum((ka[m] - kbar) ** 2 * S[m]) / np.sum(S[m]))
    zeta = (f0.z - c * sim.time + zeta_c) * k0  # zeta from the rear of the window, k0^-1
    ne = -f0.rho[:, 0].real / e / n0
    return W1, np.abs(Ex_axis).max(), kbar, ksig, zeta, f0.Ez[:, 0].real / E0, Ex_axis, ne


rows = []
done = 0
n_out = int(round(z_end / D_out))
for j in range(n_out + 1):
    target = int(round(j * D_out / (c * dt)))
    if target > done:
        sim.step(target - done)
        done = target
    W1, Exm, kbar, ksig, zeta, Ez, Ex, ne = diag()
    ct = c * sim.time * 1e6
    rows.append((ct, W1, Exm, kbar, ksig))
    np.savetxt("fbpic_out/laser_fbpic.txt", np.array(rows),
               header="ct[um]  W(m=1)[J]  max|Ex|_axis[m c w0/e]  <|k|>/k0  sigma_k/k0")
    np.savetxt(f"fbpic_out/axis_{ct:07.2f}.txt", np.c_[zeta / k0 * 1e6, Ez, Ex, ne],
               header=f"FBPIC at ct = {ct:.3f} um: zeta[um] (from the rear of the window) Ez Ex [m c omega0/e] n_e/n0")
    print(f"ct = {ct:8.2f} um  W1 = {W1:.4e} J  max|Ex| = {Exm:.3f}  <k>/k0 = {kbar:.4f}", flush=True)
