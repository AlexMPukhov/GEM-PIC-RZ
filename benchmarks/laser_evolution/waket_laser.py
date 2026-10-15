"""Laser evolution benchmark: Wake-T 0.9.1 (laser envelope + quasi-static 2D plasma, Baxevanis-Stupakov model)
for the case of gem.deck and fbpic_laser.py.
The envelope a(xi, r) (A_x = Re[a exp(i k0 xi)]) is evolved by Wake-T's envelope solver with the plasma
susceptibility of the quasi-static model. Diagnostics every D_out (waket_out/laser_waket.txt):
  z[um]  W/W0 (field energy of the laser from |i k0 a + d_xi a|^2)  max|E_x| on axis [m c omega0/e]  <k>/k0  sigma_k/k0
and on-axis line-outs waket_out/axis_<z>.txt: zeta[um] (from the rear of the GEM-PIC-RZ window) Ez [m c omega0/e] |E_x| env.
usage: python3 waket_laser.py [z_end_um] [nxi] [nr] [dz_fields_um]"""
import os
import sys
import numpy as np
import scipy.constants as ct
from wake_t import GaussianPulse, PlasmaStage
from wake_t.utilities.bunch_generation import get_gaussian_bunch_from_size

lam = 0.8e-6
k0 = 2 * np.pi / lam
E0 = ct.m_e * ct.c ** 2 * k0 / ct.e
nc = ct.epsilon_0 * ct.m_e * (ct.c * k0) ** 2 / ct.e ** 2
a0, w0, L1e = 1.5, 28.0 / k0, 16.0 / k0
n0 = 0.01 * nc
z_ramp = 10e-6
Lw, zeta_c = 196.0 / k0, 146.0 / k0
r_box = 15e-6
z_end = float(sys.argv[1]) * 1e-6 if len(sys.argv) > 1 else 300e-6
nxi = int(sys.argv[2]) if len(sys.argv) > 2 else 490
nr = int(sys.argv[3]) if len(sys.argv) > 3 else 150
dz_f = float(sys.argv[4]) * 1e-6 if len(sys.argv) > 4 else 1e-6
D_out = 10e-6
tau_fwhm = L1e / ct.c * np.sqrt(2 * np.log(2))  # FWHM of the intensity for a field exp(-(xi/L)^2)

# window: xi = z - ct; laser centre at xi = 0, window [-zeta_c, Lw - zeta_c] as in GEM-PIC-RZ
xi_min, xi_max = -zeta_c, Lw - zeta_c
laser = GaussianPulse(xi_c=0.0, a_0=a0, w_0=w0, tau=tau_fwhm, z_foc=0.0, l_0=lam)
os.makedirs("waket_out", exist_ok=True)


def density(z):
    # the plasma of GEM-PIC-RZ at level s starts at the window front; the Wake-T box sees n(z) at the box,
    # (z: position of the box), quasi-static like GEM-PIC-RZ: n(z) for the whole window
    # (floor 1e-4 n0: Wake-T normalizes with the local density, which must not vanish)
    return n0 * np.clip(z / z_ramp, 1e-4, 1.0) if z_ramp > 0 else n0


rows = []
# a passive test electron far behind (Wake-T needs a bunch to track)
bunch = get_gaussian_bunch_from_size(1e-6, 1e-6, 1e-6, 1e-6, 1000.0, 0.001, 1e-15, xi_min + 2e-6, 1e-18, 10)
z = 0.0


def diag():
    a = laser.get_envelope()  # (nxi, nr) complex, xi increasing
    xi = np.linspace(xi_min, xi_max, a.shape[0])
    r = (np.arange(a.shape[1]) + 0.5) * r_box / a.shape[1]
    dxi = xi[1] - xi[0]
    # E_x = -d_t A_x = Re[(i k0 a + d_xi a) e^{i k0 xi}] (c = 1): field envelope
    da = np.gradient(a, dxi, axis=0)
    Eenv = np.abs(1j * k0 * a + da) / k0      # in m c omega0 / e
    W = np.sum(Eenv ** 2 * r[None, :]) * dxi * (r[1] - r[0])
    # spectrum: |FFT of (a e^{i k0 xi})|^2 along xi, weighted with r dr
    F = np.fft.fft(a * np.exp(1j * k0 * xi)[:, None], axis=0)
    k = 2 * np.pi * np.fft.fftfreq(xi.size, dxi)
    S = np.sum(np.abs(F) ** 2 * r[None, :], axis=1)
    ka = np.abs(k) / k0
    m = ka > 0.3
    kbar = np.sum(ka[m] * S[m]) / np.sum(S[m])
    ksig = np.sqrt(np.sum((ka[m] - kbar) ** 2 * S[m]) / np.sum(S[m]))
    return W, Eenv[:, 0].max(), kbar, ksig, xi, Eenv[:, 0]


n_out = int(round(z_end / D_out))
for j in range(n_out + 1):
    if j > 0:
        stage = PlasmaStage(D_out, lambda zz, z0=z: density(z0 + zz), wakefield_model="quasistatic_2d",
                            n_out=1, xi_min=xi_min, xi_max=xi_max, r_max=r_box, r_max_plasma=r_box,
                            n_r=nr, n_xi=nxi, ppc=4, dz_fields=dz_f, laser=laser, laser_evolution=True,
                            laser_envelope_substeps=4)
        stage.track(bunch, show_progress_bar=False)
        z += D_out
        wf = stage.wakefield
        Ez = wf.e_z[2:-2, 2] / E0 if hasattr(wf, "e_z") else None
    else:
        Ez = None
    if laser.a_env is None:  # before the first stage: initialize the envelope on the same grid
        laser.set_envelope_solver_params(xi_min, xi_max, r_box, nxi, nr, dz_f / ct.c, 4)
        laser.initialize_envelope()
    W, Exm, kbar, ksig, xi, Eax = diag()
    rows.append((z * 1e6, W, Exm, kbar, ksig))
    out = np.array(rows)
    out[:, 1] /= out[0, 1]
    np.savetxt("waket_out/laser_waket.txt", out, header="z[um]  W/W0  max|Ex|env_axis[m c w0/e]  <|k|>/k0  sigma_k/k0")
    cols = [(xi + zeta_c) * 1e6, Eax] + ([Ez] if Ez is not None and Ez.size == xi.size else [])
    np.savetxt(f"waket_out/axis_{z * 1e6:07.2f}.txt", np.c_[tuple(cols)],
               header="zeta[um] (from the rear of the GEM-PIC-RZ window)  |E_x| envelope on axis  [Ez]  [m c omega0/e]")
    print(f"z = {z * 1e6:7.2f} um  W/W0 = {out[-1, 1]:.4f}  max|Ex| = {Exm:.3f}  <k>/k0 = {kbar:.4f}", flush=True)
