"""Density down-ramp injection benchmark: FBPIC (full PIC, spectral, quasi-cylindrical, moving window) for the case
of gem.deck: a0 = 1.8, w0 = 8.0214 um, field 1/e half length 3.1831 um (tau = 10.618 fs), lambda = 0.8 um;
electron density 0 -> 2 n0 over 0 ... 20 um, plateau to 60 um, down to n0 at 70 um, then n0 = 3.4838e24 m^-3
(0.002 n_c); immobile ions, plasma column radius 19.1 um.
Coordinates: window [-30.56, 7.64] um around the initial laser centre, as in GEM-PIC-RZ (laser centre 30.56 um from
the rear of a 38.2 um window). A level s of GEM-PIC-RZ is the plane z = s that the laser centre crosses at t = s/c
(GT2), so the lab profile n(z) is the profile of the deck without a shift.
Output every D_out: fbpic_out/elec_<ct>.txt (plasma electrons with u_z > 1: x y zeta[um] ux uy uz q[pC], zeta
from the rear of the window as in GEM-PIC-RZ) and fbpic_out/axis_<ct>.txt (zeta Ez Ex [m c omega0/e]).
usage: python3 fbpic_downramp.py [z_end_um]"""
import os
import sys
import numpy as np
from scipy.constants import c, e, m_e, m_p, epsilon_0
from fbpic.main import Simulation
from fbpic.lpa_utils.laser import add_laser_pulse, GaussianLaser

lam, a0, w0, tau = 0.8e-6, 1.8, 8.0214e-6, 10.618e-15
k0 = 2 * np.pi / lam
E0 = m_e * c * c * k0 / e
n0 = 3.4838e24
zmin, zmax = -30.56e-6, 7.64e-6
dz = lam / 24
Nz = int(round((zmax - zmin) / dz))
rmax, dr = 25e-6, 0.2e-6
Nr = int(round(rmax / dr))
dt = (zmax - zmin) / Nz / c
z_end = float(sys.argv[1]) * 1e-6 if len(sys.argv) > 1 else 200e-6
D_out = 25e-6
os.makedirs("fbpic_out", exist_ok=True)


def dens_func(z, r):
    f = np.interp(z, [0.0, 20e-6, 60e-6, 70e-6], [0.0, 2.0, 2.0, 1.0])  # (constant beyond 70 um)
    return f * (r < 19.1e-6)


sim = Simulation(Nz, zmax, Nr, rmax, 2, dt, zmin=zmin, n_order=-1, use_cuda=False,
                 boundaries={"z": "open", "r": "open"}, verbose_level=1)
elec = sim.add_new_species(q=-e, m=m_e, n=n0, dens_func=dens_func, p_nz=2, p_nr=2, p_nt=4, p_zmin=zmin,
                           p_rmax=19.1e-6)
add_laser_pulse(sim, GaussianLaser(a0, w0, tau, z0=0.0, zf=0.0, lambda0=lam))
sim.set_moving_window(v=c)
done = 0
for j in range(int(round(z_end / D_out)) + 1):
    target = int(round(j * D_out / (c * dt)))
    if target > done:
        sim.step(target - done)
        done = target
    ct = c * sim.time * 1e6
    zeta0 = c * sim.time - 30.56e-6  # lab position of the rear of the GEM-PIC-RZ window
    sp = elec
    sel = sp.uz[:sp.Ntot] > 1.0
    q = sp.w[:sp.Ntot][sel] * e * 1e12
    np.savetxt(f"fbpic_out/elec_{ct:07.2f}.txt",
               np.c_[sp.x[:sp.Ntot][sel] * 1e6, sp.y[:sp.Ntot][sel] * 1e6, (sp.z[:sp.Ntot][sel] - zeta0) * 1e6,
                     sp.ux[:sp.Ntot][sel], sp.uy[:sp.Ntot][sel], sp.uz[:sp.Ntot][sel], q],
               header=f"FBPIC ct = {ct:.3f} um: plasma electrons with u_z > 1: x y zeta[um] ux uy uz q[pC]")
    f0, f1 = sim.fld.interp[0], sim.fld.interp[1]
    np.savetxt(f"fbpic_out/axis_{ct:07.2f}.txt",
               np.c_[(f0.z - zeta0) * 1e6, f0.Ez[:, 0].real / E0, 2 * f1.Er[:, 0].real / E0],
               header=f"FBPIC ct = {ct:.3f} um: zeta[um] Ez Ex [m c omega0/e] on the axis")
    E = (np.sqrt(1 + sp.ux[:sp.Ntot][sel] ** 2 + sp.uy[:sp.Ntot][sel] ** 2 + sp.uz[:sp.Ntot][sel] ** 2) - 1) * 0.511
    print(f"ct = {ct:8.2f} um  u_z>1 {q.sum():.3f} pC  E>1 MeV {q[E > 1].sum():.3f} pC  E>5 MeV {q[E > 5].sum():.3f} pC",
          flush=True)
