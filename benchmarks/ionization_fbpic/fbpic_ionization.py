"""Ionization injection benchmark: FBPIC (full PIC, spectral, quasi-cylindrical, moving window, ADK ionization) for
the case of gem.deck: a0 = 2.3, w0 = 8.0214 um, field 1/e half length 3.1831 um (tau = 10.618 fs), lambda = 0.8 um,
free electrons n0 = 3.4838e24 m^-3 (0.002 n_c) behind a 20 um linear ramp (immobile ion background), and N5+ ions
(1 % of n0, mobile) ionized to N6+, N7+ by ADK; their electrons are a separate species.
Coordinates: window [-30.56, 7.64] um around the initial laser centre, as in GEM-PIC-RZ (laser centre 30.56 um from
the rear of a 38.2 um window). A level s of GEM-PIC-RZ is the plane z = s that the laser centre crosses at t = s/c
(GT2), so the lab profile n(z) is the profile of the deck without a shift.
Output every D_out: fbpic_out/elec_N_<ct>.txt (all electrons from nitrogen: x y zeta[um] ux uy uz q[pC], zeta from
the rear of the window as in GEM-PIC-RZ), fbpic_out/axis_<ct>.txt (zeta Ez Ex [m c omega0/e]) and
fbpic_out/ionization.txt (ct, fractions of N5+, N6+, N7+ over all nitrogen macroparticles).
usage: python3 fbpic_ionization.py [z_end_um]"""
import os
import sys
import numpy as np
from scipy.constants import c, e, m_e, m_p, epsilon_0
from fbpic.main import Simulation
from fbpic.lpa_utils.laser import add_laser_pulse, GaussianLaser

lam, a0, w0, tau = 0.8e-6, 2.3, 8.0214e-6, 10.618e-15
k0 = 2 * np.pi / lam
E0 = m_e * c * c * k0 / e
n0 = 3.4838e24
nN = 0.01 * n0
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
    return np.clip(z / 20e-6, 0.0, 1.0) * (r < 19.1e-6)


sim = Simulation(Nz, zmax, Nr, rmax, 2, dt, zmin=zmin, n_order=-1, use_cuda=False,
                 boundaries={"z": "open", "r": "open"}, verbose_level=1)
sim.add_new_species(q=-e, m=m_e, n=n0, dens_func=dens_func, p_nz=2, p_nr=1, p_nt=4, p_zmin=zmin, p_rmax=19.1e-6)
ions = sim.add_new_species(q=5 * e, m=14.007 * 1.66053907e-27, n=nN, dens_func=dens_func, p_nz=2, p_nr=2, p_nt=4,
                           p_zmin=zmin, p_rmax=19.1e-6)
elec_N = sim.add_new_species(q=-e, m=m_e)
ions.make_ionizable("N", target_species=elec_N, level_start=5)
add_laser_pulse(sim, GaussianLaser(a0, w0, tau, z0=0.0, zf=0.0, lambda0=lam))
sim.set_moving_window(v=c)
rows = []
done = 0
for j in range(int(round(z_end / D_out)) + 1):
    target = int(round(j * D_out / (c * dt)))
    if target > done:
        sim.step(target - done)
        done = target
    ct = c * sim.time * 1e6
    zeta0 = c * sim.time - 30.56e-6  # lab position of the rear of the GEM-PIC-RZ window
    sp = elec_N
    w = sp.w[:sp.Ntot]
    q = w * e * 1e12
    np.savetxt(f"fbpic_out/elec_N_{ct:07.2f}.txt",
               np.c_[sp.x[:sp.Ntot] * 1e6, sp.y[:sp.Ntot] * 1e6, (sp.z[:sp.Ntot] - zeta0) * 1e6,
                     sp.ux[:sp.Ntot], sp.uy[:sp.Ntot], sp.uz[:sp.Ntot], q],
               header=f"FBPIC ct = {ct:.3f} um: electrons from N: x y zeta[um] ux uy uz q[pC]")
    f0, f1 = sim.fld.interp[0], sim.fld.interp[1]
    np.savetxt(f"fbpic_out/axis_{ct:07.2f}.txt",
               np.c_[(f0.z - zeta0) * 1e6, f0.Ez[:, 0].real / E0, 2 * f1.Er[:, 0].real / E0],
               header=f"FBPIC ct = {ct:.3f} um: zeta[um] Ez Ex [m c omega0/e] on the axis")
    lev = ions.ionizer.ionization_level[:ions.Ntot]
    wi = ions.w[:ions.Ntot]
    frac = [np.sum(wi[lev == k]) / np.sum(wi) for k in (5, 6, 7)]
    rows.append([ct] + frac)
    np.savetxt("fbpic_out/ionization.txt", np.array(rows), header="ct[um] f(N5+) f(N6+) f(N7+) (all macroparticles)")
    E = (np.sqrt(1 + sp.ux[:sp.Ntot] ** 2 + sp.uy[:sp.Ntot] ** 2 + sp.uz[:sp.Ntot] ** 2) - 1) * 0.511
    print(f"ct = {ct:8.2f} um  freed {q.sum():.3f} pC  E>1 MeV {q[E > 1].sum():.3f} pC  E>5 MeV {q[E > 5].sum():.3f} pC",
          flush=True)
