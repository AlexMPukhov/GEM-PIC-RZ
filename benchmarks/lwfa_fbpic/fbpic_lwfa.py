"""LWFA benchmark: FBPIC (full PIC, spectral, moving window) for the same case as gem.deck.
a0 = 2, w0 = 8.0214 um, field 1/e half length 3.1831 um (tau = 10.618 fs), lambda = 0.8 um,
n0 = 3.4838e24 m^-3 (0.002 n_c) behind a 20 um linear ramp, immobile ions.
Coordinates: the laser centre starts at z = 0 (at focus); the window is [-30.56, 7.64] um as in
GEM-PIC-RZ (zeta_c = 30.56 um in a 38.2 um window). GEM-PIC-RZ loads the plasma at the front
of the window with the density n(s), so its lab profile is n(z - 7.64 um) relative to the laser
centre; FBPIC gets the same lab profile.
usage: python3 fbpic_lwfa.py [D_um ...]   (writes fbpic_D<D>.txt: zeta[um] Ez[V/m] Ex[V/m] n_e/n0 on axis)"""
import sys
import numpy as np
from scipy.constants import c, e, m_e
from fbpic.main import Simulation
from fbpic.lpa_utils.laser import add_laser_pulse, GaussianLaser

lam, a0, w0, tau = 0.8e-6, 2.0, 8.0214e-6, 10.618e-15
n0 = 3.4838e24
zmin, zmax = -30.56e-6, 7.64e-6
dz = float(sys.argv[2]) * 1e-6 if len(sys.argv) > 2 else lam / 32
dr = float(sys.argv[3]) * 1e-6 if len(sys.argv) > 3 else 0.2e-6
Nz = int(round((zmax - zmin) / dz))
rmax = 25e-6
Nr = int(round(rmax / dr))
Nm = 2
dt = (zmax - zmin) / Nz / c
D = [float(x) * 1e-6 for x in sys.argv[1].split(",")] if len(sys.argv) > 1 else [100e-6]


def dens_func(z, r):
    return np.clip((z - 7.64e-6) / 20e-6, 0.0, 1.0) * (r < 19.1e-6)


sim = Simulation(Nz, zmax, Nr, rmax, Nm, dt, zmin=zmin, n_order=-1, use_cuda=False,
                 boundaries={"z": "open", "r": "open"}, verbose_level=1)
sim.add_new_species(q=-e, m=m_e, n=n0, dens_func=dens_func, p_nz=2, p_nr=2, p_nt=4,
                    p_zmin=zmin, p_rmax=19.1e-6)
add_laser_pulse(sim, GaussianLaser(a0, w0, tau, z0=0.0, zf=0.0, lambda0=lam))
sim.set_moving_window(v=c)
done = 0
for d in D:
    n = int(round(d / (c * dt))) - done
    sim.step(n)
    done += n
    f0, f1 = sim.fld.interp[0], sim.fld.interp[1]
    z = f0.z - c * sim.time  # relative to the initial laser centre position moved with c: zeta - zeta_c
    Ez = f0.Ez[:, 0].real
    # E_x on the axis from mode 1: E_x = 2 Re(E_r,1) at theta = 0 (FBPIC: F = F_0 + 2 Re sum F_m e^{-i m theta}; checked: peak = a0 m c omega / e)
    Ex = 2.0 * f1.Er[:, 0].real
    ne = -f0.rho[:, 0].real / e / n0
    np.savetxt(f"fbpic_D{d * 1e6:.0f}.txt", np.c_[z * 1e6 + 30.56, Ez, Ex, ne],
               header=f"FBPIC at ct = {c * sim.time * 1e6:.3f} um: zeta[um] (window of GEM-PIC-RZ) Ez[V/m] Ex[V/m] n_e/n0 "
                      f"(dz = {dz * 1e6:.4f} um, dr = {dr * 1e6:.3f} um)")
    print("written", d, c * sim.time)
