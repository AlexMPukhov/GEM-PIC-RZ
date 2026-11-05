"""Stage 1 of the chained example (README.md): density down-ramp injection, simulated with FBPIC (full PIC).

GEM-PIC-RZ traps electrons in a down-ramp, but the trapped charge is not reliable there; injection is therefore
computed with a full PIC code, and GEM-PIC-RZ takes over the injected beam for the long acceleration (stage 2).
The case is that of benchmarks/downramp_fbpic (a0 = 1.8, w0 = 8.02 um, tau = 10.6 fs, lambda = 0.8 um; electron
density 0 -> 2 n0 over 0 ... 20 um, plateau to 60 um, down to n0 at 70 um, n0 = 3.48e24 m^-3 = 0.002 n_c), with the
parabolic guiding channel of stage 2: n = n0 (1 + 0.5 r^2 / w0^2), matched to the waist (Delta n(w0) = 1 / (pi r_e
w0^2) = 0.5 n0). Immobile ions, plasma column radius 19.1 um.

Output (at ct = z_end):
  stage1/hdf5/data<iteration>.h5   openPMD particles: the plasma electrons with u_z > 10 (species "electrons")
  stage1/laser_centre.txt          lab position of the laser centre (centroid of E_x^2 on the axis) and ct [m]
The hand-off script handoff.py turns both into the deck entries of stage 2 (position of the beam behind the laser).

usage: python3 stage1_fbpic.py [z_end_um] [--cpu]       (FBPIC 0.27, a CUDA GPU by default)"""
import os
import sys
import numpy as np
from scipy.constants import c, e, m_e
from fbpic.main import Simulation
from fbpic.lpa_utils.laser import add_laser_pulse, GaussianLaser
from fbpic.openpmd_diag import ParticleDiagnostic

lam, a0, w0, tau = 0.8e-6, 1.8, 8.0214e-6, 10.618e-15
n0 = 3.4838e24
zmin, zmax = -30.56e-6, 7.64e-6        # window around the initial laser centre (as the GEM-PIC-RZ window)
dz = lam / 24
Nz = int(round((zmax - zmin) / dz))
rmax, dr = 25e-6, 0.2e-6
Nr = int(round(rmax / dr))
dt = (zmax - zmin) / Nz / c
args = [a for a in sys.argv[1:] if not a.startswith("--")]
z_end = float(args[0]) * 1e-6 if args else 200e-6
use_cuda = "--cpu" not in sys.argv
os.makedirs("stage1", exist_ok=True)


def dens_func(z, r):
    f = np.interp(z, [0.0, 20e-6, 60e-6, 70e-6], [0.0, 2.0, 2.0, 1.0])  # constant n0 beyond 70 um
    return f * (1.0 + 0.5 * r ** 2 / w0 ** 2) * (r < 19.1e-6)


sim = Simulation(Nz, zmax, Nr, rmax, 2, dt, zmin=zmin, n_order=-1, use_cuda=use_cuda,
                 boundaries={"z": "open", "r": "open"}, verbose_level=1)
elec = sim.add_new_species(q=-e, m=m_e, n=n0, dens_func=dens_func, p_nz=2, p_nr=2, p_nt=4, p_zmin=zmin,
                           p_rmax=19.1e-6)
add_laser_pulse(sim, GaussianLaser(a0, w0, tau, z0=0.0, zf=0.0, lambda0=lam))
sim.set_moving_window(v=c)
N = int(round(z_end / (c * dt)))
sim.diags = [ParticleDiagnostic(N, {"electrons": elec}, select={"uz": [10.0, None]}, comm=sim.comm,
                                write_dir="stage1")]
sim.step(N + 1)  # the diagnostics of iteration N are written at the start of step N + 1

# laser centre: centroid of E_x^2 on the axis (mode 1: E_x = 2 Re E_r)
f1 = sim.fld.interp[1]
Ex = 2.0 * np.asarray(f1.Er[:, 0]).real
z = np.asarray(f1.z)
zl = float((z * Ex ** 2).sum() / (Ex ** 2).sum())
ct = c * N * dt
with open("stage1/laser_centre.txt", "w") as fo:
    fo.write("# lab position of the laser centre z_laser [m] and ct [m] at iteration %d\n%.9e %.9e\n" % (N, zl, ct))
print("stage 1 done: ct = %.1f um, laser centre at z = %.2f um (ct - z = %.2f um)" % (ct * 1e6, zl * 1e6,
                                                                                  (ct - zl) * 1e6))
