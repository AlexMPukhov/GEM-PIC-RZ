#!/usr/bin/env python3
"""Deck generator: electron-driven dielectric wakefield structure (quartz-lined copper tube) with an external FODO
channel whose quadrupole strength follows the energy of the drive beam, optionally with a slice of a 400 GeV proton
bunch in the wake (proton microbuncher), a thin chicane and a second dielectric structure (proton "klystron").

Geometry and beam (A. Petrenko, I. Sheinman, RuPAC 2016, TUPSA041): vacuum channel Rc = 1 mm, quartz eps = 3.75 up to
Rw = 1.061 mm (61 um), copper wall at Rw; electron bunch Q = 5 nC, W0 = 600 MeV, sigma_z = 0.1 mm,
sigma_x = sigma_y = 0.05 mm, normalized emittance 1 mm mrad; FODO period 11 cm with 4.7 cm quadrupoles, pole field
1.5 T at the first quadrupole, later quadrupoles scaled with the energy of the bunch centre.

usage examples:
  python3 make_deck.py > drive.deck                                # electrons only, 5 m, linear energy taper guess
  python3 make_deck.py --energy-table run1/bunches.txt > drive2.deck   # taper from the energies of a previous run
  python3 make_deck.py --protons > microbuncher.deck               # + a slice of the proton bunch in the wake
  python3 make_deck.py --protons --klystron > klystron.deck        # + chicane + second structure
  python3 make_deck.py --help
All lengths in the deck carry units; GEM-PIC-RZ converts them (grid.lambda0 only fixes the code units).
"""
import argparse
import math
import sys

MC2_E = 0.51099895  # MeV
MC2_P = 938.27208816


def brho(E_kin_MeV, mc2):
    """magnetic rigidity [T m] for the kinetic energy"""
    pc = math.sqrt(E_kin_MeV * (E_kin_MeV + 2 * mc2))  # MeV
    return pc / 299.792458


def energy_profile(args):
    """kinetic energy of the bunch centre [MeV] as a function of z [m]: a table from bunches.txt of a previous run
    (columns: see the header of the file; the energy column is found by its name), or W0 - dWdz z"""
    if args.energy_table:
        # header: "# step  z[mm]  bunch  N  Q[pC]  E[MeV] ..."; E is the mean kinetic energy of the bunch
        z, E = [], []
        cols = None
        for line in open(args.energy_table):
            if line.startswith("#"):
                h = line.lstrip("#").split("(")[0].split()
                if "z[mm]" in h and "E[MeV]" in h:
                    cols = h
                continue
            v = line.split()
            if not v or cols is None or v[cols.index("bunch")] != args.bunch_name:
                continue
            z.append(float(v[cols.index("z[mm]")]) * 1e-3)
            E.append(float(v[cols.index("E[MeV]")]))
        if len(z) < 2:
            sys.exit("could not read z[mm] and E[MeV] from " + args.energy_table)

        def f(zz):
            if zz <= z[0]:
                return E[0]
            for i in range(1, len(z)):
                if zz <= z[i]:
                    t = (zz - z[i - 1]) / (z[i] - z[i - 1])
                    return E[i - 1] + t * (E[i] - E[i - 1])
            return E[-1]
        return f
    return lambda zz: max(args.W0 - args.dWdz * max(zz - args.z_structure, 0.0), 0.05 * args.W0)


def quadrupoles(args, f_E, z0, z1, first_sign=+1):
    """FODO between z0 and z1 [m]: F (length lq), drift, D, drift; gradient G(z) = G0 E(z_mid) / E(z0), signs
    alternating. G0 = B_pole / r_pole."""
    out = []
    G0 = args.B_pole / (args.r_pole * 1e-3)
    E0 = None  # normalization: the first quadrupole has G0 (pole field B_pole)
    half = 0.5 * args.period
    z = z0 + 0.5 * (half - args.lq)  # quadrupoles centred in the half cells
    sign = first_sign
    while z + args.lq <= z1 + 1e-12:
        zm = z + 0.5 * args.lq
        if E0 is None:
            E0 = f_E(zm)
        G = sign * G0 * (f_E(zm) / E0 if args.taper else 1.0)
        out.append("quadrupole { z_start = %.6f m; length = %.4f m; gradient = %.3f T/m }" % (z, args.lq, G))
        z += half
        sign = -sign
    return out, G0


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--length", type=float, default=5.0, help="length of the (first) structure [m] (5)")
    p.add_argument("--z-structure", dest="z_structure", type=float, default=0.01,
                   help="start of the dielectric [m]: the bunch starts in vacuum with its own field (0.01)")
    p.add_argument("--W0", type=float, default=600.0, help="electron kinetic energy [MeV] (600)")
    p.add_argument("--charge", type=float, default=5.0, help="electron charge [nC] (5)")
    p.add_argument("--dWdz", type=float, default=50.0,
                   help="guess of the energy loss of the bunch centre [MeV/m] for the taper (50); better: --energy-table")
    p.add_argument("--energy-table", help="bunches.txt of a previous run: taper from its energies")
    p.add_argument("--bunch-name", default="driver", help="bunch to read from the table (driver)")
    p.add_argument("--no-taper", dest="taper", action="store_false", help="constant quadrupole strength")
    p.add_argument("--period", type=float, default=0.11, help="FODO period [m] (0.11)")
    p.add_argument("--lq", type=float, default=0.047, help="quadrupole length [m] (0.047)")
    p.add_argument("--B-pole", dest="B_pole", type=float, default=1.5, help="pole field of the first quadrupole [T]")
    p.add_argument("--r-pole", dest="r_pole", type=float, default=1.5,
                   help="pole-tip radius [mm] (1.5: an assumption, G = B_pole/r_pole = 1000 T/m)")
    p.add_argument("--x0", type=float, default=10.0, help="transverse offset of the electron bunch [um] (seeds BBU)")
    p.add_argument("--window", type=float, default=10.0, help="window length [mm] (10; with protons: 30)")
    p.add_argument("--protons", action="store_true", help="add a slice of a 400 GeV proton bunch in the wake")
    p.add_argument("--p-charge", dest="p_charge", type=float, default=4.0,
                   help="charge of the proton slice [nC] (4: SPS 3e11 p, sigma_z 12 cm, near the centre, 25 mm)")
    p.add_argument("--p-length", dest="p_length", type=float, default=25.0, help="proton slice length [mm] (25)")
    p.add_argument("--klystron", action="store_true",
                   help="thin chicane after the first structure, then a second structure of the same length")
    p.add_argument("--r56", type=float, default=0.16,
                   help="chicane R56 [m] (0.16: ~quarter wavelength of 1 mm for delta ~ 1.6e-3 of 400 GeV protons)")
    p.add_argument("--gap", type=float, default=0.5, help="gap between the structures [m] (0.5, chicane in its middle)")
    p.add_argument("--p-quads", dest="p_quads", type=float, default=0.0,
                   help="second structure: FODO gradient for the protons [T/m] (0: none; same lattice geometry)")
    p.add_argument("--dz", type=float, default=2.0, help="field step [mm] (2)")
    p.add_argument("--every", type=int, default=250, help="output every N steps (250 = every 0.5 m)")
    p.add_argument("--backend", default="openmp", help="numerics.backend (openmp, cuda, serial)")
    args = p.parse_args()
    if args.protons and args.window < 30.0:
        args.window = 30.0

    f_E = energy_profile(args)
    L1 = args.z_structure + args.length
    z_end = L1 + 0.05
    if args.klystron:
        z_end = L1 + args.gap + args.length + 0.05
    d = []
    w = d.append
    w("# GEM-PIC-RZ deck generated by examples/dwa_proton_microbuncher/make_deck.py")
    w("#   arguments: " + (" ".join(sys.argv[1:]) or "(defaults)"))
    w("# Quartz-lined copper tube (Rc = 1 mm, eps = 3.75, Rw = 1.061 mm): synchronous TM01 at 996 um (301 GHz),")
    w("# TM02 194 um, TM03 100 um; the wake of 5 nC behind the bunch ~ 127 MV/m (modal theory, tests/dielectric.py).")
    w("# Electron drive beam 600 MeV, 5 nC, external FODO (period %.3f m, quadrupoles %.3f m, G0 = %.0f T/m)."
      % (args.period, args.lq, args.B_pole / (args.r_pole * 1e-3)))
    w("grid {")
    w("  lambda0 = 10 um;  m_max = 1                       # code units only; m = 1: dipole wakes (BBU / hosing)")
    w("  radial_grid = table;  r_max = 1.07 mm")
    w("  r_points  = 0, 0.3, 0.45, 0.85, 0.97 mm         # fine at the beam, coarser in the channel,")
    w("  dr_points = 5, 5, 15, 15, 3 um                   #   3 um in the 61 um quartz layer (k_perp dr << 1)")
    w("  dzeta = 4 um;  window_length = %.1f mm" % args.window)
    w("  dz = %.3f mm;  z_end = %.4f m" % (args.dz, z_end))
    w("}")
    w("wall { radius = 1.061 mm }                         # copper jacket")
    w("dielectric { epsilon = 3.75; r_in = 1 mm; r_out = 1.061 mm; z_start = %.4f m; z_end = %.4f m }"
      % (args.z_structure, L1))
    if args.klystron:
        z2 = L1 + args.gap
        w("dielectric { epsilon = 3.75; r_in = 1 mm; r_out = 1.061 mm; z_start = %.4f m; z_end = %.4f m }"
          % (z2, z2 + args.length))
    zc_e = args.window - 1.0
    w("bunch { name = driver; charge = %.3f nC; energy = %.3f MeV; energy_spread = 0.1 %%" % (args.charge, args.W0))
    w("        sigma_x = 50 um; sigma_y = 50 um; emittance_x = 1 um; emittance_y = 1 um   # normalized, 1 mm mrad")
    w("        sigma_z = 0.1 mm; zeta_center = %.3f mm; x0 = %.1f um" % (zc_e, args.x0))
    w("        macroparticles = 400000; azimuthal_copies = 8 }   # copies: m = 1 driven by the offset, not by noise")
    if args.protons:
        w("# slice of a 400 GeV SPS proton bunch (3e11 p, sigma_z = 12 cm: ~1e9 p/mm near the centre) behind the driver;")
        w("# its energy modulation by the wake is the microbunching")
        w("bunch { name = protons; species = proton; charge = %.3f nC; energy = 400 GeV; energy_spread = 0.03 %%"
          % args.p_charge)
        w("        sigma_x = 0.2 mm; sigma_y = 0.2 mm; emittance_x = 3 um; emittance_y = 3 um")
        w("        distribution_z = flattop; length = %.2f mm; zeta_center = %.3f mm" %
          (args.p_length, zc_e - 1.0 - 0.5 * args.p_length))
        w("        macroparticles = 1000000; azimuthal_copies = 8 }")
    quads, G0 = quadrupoles(args, f_E, args.z_structure, L1)
    w("# FODO channel of the first structure (gradients follow the energy of the bunch centre%s)"
      % ("" if args.taper else ": constant"))
    d.extend(quads)
    if args.klystron:
        zc = L1 + 0.5 * args.gap
        w("# thin chicane: protons with delta = gamma/gamma_ref - 1 move by r56 delta (the energy modulation becomes")
        w("# density modulation at the 1 mm wavelength); the electrons (delta ~ -1) leave the window")
        w("chicane { z = %.4f m; r56 = %.4f m; delay = 0 mm; energy = 400 GeV }" % (zc, args.r56))
        if args.p_quads > 0:
            class A: pass
            a2 = A(); a2.__dict__.update(vars(args)); a2.B_pole = args.p_quads * args.r_pole * 1e-3; a2.taper = False
            q2, _ = quadrupoles(a2, lambda z: 1.0, L1 + args.gap, L1 + args.gap + args.length)
            w("# FODO channel of the second structure for the protons (constant gradient)")
            d.extend(q2)
    nmin = max(1, int(round(args.dz * 1e-3 / 1e-3)))
    w("numerics { backend = %s; iterations = 2; beam_dz_max = 1 mm }   # beam push: 1 mm sub-steps (<= 0.1/k_beta)"
      % args.backend)
    w("output { dir = run; every = %d; particles = true }" % args.every)
    print("\n".join(d))
    sys.stderr.write("deck: %d quadrupoles, G0 = %.0f T/m, z_end = %.3f m, %d steps\n"
                     % (len(quads), G0, z_end, int(round(z_end / (args.dz * 1e-3)))))


if __name__ == "__main__":
    main()
