#!/usr/bin/env python3
"""Plasma and laser numbers for setting up a GEM-PIC-RZ deck, and starting values for its grid.

  python3 plasma_numbers.py --density 1e18 --wavelength 0.8 --a0 2 --waist 15 --tau 30 --length 2
  python3 plasma_numbers.py --density 1e17                       # beam-driven wake: no laser

density in cm^-3, wavelength and waist in um, tau = intensity FWHM in fs, length in mm. Python standard library only.
The grid lines are STARTING values that follow the rules of references/physics_rules.md; a production run still
needs a convergence check (halve dz, dzeta and dr in turn).
"""
import argparse
import math

p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument("--density", type=float, required=True, help="plasma electron density [cm^-3]")
p.add_argument("--wavelength", type=float, default=0.8, help="laser wavelength [um] (default 0.8)")
p.add_argument("--a0", type=float, help="normalized laser amplitude")
p.add_argument("--waist", type=float, help="laser waist w0, 1/e field radius [um]")
p.add_argument("--tau", type=float, help="laser duration, FWHM of the intensity [fs]")
p.add_argument("--length", type=float, help="propagation distance of the run = grid.z_end [mm]: sizes the radial box")
a = p.parse_args()

c = 299792458.0
lam = a.wavelength                                   # um
n_c = 1.11485e21 / lam**2                            # cm^-3
n = a.density
kp_inv = lam / (2 * math.pi) * math.sqrt(n_c / n)    # c / omega_p [um]
lam_p = 2 * math.pi * kp_inv
E0 = 96.0 * math.sqrt(n)                             # cold wave-breaking field [V/m]

print(f"critical density n_c        = {n_c:.4g} cm^-3   (lambda0 = {lam} um)")
print(f"n / n_c                     = {n / n_c:.4g}")
print(f"c / omega_p                 = {kp_inv:.4g} um")
print(f"plasma wavelength           = {lam_p:.4g} um  ({lam_p / c * 1e9:.4g} fs)")
print(f"wave-breaking field E0      = {E0 / 1e9:.4g} GV/m  (= {math.sqrt(n / n_c):.4g} m c omega0 / e, the unit of "
      f"the code's E_z)")
print(f"code units                  : length 1/k0 = {lam / (2 * math.pi):.5g} um, field m c omega0/e = "
      f"{3.2107e12 / lam:.4g} V/m")

laser = a.a0 is not None
if laser:
    a0 = a.a0
    w_match = 2 * math.sqrt(a0) * kp_inv
    w0 = a.waist if a.waist else w_match
    tau = a.tau if a.tau else lam_p / 3 / c * 1e9                # fs: c tau_fwhm = lambda_p / 3
    L1e = tau * 1e-15 * c * 1e6 / 1.1774                         # 1/e half length of the field [um]
    zR = math.pi * w0**2 / lam
    P = 21.5e-3 * (a0 * w0 / lam) ** 2                            # TW, linear polarization
    Pc = 17.4e-3 * n_c / n                                        # TW
    print()
    print(f"laser: a0 = {a0}, w0 = {w0:.4g} um" + ("" if a.waist else " (matched spot, not given)") +
          f", tau_fwhm = {tau:.4g} fs" + ("" if a.tau else " (c tau_fwhm = lambda_p / 3, not given)"))
    print(f"matched spot 2 sqrt(a0)/k_p = {w_match:.4g} um")
    print(f"Rayleigh length z_R         = {zR / 1e3:.4g} mm")
    print(f"power                       = {P:.4g} TW,  P / P_c = {P / Pc:.4g}  (relativistic self-focusing above 1)")
    print(f"energy (Gaussian)           = {P * 1e12 * tau * 1e-15 * 1.0645:.4g} J")
    print(f"length_1e = c tau_fwhm/1.177= {L1e:.4g} um   (1/e half length of the field)")
    print(f"dephasing length (linear)   = {lam_p**3 / lam**2 / 1e3:.4g} mm")
    print(f"wake estimate (1D, resonant): E_z ~ (a0^2/2)/sqrt(1 + a0^2/2) E0 = "
          f"{a0**2 / 2 / math.sqrt(1 + a0**2 / 2):.3g} E0 = "
          f"{a0**2 / 2 / math.sqrt(1 + a0**2 / 2) * math.sqrt(n / n_c):.3g} m c omega0/e (upper estimate)")
    dzeta = lam / 25
    nfront = math.ceil(3.05 * L1e / dzeta)                        # slices ahead of the laser centre
    nbehind = math.ceil(1.7 * lam_p / dzeta)
    nzeta = nfront + nbehind
    dr = min(w0 / 20, 0.1 * kp_inv)
    ramp = 3 * lam_p
    print()
    print("starting values for the deck (see references/physics_rules.md):")
    print(f"  grid.dzeta         = {dzeta * 1e3:.4g} nm       (lambda0 / 25)")
    print(f"  grid.nzeta         = {nzeta}          (window {nzeta * dzeta:.4g} um: 1.7 lambda_p behind the laser centre, "
          f"3 length_1e ahead)")
    print(f"  laser.zeta_center  = {nbehind * dzeta:.4g} um")
    print(f"  grid.dr (uniform)  = {dr:.3g} um      (min(w0 / 20, 0.1 c/omega_p); finer near the axis only for a thin "
          f"witness or trapped bunch)")
    print(f"  grid.dz            = {1.8 * kp_inv:.3g} um      (1.8 c/omega_p as in examples/lwfa_witness.deck; halve it "
          f"to check)")
    print(f"  plasma up-ramp     = {ramp:.3g} um       (about 3 lambda_p, linear), laser.z_focus at its end")
    print(f"  numerics.trap_gamma= {2 * (1 + a0**2 / 2):.3g}           (twice 1 + a0^2/2: quivering electrons are not "
          f"counted as trapped)")
    print(f"  plasma.smooth_length = {max(0.1 * kp_inv, dr):.3g} um    (for a0 >~ 2: 0.1 c/omega_p, not below the axis cell)")

    def box(w, why):
        print(f"  box ({why}): grid.r_absorber = {2.5 * w:.4g} um, grid.r_max = {3 * w:.4g} um, "
              f"plasma.r_max = {2.4 * w:.4g} um  ({math.ceil(3 * w / dr)} radial cells)")
    box(w0, "spot stays near w0: short run or guided pulse")
    if a.length:
        Lum = a.length * 1e3 - ramp
        wend = w0 * math.sqrt(1 + (max(Lum, 0) / zR) ** 2)
        print(f"  run length {a.length} mm = {a.length * 1e3 / zR:.2g} z_R: without guiding the spot grows to "
              f"{wend:.3g} um at the end (vacuum)")
        if wend > 1.3 * w0:
            box(wend, "unguided pulse, vacuum diffraction")
            print("  -> P/P_c near 1 does not guarantee guiding. Use the larger box, or the smaller one and check in the")
            print("     result that |E_x| at r_absorber stays below ~1 % of the peak (plot_run.py prints it).")
    else:
        print("  (give --length <mm> to size the box for the diffraction over the run)")
else:
    print()
    print("beam-driven wake, no laser (see references/physics_rules.md):")
    print(f"  grid.m_max = 0 for an axisymmetric (centred) driver")
    print(f"  grid.dzeta        <= {lam_p / 200:.3g} um     (lambda_p / 200, and sigma_z / 10 of the shortest bunch)")
    print(f"  grid.dr           <= {0.1 * kp_inv:.3g} um     (0.1 c/omega_p in the blowout, sigma_r / 10 at a thin bunch)")
    print(f"  grid.dz           : 1/15 of the shortest betatron length 2 pi sqrt(2 gamma) c/omega_p = "
          f"{2 * math.pi * math.sqrt(2.0) * kp_inv / 1e3:.3g} mm x sqrt(gamma)")
