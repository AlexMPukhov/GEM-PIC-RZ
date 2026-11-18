#!/usr/bin/env python3
"""Write a Gaussian laser pulse as an envelope file in the layout of LASY (openPMD mesh 'laserEnvelope',
thetaMode, axes t, r), for  laser { file = ... }  of GEM-PIC-RZ. Needs numpy and h5py only (not LASY).

  python3 make_lasy_gaussian.py                       # lasy_gaussian.h5: 0.8 um, w0 = 8 um, 23.5 fs FWHM
  python3 make_lasy_gaussian.py --waist 15 --tau 30 --out my_laser.h5

It is the stand-in for a real LASY file (measured profile, axiparabola, flying focus, ...): with LASY installed,
`Laser(...).write_to_file(...)` gives the same layout. The amplitude in the file is arbitrary here; set `a0` in the
laser block of the deck to rescale it.
"""
import argparse
import numpy as np
import h5py

p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument("--wavelength", type=float, default=0.8, help="um")
p.add_argument("--waist", type=float, default=8.0, help="w0, 1/e field radius at the plane of the file [um]")
p.add_argument("--tau", type=float, default=23.548, help="FWHM of the intensity [fs]")
p.add_argument("--out", default="lasy_gaussian.h5")
a = p.parse_args()

c = 299792458.0
lam = a.wavelength * 1e-6
omega = 2 * np.pi * c / lam
w0 = a.waist * 1e-6
tau = a.tau * 1e-15 / np.sqrt(2 * np.log(2))          # 1/e half duration of the field
t = np.arange(-3 * tau, 3 * tau + 1e-20, lam / c / 8)  # +- 3 tau, 8 points per period (the envelope is smooth)
r = np.arange(0, 5 * w0 + 1e-12, w0 / 40)              # out to 5 w0
env = 1e12 * np.exp(-(r[None, :] / w0) ** 2) * np.exp(-(t[:, None] / tau) ** 2)   # [V/m], linear along x

with h5py.File(a.out, "w") as f:
    for k, v in dict(openPMD="1.1.0", basePath="/data/%T/", meshesPath="meshes/", iterationEncoding="fileBased",
                     iterationFormat=a.out[:-3] + "_%05T", software="gempic-rz-simulation/make_lasy_gaussian.py").items():
        f.attrs[k] = np.bytes_(v)
    f.attrs["openPMDextension"] = np.uint32(0)
    d = f.create_dataset("/data/0/meshes/laserEnvelope", data=env[None].astype(np.complex128))  # [mode][t][r]
    d.attrs["geometry"] = np.bytes_("thetaMode")
    d.attrs["axisLabels"] = np.array([b"t", b"r"])
    d.attrs["dataOrder"] = np.bytes_("C")
    d.attrs["gridSpacing"] = np.array([t[1] - t[0], r[1] - r[0]])
    d.attrs["gridGlobalOffset"] = np.array([t[0], 0.0])
    d.attrs["gridUnitSI"] = 1.0
    d.attrs["unitSI"] = 1.0
    d.attrs["position"] = np.zeros(2)
    d.attrs["timeOffset"] = 0.0
    d.attrs["unitDimension"] = np.array([1.0, 1, -3, -1, 0, 0, 0])
    d.attrs["angularFrequency"] = omega
    d.attrs["polarization"] = np.array([1.0, 0.0])
    d.attrs["envelopeField"] = np.bytes_("electric_field")
print(f"{a.out}: {len(t)} x {len(r)} (t, r), lambda = {a.wavelength} um, w0 = {a.waist} um, tau_fwhm = {a.tau} fs")
