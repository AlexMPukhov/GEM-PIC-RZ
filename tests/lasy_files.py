#!/usr/bin/env python3
"""Envelope files in the layout of LASY (openPMD, mesh 'laserEnvelope') for tests/lasy.sh, written with h5py only.

  gauss_rt.h5    Gaussian pulse at focus, w0 = 8 um, field 1/e half duration 20 fs, linear x, thetaMode
  gauss_circ.h5  the same, circular polarization (1, i)/sqrt(2)
  gauss_xyt.h5   the same as gauss_rt on a cartesian (t, y, x) grid
  ffm.h5         luminal Bessel pulse of GEM-PIC-RZ's flying_focus block (a0 = 2, 50 mrad, 50 fs), at z = 0
"""
import numpy as np
import h5py

lam = 0.8e-6
c = 299792458.0
w0SI = 2 * np.pi * c / lam
k0inv = lam / (2 * np.pi)
E0 = 0.51099895e6 / k0inv  # m c omega0 / e [V/m]


def write(name, env, axes, spacing, offset, pol, cartesian=False):
    f = h5py.File(name, "w")
    for k, v in dict(openPMD="1.1.0", basePath="/data/%T/", meshesPath="meshes/", iterationEncoding="fileBased",
                     iterationFormat=name[:-3] + "_%05T", software="tests/lasy_files.py").items():
        f.attrs[k] = np.bytes_(v)
    f.attrs["openPMDextension"] = np.uint32(0)
    d = f.create_dataset("/data/0/meshes/laserEnvelope", data=env.astype(np.complex128))
    d.attrs["geometry"] = np.bytes_("cartesian" if cartesian else "thetaMode")
    d.attrs["axisLabels"] = np.array([a.encode() for a in axes])
    d.attrs["dataOrder"] = np.bytes_("C")
    d.attrs["gridSpacing"] = np.array(spacing, dtype=float)
    d.attrs["gridGlobalOffset"] = np.array(offset, dtype=float)
    d.attrs["gridUnitSI"] = 1.0
    d.attrs["unitSI"] = 1.0
    d.attrs["position"] = np.zeros(len(axes))
    d.attrs["timeOffset"] = 0.0
    d.attrs["unitDimension"] = np.array([1.0, 1, -3, -1, 0, 0, 0])
    d.attrs["angularFrequency"] = w0SI
    d.attrs["polarization"] = np.array(pol)
    d.attrs["envelopeField"] = np.bytes_("electric_field")
    f.close()


# ---- Gaussian pulses
w0, tau = 8e-6, 20e-15
t = np.arange(-60e-15, 60e-15 + 1e-20, 0.3e-15)
r = np.arange(0, 40e-6 + 1e-12, 0.2e-6)
env = 1e12 * np.exp(-(r[None, :] / w0) ** 2) * np.exp(-(t[:, None] / tau) ** 2)
write("gauss_rt.h5", env[None], ["t", "r"], [t[1] - t[0], r[1] - r[0]], [t[0], 0.0], [1.0, 0.0])
write("gauss_circ.h5", env[None], ["t", "r"], [t[1] - t[0], r[1] - r[0]], [t[0], 0.0],
      np.array([1.0, 1j]) / np.sqrt(2))
x = np.arange(-30e-6, 30e-6 + 1e-12, 0.4e-6)
tc = t[::2]
envc = 1e12 * np.exp(-(x[None, :, None] ** 2 + x[None, None, :] ** 2) / w0**2) * np.exp(-(tc[:, None, None] / tau) ** 2)
write("gauss_xyt.h5", envc, ["t", "y", "x"], [tc[1] - tc[0], x[1] - x[0], x[1] - x[0]], [tc[0], x[0], x[0]],
      [1.0, 0.0], cartesian=True)


# ---- luminal Bessel pulse (Simulation::setupFlyingFocus): E_x = Re sum_j c_j J0(k_j r) e^{i w_j (zeta - zeta_c)}
def j0(xx):  # J0(x) = (1/pi) int_0^pi cos(x sin th) dth, trapezoid (exact to rounding for a periodic integrand)
    th = (np.arange(96) + 0.5) * np.pi / 96
    return np.cos(np.multiply.outer(xx, np.sin(th))).mean(axis=-1)


a0, theta = 2.0, 0.05
L = 50e-15 / np.sqrt(2 * np.log(2)) * w0SI
kap = 1 - np.cos(theta)
w = np.linspace(1 - 6 / L, 1 + 6 / L, 241)
kj = np.sqrt(np.maximum(0, 2 * w * kap - kap**2))
cj = a0 * (w[1] - w[0]) * L / (2 * np.sqrt(np.pi)) * np.exp(-0.25 * (w - 1) ** 2 * L**2)
rb = np.arange(0, 110e-6 + 1e-12, 0.5e-6)
tb = np.arange(-160e-15, 160e-15 + 1e-20, 1e-15)
J = j0(np.outer(rb / k0inv, kj))                         # [r][w]
P = np.exp(-1j * np.outer(tb * w0SI, w - 1))             # [t][w]: zeta - zeta_c = -t, envelope = G e^{+i t}
envb = (P * cj) @ J.T * E0
write("ffm.h5", envb[None], ["t", "r"], [tb[1] - tb[0], rb[1] - rb[0]], [tb[0], 0.0], [1.0, 0.0])
