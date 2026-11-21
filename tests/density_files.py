#!/usr/bin/env python3
"""openPMD density meshes for tests/density_file.sh (plasma.density_file), written with h5py only.

  jet2d_theta.h5  f(z) g(r) = min(z/400, 1) (0.6 + 0.8 min(r/75, 1)) [lengths in 1/k0, lambda0 = 0.8 um]: a factor
                  on `density`, mesh 'density', geometry thetaMode ([mode][r][z], one mode), regular axes
  jet2d_rect.h5   the same as absolute density in m^-3 (n0 = 0.002 n_c), mesh 'n_gas' beside a second mesh, data
                  [z][r], rectilinear: the coordinates are data sets named in gridCoordinates
                  (docs/openpmd_rectilinear.md)
  gauss_jet.h5    smooth jet exp(-((z - 150)/200)^2 - (r/120)^2) on 81 x 41 points, a factor, thetaMode
"""
import numpy as np
import h5py

lam = 0.8e-6
k0inv = lam / (2 * np.pi)
nc_m3 = 1.1148e21 / 0.8**2 * 1e6


def root(name):
    f = h5py.File(name, "w")
    for k, v in dict(openPMD="1.1.0", basePath="/data/%T/", meshesPath="meshes/", iterationEncoding="fileBased",
                     iterationFormat=name[:-3] + "_%T", software="tests/density_files.py").items():
        f.attrs[k] = np.bytes_(v)
    f.attrs["openPMDextension"] = np.uint32(0)
    it = f.create_group("/data/0")
    it.attrs["time"] = 0.0
    it.attrs["dt"] = 1.0
    it.attrs["timeUnitSI"] = 1.0
    return f


def common(d, axes, unit_dim, unit_si, grid_unit):
    d.attrs["axisLabels"] = np.array([a.encode() for a in axes])
    d.attrs["dataOrder"] = np.bytes_("C")
    d.attrs["gridUnitSI"] = grid_unit
    d.attrs["unitSI"] = unit_si
    d.attrs["unitDimension"] = np.array(unit_dim, dtype=float)
    d.attrs["timeOffset"] = 0.0


def shape2d(z, r):  # [r][z]
    return (0.6 + 0.8 * np.minimum(r / 75.0, 1.0))[:, None] * np.minimum(z / 400.0, 1.0)[None, :]


# regular thetaMode mesh, lengths in 1/k0 (gridUnitSI), a factor on the density
z = 200.0 * np.arange(6)
r = 75.0 * np.arange(7)
f = root("jet2d_theta.h5")
d = f.create_dataset("/data/0/meshes/density", data=shape2d(z, r)[None, :, :])
d.attrs["geometry"] = np.bytes_("thetaMode")
d.attrs["geometryParameters"] = np.bytes_("m=0;imag=+")
common(d, ["r", "z"], [0] * 7, 1.0, k0inv)
d.attrs["gridSpacing"] = np.array([75.0, 200.0])
d.attrs["gridGlobalOffset"] = np.array([0.0, 0.0])
d.attrs["position"] = np.array([0.0, 0.0])
f.close()

# rectilinear mesh in SI: coordinates in m, density in m^-3, data [z][r]
z = np.array([0.0, 400.0, 1000.0])
r = np.array([0.0, 75.0, 150.0, 400.0])
f = root("jet2d_rect.h5")
f.create_dataset("/data/0/gridCoordinates/z_nozzle", data=z * k0inv)
f.create_dataset("/data/0/gridCoordinates/r_nozzle", data=r * k0inv)
d = f.create_dataset("/data/0/meshes/n_gas", data=(0.002 * nc_m3 * shape2d(z, r)).T.copy())
d.attrs["geometry"] = np.bytes_("other")
d.attrs["geometryParameters"] = np.bytes_("thetaMode;m=0;imag=+")
common(d, ["z", "r"], [-3, 0, 0, 0, 0, 0, 0], 1.0, 1.0)
d.attrs["gridCoordinates"] = np.array([b"/data/0/gridCoordinates/z_nozzle", b"/data/0/gridCoordinates/r_nozzle"])
d.attrs["gridSpacing"] = np.array([400.0 * k0inv, 75.0 * k0inv])  # informative: smallest spacings
d.attrs["gridGlobalOffset"] = np.array([0.0, 0.0])
d.attrs["position"] = np.array([0.0, 0.0])
t = f.create_dataset("/data/0/meshes/temperature", data=np.full((3, 4), 300.0))  # a second mesh: file_mesh needed
t.attrs["geometry"] = np.bytes_("other")
common(t, ["z", "r"], [0, 0, 0, 0, 1, 0, 0], 1.0, 1.0)
t.attrs["gridSpacing"] = np.array([1.0, 1.0])
t.attrs["gridGlobalOffset"] = np.array([0.0, 0.0])
t.attrs["position"] = np.array([0.0, 0.0])
f.close()

# smooth gas jet
z = np.linspace(0.0, 400.0, 81)
r = np.linspace(0.0, 160.0, 41)
f = root("gauss_jet.h5")
jet = np.exp(-(r[:, None] / 120.0) ** 2) * np.exp(-((z[None, :] - 150.0) / 200.0) ** 2)
d = f.create_dataset("/data/0/meshes/density", data=jet[None, :, :])
d.attrs["geometry"] = np.bytes_("thetaMode")
d.attrs["geometryParameters"] = np.bytes_("m=0;imag=+")
common(d, ["r", "z"], [0] * 7, 1.0, k0inv)
d.attrs["gridSpacing"] = np.array([4.0, 5.0])
d.attrs["gridGlobalOffset"] = np.array([0.0, 0.0])
d.attrs["position"] = np.array([0.0, 0.0])
f.close()
