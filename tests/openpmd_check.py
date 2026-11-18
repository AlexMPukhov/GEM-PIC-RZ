"""Check an openPMD file of GEM-PIC-RZ against the native HDF5 file of the same step.
usage: python3 openpmd_check.py <output dir> <step> <h5_extract> [<tools dir>]"""
import subprocess
import sys

import h5py
import numpy as np

out, step, h5x = sys.argv[1], int(sys.argv[2]), sys.argv[3]
nat = f"{out}/data_{step:05d}.h5"
f = h5py.File(f"{out}/openpmd/openpmd_{step}.h5", "r")
fails = []


def need(ok, msg):
    if not ok:
        fails.append(msg)


a = f.attrs
need(a["openPMD"].decode() == "1.1.0", "openPMD version")
need(a["iterationEncoding"].decode() == "fileBased", "iterationEncoding")
it = f[f"/data/{step}"]
g = h5py.File(nat, "r")
E0, B0, s, L = g.attrs["E0_V_per_m"], g.attrs["B0_T"], g.attrs["z_m"], g.attrs["k0inv_m"]
need(abs(it.attrs["time"] * it.attrs["timeUnitSI"] - s / 299792458.0) < 1e-6 * s / 3e8 + 1e-30, "time = s / c")

E = it["meshes/E"]
need(E.attrs["geometry"].decode() == "thetaMode", "geometry")
need([x.decode() for x in E.attrs["axisLabels"]] == ["r", "z"], "axisLabels")
dr, dz = E.attrs["gridSpacing"] * E.attrs["gridUnitSI"]
z0 = E.attrs["gridGlobalOffset"][1] * E.attrs["gridUnitSI"]
need(abs(z0 - s) < 1e-12, "gridGlobalOffset z = s")
nm = E["r"].shape[0]


def field(rec, comp, r, th):  # at radius r [m], angle th, along z; SI
    d = it[f"meshes/{rec}/{comp}"]
    A = d[...] * d.attrs["unitSI"]
    j = int(r / dr)
    w = r / dr - j
    col = (1 - w) * A[:, j, :] + w * A[:, j + 1, :]
    F = col[0].copy()
    for m in range(1, (nm - 1) // 2 + 1):
        F += col[2 * m - 1] * np.cos(m * th) + col[2 * m] * np.sin(m * th)
    return F


zgrid = z0 + dz * np.arange(E["r"].shape[2])
for r_um, deg in [(2.0, 30.0), (5.0, 250.0)]:
    ref = np.loadtxt(subprocess.run([h5x, nat, "lineout", str(r_um), str(deg)], stdout=subprocess.PIPE,
                                    universal_newlines=True, check=True).stdout.splitlines())
    zeta = ref[:, 0]
    th = np.deg2rad(deg)
    Fr, Ft, Fz = (field(rec, c, r_um * 1e-6, th) for rec in ["E"] for c in "rtz")
    Ex = Fr * np.cos(th) - Ft * np.sin(th)
    Ey = Fr * np.sin(th) + Ft * np.cos(th)
    for name, F, col, U in [("Ex", Ex, 1, E0), ("Ey", Ey, 2, E0), ("Ez", Fz, 3, E0)]:
        fi = np.interp(zeta, (zgrid - s) / L * L * 1e6, F / U)
        err = np.max(np.abs(fi - ref[:, col])) / np.max(np.abs(ref[:, 1]))  # relative to the laser field
        print(f"  {name} at r = {r_um} um, phi = {deg} deg: max difference {err:.1e} of the peak field")
        need(err < 1e-3, f"{name} at ({r_um}, {deg})")

w = it["particles/witness"]
q = np.sum(w["weighting"][...]) * abs(w["charge"].attrs["value"]) * 1e12
print(f"  witness: {w['weighting'].shape[0]} macroparticles, {q:.4f} pC")
need(abs(q - 5.0) < 1e-6, "witness charge 5 pC")
zlab = w["position/z"][...] * w["position/z"].attrs["unitSI"] + w["positionOffset/z"].attrs["value"] * L
zn = g["particles/beam/zeta"][...] * L + s
need(abs(np.mean(zlab) - np.mean(zn)) < 1e-12, "lab z of the particles")

# position (z = const) snapshot, docs/openpmd_zplane.md: attributes and the crossing time of the particles
ia = it.attrs
need(ia["snapshotType"].decode() == "position", "snapshotType = position")
need(ia["snapshotAxis"].decode() == "z", "snapshotAxis = z")
need(abs(ia["snapshotPosition"] * ia["snapshotPositionUnitSI"] - s) <= 1e-12 * max(s, 1e-30), "snapshotPosition = s")
vax = ia["snapshotAxisVelocity"] * ia["snapshotPositionUnitSI"] / ia["timeUnitSI"]
need(abs(vax / 299792458.0 - 1) < 1e-9, "snapshotAxisVelocity = c")
t = w["time"][...] * w["time"].attrs["unitSI"]
# t = time + (snapshotPosition - x) / v for the axis label x = z (the same mapping as the meshes)
t_ref = ia["time"] * ia["timeUnitSI"] + (s - zlab) / vax
terr = np.max(np.abs(t - t_ref)) / (s / 299792458.0)
t_nat = (s - g["particles/beam/zeta"][...] * L) / 299792458.0
print(f"  particle time at the plane: max |t - (time + (s - z)/c)| = {terr:.1e} of s/c, "
      f"mean t = {np.mean(t):.6e} s (data file: {np.mean(t_nat):.6e} s)")
need(terr < 1e-9, "particle time = mapping of the axis label")
need(abs(np.mean(t) - np.mean(t_nat)) < 1e-9 * s / 3e8, "particle time = (s - zeta)/c of the data file")
# every species in every file, also the empty ones
spn = g["particles/beam/species"][...]
ions = [x for x in g["particles/beam"].attrs["ion_names"].decode().split(",") if x]
for name, code in [("witness", 0), ("trapped", -1)] + [("inj_" + x, -30 - k) for k, x in enumerate(ions)]:
    need(name in it["particles"], f"species {name} present")
    if name in it["particles"]:
        n = it[f"particles/{name}/position/x"].shape[0]
        need(n == np.count_nonzero(spn == code), f"species {name}: {n} particles")
        print(f"  species {name}: {n} macroparticles")
need("native_E_z" in it["meshes"], "native records with openpmd_grid = both")

# native records (rectilinear meshes) through the generic reader: exact native data
if len(sys.argv) > 4:
    sys.path.insert(0, sys.argv[4])
    from openpmd_rectilinear import read_mesh, theta_slice
    m = read_mesh(f"{out}/openpmd/openpmd_{step}.h5", "native_E_z")
    rc = g["grid/r_centres"][...] * L
    need(np.allclose(m.axes["r"], rc, rtol=1e-12), "native_E_z: radial nodes = cell centres")
    need(np.allclose(m.axes["z"] - s, g["grid/zeta_slices"][...] * L, rtol=0, atol=1e-18), "native_E_z: zeta nodes")
    th = np.deg2rad(30.0)
    Fz = theta_slice(m, th)
    Ez = g["fields/Ez"][...]  # [m][n][i] complex code units
    ref = Ez[0].real + 2 * sum((Ez[k] * np.exp(1j * k * th)).real for k in range(1, Ez.shape[0]))
    err = np.max(np.abs(Fz.T / E0 - ref)) / np.max(np.abs(ref))
    print(f"  native_E_z (rectilinear record) against the data file: max difference {err:.1e}")
    need(err < 1e-12, "native_E_z values")

if fails:
    print("FAILED:", ", ".join(fails))
    sys.exit(1)
