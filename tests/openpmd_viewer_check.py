"""Read the openPMD files of a GEM-PIC-RZ run with openpmd-viewer and compare with the native data files.
usage: python3 openpmd_viewer_check.py <output dir> <h5_extract>   (see tests/openpmd_viewer.sh)"""
import glob
import subprocess
import sys

import h5py
import numpy as np
from openpmd_viewer import OpenPMDTimeSeries

out, h5x = sys.argv[1], sys.argv[2]
C = 299792458.0
fails = []


def need(ok, msg):
    if not ok:
        fails.append(msg)
        print("  FAIL:", msg)


ts = OpenPMDTimeSeries(f"{out}/openpmd", check_all_files=True)
its = list(ts.iterations)
ndata = len(glob.glob(f"{out}/data_*.h5"))
print(f"  openpmd-viewer: iterations {its}, fields {ts.avail_fields}, species {ts.avail_species}")
need(len(its) == ndata and len(its) >= 3, f"{len(its)} iterations, {ndata} data files")
need(not any(f.startswith("native_") for f in ts.avail_fields or []), "no native_* records by default")
g0 = h5py.File(f"{out}/data_{its[0]:05d}.h5", "r")
ions = [x for x in g0["particles/beam"].attrs["ion_names"].decode().split(",") if x]
codes = {"trapped": -1, **{"inj_" + x: -30 - k for k, x in enumerate(ions)}}
inj = "inj_" + ions[0]
need(sorted(ts.avail_species or []) == sorted(codes), f"species {sorted(codes)}")
counts = {sp: [] for sp in codes}
for it in its:
    g = h5py.File(f"{out}/data_{it:05d}.h5", "r")
    s, L, E0 = g.attrs["z_m"], g.attrs["k0inv_m"], g.attrs["E0_V_per_m"]
    f = h5py.File(f"{out}/openpmd/openpmd_{it}.h5", "r")
    ia = f[f"/data/{it}"].attrs
    # the species list of every file (openpmd-viewer takes it from the first one)
    need(sorted(f[f"/data/{it}/particles"].keys()) == sorted(codes), f"iteration {it}: species list")
    need(ia["snapshotType"].decode() == "position" and ia["snapshotAxis"].decode() == "z", f"iteration {it}: snapshot")
    need(abs(ia["snapshotPosition"] * ia["snapshotPositionUnitSI"] - s) <= 1e-12 * max(s, 1e-30) + 1e-30,
         f"iteration {it}: snapshotPosition = s")
    # fields at (r, phi) against the code's reconstruction from the native data
    for r_um, deg in [(1.0, 30.0), (3.0, 200.0)]:
        ref = np.loadtxt(subprocess.run([h5x, f"{out}/data_{it:05d}.h5", "lineout", str(r_um), str(deg)],
                                        stdout=subprocess.PIPE, universal_newlines=True, check=True).stdout.splitlines())
        th = np.deg2rad(deg)
        peak = np.max(np.abs(ref[:, 1]))
        for col, comp in [(1, "x"), (2, "y"), (3, "z")]:
            F, info = ts.get_field("E", comp, iteration=it, m="all", theta=th)
            rr = info.r
            pos = rr >= 0
            # value at radius r_um along z (rows: r from -rmax to rmax, columns: z)
            Fz = np.array([np.interp(r_um * 1e-6, rr[pos], F[pos, k]) for k in range(F.shape[1])])
            fi = np.interp(s + ref[:, 0] * 1e-6, info.z, Fz) / E0
            err = np.max(np.abs(fi - ref[:, col])) / peak
            if comp == "x" and r_um == 1.0:
                print(f"  iteration {it}: E_{comp} at r = {r_um} um, phi = {deg} deg: max difference {err:.1e} "
                      f"of the peak")
            need(err < 1e-3, f"iteration {it}: E_{comp} at ({r_um} um, {deg} deg): {err:.1e}")
    # particles: moments against the beam list of the data file
    sp_n = g["particles/beam/species"][...]
    zeta = g["particles/beam/zeta"][...] * L
    uz = g["particles/beam/uz"][...]
    q = g["particles/beam/charge"][...] * g["particles/beam/charge"].attrs["unitSI"]
    for sp, code in codes.items():
        sel = sp_n == code
        z, uzv, w = ts.get_particle(["z", "uz", "w"], species=sp, iteration=it)
        t = f[f"/data/{it}/particles/{sp}/time"]
        t = t[...] * t.attrs["unitSI"]
        n = len(w)
        counts[sp].append(n)
        need(n == np.count_nonzero(sel) and len(z) == n and len(t) == n, f"iteration {it}, {sp}: {n} particles")
        if n == 0 or n != np.count_nonzero(sel):
            continue
        Q = np.sum(w) * 1.602176634e-19
        Qn = np.sum(np.abs(q[sel]))
        need(abs(Q - Qn) <= 1e-6 * Qn, f"iteration {it}, {sp}: charge {Q:.6e} vs {Qn:.6e} C")
        zn = np.average(zeta[sel] + s, weights=np.abs(q[sel]))
        zm = np.average(z, weights=w)  # openpmd-viewer 1.11: z in m (older versions: um)
        if abs(zm * 1e-6 - zn) < abs(zm - zn):
            zm *= 1e-6
        need(abs(zm - zn) <= 1e-6 * (abs(zn) + 1e-6), f"iteration {it}, {sp}: mean z {zm:.9e} vs {zn:.9e} m")
        um = np.average(uzv, weights=w)
        un = np.average(uz[sel], weights=np.abs(q[sel]))
        need(abs(um - un) <= 1e-5 * (abs(un) + 1e-3), f"iteration {it}, {sp}: mean uz {um:.6e} vs {un:.6e}")
        terr = np.max(np.abs(t - (s - zeta[sel]) / C)) / (s / C)
        need(terr < 1e-9, f"iteration {it}, {sp}: time = (s - zeta)/c ({terr:.1e})")
        print(f"  iteration {it}, {sp}: {n} macroparticles, {Q * 1e12:.4e} pC, <z> = {zm * 1e6:.4f} um, "
              f"<u_z> = {um:.5f}; data file: {Qn * 1e12:.4e} pC, {zn * 1e6:.4f} um, {un:.5f}")
print(f"  macroparticles per iteration: {counts}")
need(counts[inj][0] == 0 and max(counts[inj]) > 0, f"{inj} empty at first, filled later")

if fails:
    print("FAILED:", "; ".join(fails))
    sys.exit(1)
print("  openpmd-viewer: all iterations, fields and particles = native data")
