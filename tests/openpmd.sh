#!/bin/sh
# openPMD output (output.format = openpmd): tests/regression.deck for 2 steps
#   1. the files pass the openPMD validator (openPMD_check_h5, if installed)
#   2. tests/openpmd_check.py (python3 with numpy, h5py): the standard's attributes,
#      the fields rebuilt from the thetaMode data at (r, theta) (output nodes on the slices) equal the code's own
#      reconstruction (h5_extract lineout of data_NNNNN.h5) to 1e-3 of the peak,
#      lab z = zeta + s, beam charge and particle count; the native records (rectilinear meshes
#      native_*) read with tools/openpmd_rectilinear.py equal the native data
#   3. h5_extract gives the same tables from the native records of the openPMD file as from
#      data_NNNNN.h5 (axis, map, lineout, energy, spectrum: identical text)
#   4. without particles the file has no particlesPath (validator)
#   5. density normalization (uniform plasma): n_e, Zn_i in SI = deck density ahead of the laser, = the data file
# Exit code 77 (skipped) without python3 / numpy / h5py.
# usage: sh tests/openpmd.sh <gempic_rz> <h5_extract>
set -e
abspath() { (cd "$(dirname "$1")" && echo "$(pwd)/$(basename "$1")"); }
EXE=$(abspath "$1"); H5X=$(abspath "$2")
SRC=$(cd "$(dirname "$0")/.." && pwd)
python3 -c "import numpy, h5py" 2>/dev/null || { echo "python3 with numpy and h5py not found: skipped"; exit 77; }
WORK=$(mktemp -d)
cd "$WORK"
"$EXE" "$SRC/tests/regression.deck" numerics.backend=serial grid.z_end=80k0^-1 output.every=2 \
  output.format=hdf5,openpmd output.openpmd_grid=both output.openpmd_dr=0.1um output.openpmd_dzeta=0.25k0^-1 \
  output.dir=out > run.log
test -f out/openpmd/openpmd_0.h5 && test -f out/openpmd/openpmd_2.h5
if command -v openPMD_check_h5 > /dev/null 2>&1; then
  for f in out/openpmd/openpmd_0.h5 out/openpmd/openpmd_2.h5; do
    openPMD_check_h5 -i "$f" > check.log 2>&1 || true
    grep -q "Result: 0 Errors" check.log || { cat check.log; echo "openPMD validator: errors in $f"; exit 1; }
  done
  echo "  openPMD validator: 0 errors"
fi
python3 "$SRC/tests/openpmd_check.py" out 2 "$H5X" "$SRC/tools"
for cmd in axis "map 8" "lineout 2 30" energy "spectrum 0 30 0.5 1.5 21"; do
  "$H5X" out/data_00002.h5 $cmd > a.txt
  "$H5X" out/openpmd/openpmd_2.h5 $cmd > b.txt
  cmp -s a.txt b.txt || { echo "h5_extract $cmd: openPMD native records differ from data_00002.h5"; exit 1; }
done
echo "  h5_extract: native openPMD records = data file (axis, map, lineout, energy, spectrum)"
"$EXE" "$SRC/tests/regression.deck" numerics.backend=serial grid.z_end=40k0^-1 output.every=1 \
  output.format=openpmd output.openpmd_grid=native output.particles=false output.dir=np > np.log
if command -v openPMD_check_h5 > /dev/null 2>&1; then
  openPMD_check_h5 -i np/openpmd/openpmd_1.h5 > check.log 2>&1 || true
  grep -q "Result: 0 Errors" check.log || { cat check.log; echo "openPMD validator: errors without particles"; exit 1; }
fi
# density normalization: uniform plasma 0.002 n_c (no ramp); n_e and Zn_i of the uniform and native records in SI
# against the deck density (n_c = eps0 m omega^2 / e^2 from lambda0 = 0.8 um) ahead of the laser, and against the
# densities of the native data file everywhere
"$EXE" "$SRC/tests/regression.deck" numerics.backend=serial grid.z_end=40k0^-1 output.every=1 plasma.z_values=1,1 \
  output.format=hdf5,openpmd output.openpmd_grid=both output.particles=false output.dir=dn > dn.log
python3 - dn 1 <<'PY'
import sys
import h5py
import numpy as np
out, step = sys.argv[1], int(sys.argv[2])
e, me, eps0, c = 1.602176634e-19, 9.1093837015e-31, 8.8541878128e-12, 299792458.0
n0 = 0.002 * eps0 * me * (2 * np.pi * c / 0.8e-6) ** 2 / e ** 2   # deck density [m^-3]
f = h5py.File(f"{out}/openpmd/openpmd_{step}.h5", "r")
g = h5py.File(f"{out}/data_{step:05d}.h5", "r")
L = g.attrs["k0inv_m"]
it = f[f"/data/{step}/meshes"]
ok = True
for rec, nat in [("n_e", "ne_over_n0"), ("Zn_i", "Zni_over_n0")]:
    d = it[rec]
    A = d[0] * d.attrs["unitSI"]                       # [r][z] in m^-3
    dr, dz = d.attrs["gridSpacing"] * d.attrs["gridUnitSI"]
    r = dr * np.arange(A.shape[0])
    zeta = dz * np.arange(A.shape[1])
    win = zeta[-1]
    # ahead of the laser (zeta > 270 k0^-1, laser at 240 +- 25), inside the plasma (r < 140 k0^-1)
    sel = (r > 20 * L) & (r < 140 * L)
    head = (zeta > 270 * L) & (zeta < win - 2 * L)
    v = A[np.ix_(sel, head)]
    err = abs(np.mean(v) / n0 - 1)
    print(f"  {rec} (uniform record) ahead of the laser: {np.mean(v):.5e} m^-3, deck {n0:.5e} m^-3: {err:.1e}")
    ok &= err < 0.01
    m = it["native_" + rec]
    B = m[0] * m.attrs["unitSI"]
    ref = (g[f"densities/{nat}"][...] * g[f"densities/{nat}"].attrs["unitSI"]).T   # [i][k]
    e2 = np.max(np.abs(B - ref)) / n0
    print(f"  native_{rec} = data file: max difference {e2:.1e} of n0")
    ok &= e2 < 1e-6
    # uniform vs native at the native nodes (linear interpolation): same scale
    rc = g["grid/r_centres"][...] * L
    zh = g["grid/zeta_halfcells"][...] * L
    i = np.searchsorted(rc, 60 * L); k = np.searchsorted(zh, 285 * L)
    j = int(round(rc[i] / dr)); kk = int(round(zh[k] / dz))
    e3 = abs(A[j, kk] - ref[i, k]) / n0
    print(f"  uniform {rec} at ({rc[i] / L:.1f}, {zh[k] / L:.1f}) k0^-1 vs data file: {e3:.1e} of n0")
    ok &= e3 < 0.02
sys.exit(0 if ok else 1)
PY
cd - > /dev/null
rm -rf "$WORK"
echo "openPMD test passed"
