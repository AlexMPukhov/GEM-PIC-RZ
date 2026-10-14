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
  output.format=hdf5,openpmd output.openpmd_dr=0.1um output.openpmd_dzeta=0.25k0^-1 output.dir=out > run.log
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
cd - > /dev/null
rm -rf "$WORK"
echo "openPMD test passed"
