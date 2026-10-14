#!/bin/sh
# openPMD output (output.format = openpmd): tests/regression.deck for 2 steps
#   1. the files pass the openPMD validator (openPMD_check_h5, if installed)
#   2. tests/openpmd_check.py (python3 with numpy, h5py): the standard's attributes,
#      the fields rebuilt from the thetaMode data at (r, theta) (output nodes on the slices) equal the code's own
#      reconstruction (h5_extract lineout of data_NNNNN.h5) to 1e-3 of the peak,
#      lab z = zeta + s, beam charge and particle count
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
python3 "$SRC/tests/openpmd_check.py" out 2 "$H5X"
cd - > /dev/null
rm -rf "$WORK"
echo "openPMD test passed"
