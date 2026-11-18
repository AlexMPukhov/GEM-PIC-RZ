#!/bin/sh
# openPMD files read with openpmd-viewer (the reader of FBPIC users): examples/ionization_injection.deck (He / N
# mixture, a0 = 2, the K-shell electrons of the N5+ dopant are freed near the peak and trapped in the bubble) for
# 10 steps (z = 50.9 um) with output every 2 steps, format = hdf5, openpmd, default openpmd_grid (uniform only).
# The species inj_N5+ (trapped electrons freed by ionization) is empty in the first files and filled from step 8 on;
# trapped (self-injected background electrons) stays empty. tests/openpmd_viewer_check.py checks with
# openpmd_viewer.OpenPMDTimeSeries:
#   - every iteration opens, the species list is the same in all of them, no native_* records by default;
#   - E_x, E_y, E_z at (r, phi) from get_field (m = all, theta) = the code's own reconstruction from the native
#     data file (h5_extract lineout of data_NNNNN.h5) to 1e-3 of the peak;
#   - per species and iteration: number of macroparticles, charge, mean lab z, mean u_z, from get_particle,
#     = the beam list of data_NNNNN.h5 (empty species: zero-length arrays);
#   - z = const metadata: snapshotType = position, particle time = (s - zeta) / c.
# Exit code 77 (skipped) without python3 / numpy / h5py / openpmd_viewer.
# usage: sh tests/openpmd_viewer.sh <gempic_rz> <h5_extract>
set -e
abspath() { (cd "$(dirname "$1")" && echo "$(pwd)/$(basename "$1")"); }
EXE=$(abspath "$1"); H5X=$(abspath "$2")
SRC=$(cd "$(dirname "$0")/.." && pwd)
python3 -c "import numpy, h5py" 2>/dev/null || { echo "python3 with numpy and h5py not found: skipped"; exit 77; }
python3 -c "import openpmd_viewer" 2>/dev/null || { echo "openpmd_viewer not found: skipped"; exit 77; }
WORK=$(mktemp -d)
cd "$WORK"
"$EXE" "$SRC/examples/ionization_injection.deck" grid.z_end=50.9um output.every=2 output.maps=false \
  output.format=hdf5,openpmd output.dir=out > run.log
python3 "$SRC/tests/openpmd_viewer_check.py" out "$H5X"
cd - > /dev/null
rm -rf "$WORK"
echo "openpmd_viewer test passed"
