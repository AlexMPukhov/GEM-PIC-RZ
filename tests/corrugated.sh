#!/bin/sh
# Corrugated metal pipe as an equivalent thin dielectric layer (tests/corrugated.py): bore 1 mm, depth 60 um,
# period 230 um, gap 115 um -> eps = p/(p - g) = 2 in a 60 um layer on the wall. The on-axis wake of a 1 nC, 50 um
# bunch against the modal theory of the layer (wavenumber 1e-3, amplitude 1e-3); the leading-order corrugated-pipe
# formulas of Bane and Stupakov (k0, kappa) are printed (2.5 % / 7 % away at h/a = 0.06, converging as h/a -> 0).
# usage: sh tests/corrugated.sh <gempic_rz> [backend]     (exit 77: skipped, needs python3 with numpy)
set -e
abspath() { (cd "$(dirname "$1")" && echo "$(pwd)/$(basename "$1")"); }
EXE=$(abspath "$1")
BACKEND=${2:-serial}
SRC=$(cd "$(dirname "$0")/.." && pwd)
python3 -c "import numpy" 2>/dev/null || { echo "corrugated test skipped: needs python3 with numpy"; exit 77; }
WORK=$(mktemp -d)
cd "$WORK"
export OMP_PROC_BIND=false
python3 "$SRC/tests/corrugated.py" theory
"$EXE" "$SRC/tests/corrugated.deck" numerics.backend=$BACKEND > corr.log
python3 "$SRC/tests/corrugated.py" check corr/axis_00080.txt 1000 60 230 115 1000 50 4800 39000 10 0.005 0.01
echo "corrugated test passed"
