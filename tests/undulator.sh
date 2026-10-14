#!/bin/sh
# Undulator radiation (tests/undulator.deck): coherent spontaneous emission of a
# short bunch in a planar undulator, against the far-field radiation integral
# (Jackson 14.65) of the same trajectory with the form factor of the bunch:
#   1. radiated energy (modes m >= 1) 1.103e-8 code units (+- 5 %)
#   2. peak of the angle-integrated spectrum (h5_extract spectrum) at k = 0.981 k0 (+- 0.3 %):
#      the off-axis emission is red-shifted, on axis lambda = lambda_r = lambda0
#   3. slippage of the bunch: 8.260 um (+- 0.03 um)
#   4. the FEL table fel.txt: radiation energy as in the data file, bunching factor 0.7346
# usage: sh tests/undulator.sh <gempic_rz> <h5_extract> [key=value ...]
set -e
abspath() { (cd "$(dirname "$1")" && echo "$(pwd)/$(basename "$1")"); }
EXE=$(abspath "$1"); H5X=$(abspath "$2"); shift 2
SRC=$(cd "$(dirname "$0")/.." && pwd)
WORK=$(mktemp -d)
cd "$WORK"
"$EXE" "$SRC/tests/undulator.deck" "$@" > run.log
F=undulator_out/data_00210.h5
check() {  # <label> <value> <expected> <tolerance> <unit>
  if awk -v v="$2" -v e="$3" -v t="$4" 'BEGIN { d = v - e; if (d < 0) d = -d; exit !(d <= t) }'; then
    echo "  $1 = $2 $5 (expected $3 +- $4)  ok"
  else
    echo "  $1 = $2 $5, expected $3 +- $4  FAILED"; exit 1
  fi
}
W=$("$H5X" $F energy 0 30 | awk '!/^#/ { printf "%.4e", $6 }')
check "radiated energy" "$W" 1.103e-8 0.055e-8 "(code units)"
K=$("$H5X" $F spectrum 0 30 0.90 1.06 161 | awk '!/^#/ && $3 > best { best = $3; k = $1 } END { printf "%.4f", k }')
check "peak of the spectrum k/k0" "$K" 0.981 0.003 ""
S=$(awk '!/^#/ && $3 == "b" { if (!seen) z0 = $12; seen = 1; z1 = $12 } END { printf "%.3f", z0 - z1 }' \
  undulator_out/bunches.txt)
check "slippage" "$S" 8.260 0.03 um
# FEL table (fel.txt): radiation energy = the energy of the modes m >= 1 in the file (J: x eps0 E0^2 / k0^3),
# bunching at lambda_r of the short bunch = its form factor exp(-(k sigma_z)^2 / 2) = 0.7346
WF=$(awk '!/^#/ { w = $3 } END { print w }' undulator_out/fel.txt)
WH=$("$H5X" $F energy | awk '!/^#/ { printf "%.6e", $6 * 2.943736e-7 }')
check "fel.txt radiation energy / data file" "$(awk -v a="$WF" -v b="$WH" 'BEGIN { printf "%.5f", a / b }')" 1 0.0005 ""
check "fel.txt bunching" "$(awk '!/^#/ { b = $7 } END { printf "%.4f", b }' undulator_out/fel.txt)" 0.7346 0.002 ""
cd - > /dev/null
rm -rf "$WORK"
echo "undulator test passed"
