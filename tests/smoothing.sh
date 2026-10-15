#!/bin/sh
# Radial smoothing of the plasma sources (plasma.smooth_length), on tests/regression.deck:
#   1. a = 1e-9 k0^-1 reproduces the cold-plasma reference (tests/reference)
#   2. a = 1 k0^-1 caps the axis density spike at the back of the bubble (124 n0 cold;
#      below 40 n0 smoothed) and changes E_z in the front half of the window by < 5e-3 of its maximum (2e-3)
#   3. the chosen backend agrees with the serial run (a = 1 k0^-1)
# The filter itself (commutation with the divergence, conservation, constants) is checked in
# test_deck.
# usage: sh tests/smoothing.sh <gempic_rz> <compare_columns> [backend]
set -e
abspath() { (cd "$(dirname "$1")" && echo "$(pwd)/$(basename "$1")"); }
EXE=$(abspath "$1"); CMP=$(abspath "$2"); BACKEND=${3:-serial}
SRC=$(cd "$(dirname "$0")/.." && pwd)
WORK=$(mktemp -d)
cd "$WORK"
R="$SRC/tests/regression.deck"
# 1
"$EXE" "$R" numerics.backend=$BACKEND plasma.smooth_length=1e-9k0^-1 output.dir=tiny > tiny.log
"$CMP" "$SRC/tests/reference/axis_00005.txt" tiny/axis_00005.txt 1e-6 1e-9
# 2
"$EXE" "$R" numerics.backend=serial plasma.smooth_length=1k0^-1 output.dir=ser > ser.log
grep -q "smoothing" ser/summary.txt
awk '!/^#/ { if ($4 > m) m = $4 } END { printf "axis density peak: %.1f n0 (smoothed)\n", m; exit !(m < 40) }' ser/axis_00005.txt
awk '!/^#/ { if ($4 > m) m = $4 } END { printf "axis density peak: %.1f n0 (cold)\n", m; exit !(m > 100) }' \
  "$SRC/tests/reference/axis_00005.txt"
# fields ahead of the bubble back (front half of the window): smoothing changes them little
paste "$SRC/tests/reference/axis_00005.txt" ser/axis_00005.txt | awk '!/^#/ && NR > 600 {
    d = $3 - $(NF - 2); if (d < 0) d = -d; if (d > dm) dm = d
    a = $3; if (a < 0) a = -a; if (a > am) am = a }
  END { printf "E_z ahead of the bubble back: max change %.2e of %.2e\n", dm, am; exit !(dm < 5e-3 * am) }'
# 3
if [ "$BACKEND" != serial ]; then
  "$EXE" "$R" numerics.backend=$BACKEND plasma.smooth_length=1k0^-1 output.dir=par > par.log
  "$CMP" ser/axis_00005.txt par/axis_00005.txt 1e-6 1e-9
  "$CMP" ser/bunches.txt par/bunches.txt 1e-6 1e-9
fi
cd - > /dev/null
rm -rf "$WORK"
echo "smoothing test passed"
