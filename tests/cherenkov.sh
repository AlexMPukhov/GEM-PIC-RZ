#!/bin/sh
# Numerical Cherenkov test.  The GT2 scheme has a zero-frequency mode at the
# zeta-Nyquist wavenumber (kh = pi) for every k_perp, which is resonantly
# driven by the shot noise of a zeta-stationary relativistic bunch; the beam
# current filter along zeta (numerics.beam_filter) removes this source.
# Measured with `h5_extract energy`: W_nyquist / W_total at the end of
#   1) a 1 pC, gamma = 1000 bunch in vacuum, 300 steps, filter on   (< 1e-3)
#   2) the same with the filter off                                  (> 1e-2,
#      checks that the diagnostic sees the effect)
#   3) a bunch-driven plasma wake reaching the rear end of the window,
#      100 steps, filter on                                          (< 1e-3)
# usage: sh tests/cherenkov.sh <gempic_rz> <h5_extract> [key=value ...]
set -e
abspath() { (cd "$(dirname "$1")" && echo "$(pwd)/$(basename "$1")"); }
EXE=$(abspath "$1"); H5X=$(abspath "$2"); shift 2
SRC=$(cd "$(dirname "$0")/.." && pwd)
WORK=$(mktemp -d)
cd "$WORK"
ratio() {  # <file> -> W_nyquist / W_total
  "$H5X" "$1" energy | awk '!/^#/ { printf "%.3e", $4 / $3 }'
}
check() {  # <label> <ratio> <lt|gt> <bound>
  if awk -v r="$2" -v b="$4" -v op="$3" 'BEGIN { exit !((op == "lt") ? (r < b) : (r > b)) }'; then
    echo "  $1: W_nyquist/W_total = $2 ($3 $4)  ok"
  else
    echo "  $1: W_nyquist/W_total = $2, expected $3 $4  FAILED"; exit 1
  fi
}
"$EXE" "$SRC/tests/ncr_vacuum.deck" "$@" > vac.log
check "vacuum, filter 1" "$(ratio vac/data_00300.h5)" lt 1e-3
"$EXE" "$SRC/tests/ncr_vacuum.deck" "$@" numerics.beam_filter=0 output.dir=vac0 > vac0.log
check "vacuum, filter 0" "$(ratio vac0/data_00300.h5)" gt 1e-2
"$EXE" "$SRC/tests/ncr_plasma.deck" "$@" > pla.log
check "plasma, filter 1" "$(ratio pla/data_00100.h5)" lt 1e-3
cd - > /dev/null
rm -rf "$WORK"
echo "cherenkov test passed"
