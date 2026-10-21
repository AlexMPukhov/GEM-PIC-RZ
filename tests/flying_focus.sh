#!/bin/sh
# Flying focus (benchmarks/flying_focus/vacuum.deck, 0.5 mm, feeding layer strength 1): the luminal Bessel pulse
# fed through the radial absorber layer must keep its on-axis pattern:
#   1. peak |E_x| on the axis = a0 = 2 (+- 2 %)
#   2. the envelope (centroid of E_x^2 within +- 13 um of zeta_c = 90 um) stays at zeta_c (+- 0.05 um): moves at c
# usage: sh tests/flying_focus.sh <gempic_rz> [key=value ...]
set -e
abspath() { (cd "$(dirname "$1")" && echo "$(pwd)/$(basename "$1")"); }
EXE=$(abspath "$1"); shift
SRC=$(cd "$(dirname "$0")/.." && pwd)
WORK=$(mktemp -d)
cd "$WORK"
"$EXE" "$SRC/benchmarks/flying_focus/vacuum.deck" grid.z_end=0.5mm grid.absorber_strength=1 output.dir=out \
  output.every=0 output.maps=false "$@" > run.log
F=$(ls out/axis_*.txt | tail -1)
check() {  # <label> <value> <expected> <tolerance>
  if awk -v v="$2" -v e="$3" -v t="$4" 'BEGIN { d = v - e; if (d < 0) d = -d; exit !(d <= t) }'; then
    echo "  $1 = $2 (expected $3 +- $4)  ok"
  else
    echo "  $1 = $2, expected $3 +- $4  FAILED"; exit 1
  fi
}
check "peak |E_x| on the axis" "$(awk '!/^#/ { v = $2 < 0 ? -$2 : $2; if (v > m) m = v } END { printf "%.4f", m }' $F)" 2 0.04
check "envelope centroid" "$(awk '!/^#/ && $1 > 77 && $1 < 103 { w += $2 * $2; s += $1 * $2 * $2 } END { printf "%.4f", s / w }' $F)" 90 0.05
cd - > /dev/null
rm -rf "$WORK"
echo "flying focus test passed"
