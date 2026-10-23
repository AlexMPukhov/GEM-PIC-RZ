#!/bin/sh
# bunch.axis_refine: more macroparticles near the axis, same charge density.
# A Gaussian bunch (examples/undulator_emission.deck, sigma_r = 100 um, radial cells of 20 um) with axis_refine = 9
# inside 60 um: the number of macroparticles grows by 8 x the fraction inside, 1 - exp(-r_a^2 / 2 sigma^2) = 0.1647,
# the charge stays exactly, and the rms size agrees with the unrefined bunch (+-2 %).
# usage: sh tests/axis_refine.sh <gempic_rz> [backend]
set -e
abspath() { (cd "$(dirname "$1")" && echo "$(pwd)/$(basename "$1")"); }
EXE=$(abspath "$1"); BACKEND=${2:-openmp}
SRC=$(cd "$(dirname "$0")/.." && pwd)
D="$SRC/examples/undulator_emission.deck"
WORK=$(mktemp -d)
cd "$WORK"
O="numerics.backend=$BACKEND grid.z_end=1.5mm output.format=text output.every=1 bunch.macroparticles=200000"
# shellcheck disable=SC2086
"$EXE" "$D" $O output.dir=a > a.log
# shellcheck disable=SC2086
"$EXE" "$D" $O bunch.axis_refine=9 bunch.axis_radius=60um output.dir=b > b.log
grep -q "axis refinement x9" b.log || { echo "  no axis refinement reported  FAILED"; exit 1; }
row() { awk '$1 == 0 && $3 == "b"' "$1/bunches.txt"; }
NA=$(row a | awk '{print $4}'); NB=$(row b | awk '{print $4}')
QA=$(row a | awk '{print $5}'); QB=$(row b | awk '{print $5}')
SA=$(row a | awk '{print $8}'); SB=$(row b | awk '{print $8}')
echo "  macroparticles $NA -> $NB, charge $QA -> $QB pC, sigma_x $SA -> $SB um"
awk -v a="$NA" -v b="$NB" 'BEGIN { r = (b / a - 1) / 8; exit !(r > 0.155 && r < 0.175) }' ||
  { echo "  added fraction $(awk -v a="$NA" -v b="$NB" 'BEGIN{print (b/a-1)/8}') vs 0.1647  FAILED"; exit 1; }
awk -v a="$QA" -v b="$QB" 'BEGIN { d = b / a - 1; if (d < 0) d = -d; exit !(d < 1e-9) }' || { echo "  charge changed  FAILED"; exit 1; }
awk -v a="$SA" -v b="$SB" 'BEGIN { d = b / a - 1; if (d < 0) d = -d; exit !(d < 0.02) }' || { echo "  sigma_x changed  FAILED"; exit 1; }
echo "  ok"
cd /
rm -rf "$WORK"
