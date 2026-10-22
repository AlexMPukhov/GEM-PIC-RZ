#!/bin/sh
# Beam push in a planar undulator with an intense beam (benchmarks/harmonics/fel.deck: 4 kA, gamma = 1000,
# sigma_r = 10 um, zero emittance, K = 1.5), field step dz = lambda_u / 16, 5 cm (5 periods).
# The beam's self-field turns with the wiggle. Pushed in the fields of s^0 (numerics.beam_field = level),
# the beam blows up transversely: eps_nx reaches 0.14 um at 5 cm (2.8 um at 25 cm). Pushed again in the
# mean of the fields of s^0 and s^1 (beam_field = mean, the default with undulators), eps_nx stays at
# 0.008 um, close to the converged value (level at lambda_u / 64: 0.005 um).
# Checks: mean < 0.02 um, level > 0.07 um (the test still sees the effect).
# usage: sh tests/planar_step.sh <gempic_rz> [backend]
set -e
abspath() { (cd "$(dirname "$1")" && echo "$(pwd)/$(basename "$1")"); }
EXE=$(abspath "$1"); BACKEND=${2:-openmp}
SRC=$(cd "$(dirname "$0")/.." && pwd)
D="$SRC/benchmarks/harmonics/fel.deck"
WORK=$(mktemp -d)
cd "$WORK"
O="numerics.backend=$BACKEND grid.z_end=5cm grid.dzeta=1.328nm grid.dz=0.625mm numerics.beam_dz_max=0.625mm
   bunch.macroparticles=128000 output.format=text output.every=100000"
eps() { awk '$3 == "b" { e = $10 } END { print e }' "$1/bunches.txt"; }
# shellcheck disable=SC2086
"$EXE" "$D" $O output.dir=mean > mean.log
grep -q "beam field mean" mean.log || { echo "  beam_field is not 'mean' by default with an undulator  FAILED"; exit 1; }
# shellcheck disable=SC2086
"$EXE" "$D" $O numerics.beam_field=level output.dir=level > level.log
EM=$(eps mean); EL=$(eps level)
echo "  eps_nx at 5 cm: mean $EM um, level $EL um"
awk -v m="$EM" 'BEGIN { exit !(m < 0.02) }' || { echo "  mean: eps_nx = $EM um >= 0.02 um  FAILED"; exit 1; }
awk -v l="$EL" 'BEGIN { exit !(l > 0.07) }' || { echo "  level: eps_nx = $EL um <= 0.07 um (effect not seen)  FAILED"; exit 1; }
echo "  ok"
cd /
rm -rf "$WORK"
