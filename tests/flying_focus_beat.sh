#!/bin/sh
# Flying focus of two luminal Bessel components (benchmarks/flying_focus/beat_vacuum.deck: J0 at 40 mrad + 0.2 J1 e^{i phi}
# at 100 mrad, beat period 190.7 um), fed through the absorber layer, against the analytic field (check_physics ffbeat:
# E_x in the plane y = 0 at z = 40, 80, 120 um, independent frequency quadrature):
#   1. circular, helicity +1 (J1 part in mode 0), m_max = 1
#   2. circular, helicity -1, Bessel order -1 (again mode 0), m_max = 1
#   3. linear (J1 part in modes 0 and 2), m_max = 2
# The field fed through the layer at the level s is, after the step, the incident field of s - dz/2 (a half-step lag
# of the feeding, also present for a single component): the reference is taken there (ds = -2.5 um). Tolerances:
# rms(E_x - reference) < 3e-3 a0 within 2 L of the pulse centre and r < 40 um (measured <= 1.5e-3); the odd part in x
# (the J1 component, which carries the beat) to 12 % of its rms (measured <= 8 %, numerical dispersion of the 100 mrad
# cone, ~1.5 % in kappa); the intensity maximum near the axis to 0.1 um (it moves on a circle of 1.15 um).
# usage: sh tests/flying_focus_beat.sh <gempic_rz> <check_physics> [key=value ...]
set -e
abspath() { (cd "$(dirname "$1")" && echo "$(pwd)/$(basename "$1")"); }
EXE=$(abspath "$1"); CHK=$(abspath "$2"); shift 2
SRC=$(cd "$(dirname "$0")/.." && pwd)
WORK=$(mktemp -d)
cd "$WORK"
run() {  # <label> <beat order> <key=value ...>
  label=$1; l=$2; shift 2
  echo "$label"
  "$EXE" "$SRC/benchmarks/flying_focus/beat_vacuum.deck" output.dir=out "$@" > run.log
  for f in out/map_00008.txt out/map_00016.txt out/map_00024.txt; do
    "$CHK" ffbeat $f 3e-3 0.12 0.1 -2.5 1 7.639 32 0.04 0 1 0 0.10 $l 0.2 0
  done
  rm -rf out
}
run "circular, helicity +1, J0 + J1 e^{i phi}, m_max = 1" 1 "$@"
run "circular, helicity -1, J0 + J1 e^{-i phi}, m_max = 1" -1 flying_focus.helicity=-1 flying_focus.beat_order=-1 "$@"
run "linear, J0 + J1 e^{i phi}, m_max = 2" 1 flying_focus.polarization=linear grid.m_max=2 "$@"
cd - > /dev/null
rm -rf "$WORK"
echo "flying focus beat test passed"
