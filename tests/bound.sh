#!/bin/sh
# Bound electrons (block 'bound_electrons', Maxwell-Bloch medium of a V-type atom) against independent references
# (tests/bound.py):
#  1. Rabi flopping (tests/bound_rabi.deck): a thin resonant medium, Gaussian pulses of area pi/2, pi, 3pi/2, 2pi:
#     the excited population after the pulse against the two-level Schroedinger equation (no rotating-wave
#     approximation) integrated in the code's own field (pi, 2pi: 1e-5; pi/2, 3pi/2: 3e-3 at lambda/32);
#  2. weak field in Rb vapour (tests/bound_linear.deck, D1 + D2, 840 nm): the spectral phase after 50 um against the
#     Lorentz susceptibility of the two lines (diffraction removed with a vacuum run): 0.4 %;
#  3. static limit (tests/bound_static.deck): a line at 10 omega0 against the dielectric of the same chi: the group
#     delays differ by the dispersion of the line (9 %);
#  4. ionization (tests/bound_ionize.deck, Rb vapour as a neutral ion species): the neutral fraction on the axis
#     against the V-type atom with ADK losses per level (5S 4.177 eV, 5P 2.6 eV) integrated in the code's field
#     (dzeta = lambda/64: 4e-3 at a0 = 4e-4, 42 % ionized almost only through the excited state; 6e-3 at a0 = 1e-3,
#     97 % ionized), and the mean
#     charge state of the macro-atoms leaving the window against the grid's area mean (1e-3).
# usage: sh tests/bound.sh <gempic_rz> [backend]     (exit 77: skipped, needs python3 with numpy)
set -e
abspath() { (cd "$(dirname "$1")" && echo "$(pwd)/$(basename "$1")"); }
EXE=$(abspath "$1")
BACKEND=${2:-serial}
SRC=$(cd "$(dirname "$0")/.." && pwd)
python3 -c "import numpy" 2>/dev/null || { echo "bound test skipped: needs python3 with numpy"; exit 77; }
WORK=$(mktemp -d)
cd "$WORK"
export OMP_PROC_BIND=false
PY="python3 $SRC/tests/bound.py"

echo "--- 1. Rabi flopping, pulse areas pi/2, pi, 3pi/2, 2pi"
for A in 0.5:4.394083e-04:5e-3 1:8.788166e-04:1e-3 1.5:1.318225e-03:5e-3 2:1.757633e-03:1e-3; do
  a=${A%%:*}; r=${A#*:}; a0=${r%%:*}; tol=${r#*:}
  "$EXE" "$SRC/tests/bound_rabi.deck" numerics.backend=$BACKEND laser.a0=$a0 output.dir=rabi$a > rabi$a.log
  $PY rabi rabi$a/map_00001.txt rabi$a.log 0.8 0.8 1e-3 $tol
done

echo "--- 2. weak field in Rb vapour (D1, D2)"
"$EXE" "$SRC/tests/bound_linear.deck" numerics.backend=$BACKEND laser.a0=1e-6 > lin.log
sed '/^bound_electrons/d' "$SRC/tests/bound_linear.deck" > vlin.deck
"$EXE" vlin.deck numerics.backend=$BACKEND laser.a0=1e-6 output.dir=vlin > vlin.log
$PY linear linear/axis_00000.txt linear/axis_00100.txt vlin/axis_00000.txt vlin/axis_00100.txt 50 0.84 4e17 0.0125 \
  0.01 0.794979 0.342 0.780241 0.695

echo "--- 3. static limit against the dielectric"
"$EXE" "$SRC/tests/bound_static.deck" numerics.backend=$BACKEND output.dir=sb > sb.log
sed 's/^bound_electrons.*/dielectric { epsilon = 1.01; r_out = 60 um }/' "$SRC/tests/bound_static.deck" > sd.deck
"$EXE" sd.deck numerics.backend=$BACKEND output.dir=sd > sd.log
sed '/^bound_electrons/d' "$SRC/tests/bound_static.deck" > sv.deck
"$EXE" sv.deck numerics.backend=$BACKEND output.dir=sv > sv.log
$PY static sb/axis_00000.txt sb/axis_00060.txt sd/axis_00000.txt sd/axis_00060.txt sv/axis_00000.txt sv/axis_00060.txt \
  30 0.8 0.08 1 1.7246e21 0.15
echo "--- 4. ionization: depopulation of the bound electrons, ADK per level (Rb, 33 fs, a0 = 4e-4 and 1e-3)"
for a0 in 4e-4 1e-3; do
  "$EXE" "$SRC/tests/bound_ionize.deck" numerics.backend=$BACKEND laser.a0=$a0 grid.dzeta=0.0125um output.dir=ion$a0     > ion$a0.log
  $PY ionize ion$a0/map_00001.txt ion$a0.log 0.8 4.177128 1e-2 5e-3 0.794979 0.342 0.780241 0.695
done
echo "bound test passed"
