#!/bin/sh
# Free-hand grids along zeta and z (grid.zeta_grid = table, grid.z_grid = table):
#   1. tables with a constant step reproduce the uniform run (tests/regression.deck)
#   2. tests/flexgrid.deck (non-uniform zeta and z): the given backend agrees
#      with the serial code (to rounding)
#   3. checkpoint / restart on the non-uniform grids reproduces the
#      uninterrupted run; a restart may continue with a different z grid
# usage: sh tests/flexgrid.sh <gempic_rz> <compare_columns> [backend]
set -e
abspath() { (cd "$(dirname "$1")" && echo "$(pwd)/$(basename "$1")"); }
EXE=$(abspath "$1"); CMP=$(abspath "$2"); BACKEND=${3:-openmp}
SRC=$(cd "$(dirname "$0")/.." && pwd)
WORK=$(mktemp -d)
cd "$WORK"
# 1. constant-step tables = uniform grid
sed -e 's/^  dzeta .*/  zeta_grid = table\n  zeta_points = 0, 100 k0^-1\n  dzeta_points = 0.25, 0.25 k0^-1/' \
    -e 's/^  dz  .*/  z_grid = table\n  z_points = 0 k0^-1\n  dz_points = 40 k0^-1/' "$SRC/tests/regression.deck" > table.deck
U="numerics.backend=serial grid.z_end=80k0^-1 output.every=1 output.map_stride=8"
"$EXE" "$SRC/tests/regression.deck" $U output.dir=uni > uni.log
"$EXE" table.deck $U output.dir=tab > tab.log
for f in axis_00002.txt map_00002.txt bunches.txt; do "$CMP" uni/$f tab/$f 1e-9 1e-12; done
# 2. non-uniform grids: backend against serial
F="output.map_stride=8"
"$EXE" "$SRC/tests/flexgrid.deck" $F numerics.backend=serial output.dir=ser > ser.log
"$EXE" "$SRC/tests/flexgrid.deck" $F numerics.backend=$BACKEND output.dir=par > par.log
# (beam files are not compared: the GPU code paths append particles in any order)
for f in axis_00004.txt map_00004.txt bunches.txt; do "$CMP" ser/$f par/$f 1e-6 1e-9; done
# 3. restart (serial, exact): checkpoint at step 2, continue to the end
"$EXE" "$SRC/tests/flexgrid.deck" $F numerics.backend=serial output.dir=ck grid.z_end=40k0^-1 \
  output.checkpoint_every=2 > ck1.log
test -f ck/checkpoint_00002.h5
"$EXE" "$SRC/tests/flexgrid.deck" $F numerics.backend=serial output.dir=ck --restart latest > ck2.log
for f in axis_00004.txt map_00004.txt beam_00004.txt bunches.txt; do "$CMP" ser/$f ck/$f 0 0; done
#    ... with a different z grid after the restart: uniform dz = 10 k0^-1 from z = 42.1 to 120 k0^-1
sed -e 's/^  z_grid .*//' -e 's/^  z_points          = 0, 200 k0^-1//' -e 's/^  dz_points .*/  dz = 10 k0^-1/' \
  "$SRC/tests/flexgrid.deck" > dz10.deck
"$EXE" dz10.deck $F numerics.backend=serial output.dir=ck --restart ck/checkpoint_00002.h5 > ck3.log
grep -q "done: 10 steps" ck3.log || { echo "restart with dz = 10 k0^-1: expected 10 steps"; tail -2 ck3.log; exit 1; }
cd - > /dev/null
rm -rf "$WORK"
echo "flexgrid test passed ($BACKEND)"
