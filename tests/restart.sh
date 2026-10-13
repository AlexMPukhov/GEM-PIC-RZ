#!/bin/sh
# Checkpoint / restart: a run interrupted at a checkpoint and restarted must
# reproduce the uninterrupted run (tests/regression.deck, 4 steps):
#   1. checkpoint at step 2, restart from it (--restart <file>)
#   2. stop at the wall-time limit (exit code 3), continue with --restart latest
# usage: sh tests/restart.sh <gempic_rz> <compare_columns> [backend]
#   serial, openmp: results must be identical; cuda: equal to rounding
set -e
abspath() { (cd "$(dirname "$1")" && echo "$(pwd)/$(basename "$1")"); }
EXE=$(abspath "$1"); CMP=$(abspath "$2"); BACKEND=${3:-serial}
if [ "$BACKEND" = cuda ]; then TOL="1e-6 1e-9"; else TOL="0 0"; fi
SRC=$(cd "$(dirname "$0")/.." && pwd)
WORK=$(mktemp -d)
cd "$WORK"
DECK="$SRC/tests/regression.deck"
RUN="numerics.backend=$BACKEND output.every=1 output.map_stride=8 grid.z_end=160k0^-1"
compare() {
  for f in axis_00003.txt axis_00004.txt map_00004.txt beam_00004.txt bunches.txt; do "$CMP" full/$f $1/$f $TOL; done
}
"$EXE" "$DECK" $RUN output.dir=full > full.log
# 1. checkpoint at step 2 (the first run ends there), restart from the file
"$EXE" "$DECK" $RUN output.dir=a grid.z_end=80k0^-1 output.checkpoint_every=2 > a1.log
test -f a/checkpoint_00002.h5
"$EXE" "$DECK" $RUN output.dir=a --restart a/checkpoint_00002.h5 > a2.log
compare a
# 2. wall-time limit: stops after the first step with a checkpoint (exit code 3)
rc=0
"$EXE" "$DECK" $RUN output.dir=w output.walltime_hours=1e-5 > w1.log || rc=$?
if [ $rc -ne 3 ]; then echo "expected exit code 3 at the wall-time limit, got $rc"; exit 1; fi
grep -q "wall-time limit" w1.log
"$EXE" "$DECK" $RUN output.dir=w --restart latest > w2.log
compare w
cd - > /dev/null
rm -rf "$WORK"
echo "restart test passed ($BACKEND)"
