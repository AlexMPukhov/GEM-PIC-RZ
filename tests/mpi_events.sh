#!/bin/sh
# Events of the MPI pipeline that act on the whole beam (synchronization points, pipeline.cpp: drain, syncEvent):
#   bunch.z_start: the witness of tests/regression.deck loaded at the end of step 2 (z = 80 k0^-1) in the plasma
#   wake. Reference: the serial run (Simulation::loadBunches). Checked: the pipelined step on one rank = serial
#   (to all digits), 2-4 ranks = serial (rounding: 1e-6), also with beam_field = mean, with a checkpoint before
#   the loading (3 ranks, restart on 2) and one after it (restart: the bunch is not loaded again), and on the
#   device (GPU emulation build) on 3 ranks.
# usage: sh tests/mpi_events.sh <gempic_rz> <compare_columns> <mpiexec> [extra mpiexec flags]
set -e
abspath() { (cd "$(dirname "$1")" && echo "$(pwd)/$(basename "$1")"); }
EXE=$(abspath "$1"); CMP=$(abspath "$2"); MPIEXEC=${3:-mpirun}; shift 3 2>/dev/null || true
FLAGS="$*"
if [ "$(id -u)" = 0 ]; then export OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1; fi
if $MPIEXEC --version 2>&1 | grep -qiE "open mpi|openrte|open-mpi"; then FLAGS="$FLAGS --oversubscribe"; fi
run() { np=$1; shift; $MPIEXEC $FLAGS -n $np "$EXE" "$@"; }
runx() { np=$1; x=$2; shift 2; $MPIEXEC $FLAGS -n $np "$x" "$@"; }
SRC=$(cd "$(dirname "$0")/.." && pwd)
WORK=$(mktemp -d)
cd "$WORK"
DECK="$SRC/tests/regression.deck"
compare() {  # <reference dir> <dir> <tolerances>
  for f in axis_00002.txt axis_00003.txt axis_00004.txt map_00002.txt map_00004.txt bunches.txt; do
    "$CMP" $1/$f $2/$f $3 $4
  done
}
# ---- bunch.z_start
Z="numerics.backend=serial output.every=1 output.map_stride=8 grid.z_end=160k0^-1 bunch.z_start=80k0^-1"
"$EXE" "$DECK" $Z output.dir=ser > ser.log
grep -q "bunch 'witness' loaded at z" ser.log
GEMPIC_PIPELINE=1 "$EXE" "$DECK" $Z output.dir=p1 > p1.log
grep -q "bunch 'witness' loaded at z" p1.log
compare ser p1 0 0
for np in 2 3 4; do
  run $np "$DECK" $Z output.dir=mp$np > mp$np.log
  grep -q "bunch 'witness' loaded at z" mp$np.log
  compare ser mp$np 1e-6 1e-9
done
echo "  z_start: 1 rank pipelined = serial (exact), 2-4 ranks = serial"
ZM="$Z numerics.beam_field=mean numerics.beam_filter=2"
"$EXE" "$DECK" $ZM output.dir=mser > mser.log
GEMPIC_PIPELINE=1 "$EXE" "$DECK" $ZM output.dir=mp1 > mp1.log
compare mser mp1 0 0
for np in 2 3 4; do
  run $np "$DECK" $ZM output.dir=mmp$np > mmp$np.log
  compare mser mmp$np 1e-6 1e-9
done
echo "  z_start, beam_field = mean: 1 rank pipelined = serial (exact), 2-4 ranks = serial"
# checkpoint before the loading (step 1) and after it (step 3)
run 3 "$DECK" $Z output.dir=ck1 grid.z_end=40k0^-1 output.checkpoint_every=1 > ck1a.log
run 2 "$DECK" $Z output.dir=ck1 --restart latest > ck1b.log
grep -q "bunch 'witness' loaded at z" ck1b.log
compare ser ck1 1e-6 1e-9
run 3 "$DECK" $Z output.dir=ck3 grid.z_end=120k0^-1 output.checkpoint_every=3 > ck3a.log
run 2 "$DECK" $Z output.dir=ck3 --restart latest > ck3b.log
if grep -q "loaded at z" ck3b.log; then echo "FAIL: the bunch was loaded again after the restart"; exit 1; fi
for f in axis_00004.txt map_00004.txt; do "$CMP" ser/$f ck3/$f 1e-6 1e-9; done
echo "  z_start across checkpoints (3 -> 2 ranks): before and after the loading = serial"
if [ -x "$(dirname "$EXE")/gempic_rz_gpuemu" ]; then
  E="$(echo $Z | sed s/backend=serial/backend=openmp/)"
  GEMPIC_PIPELINE_ENGINE=1 OMP_NUM_THREADS=2 runx 3 "$(dirname "$EXE")/gempic_rz_gpuemu" "$DECK" $E output.dir=dv > dv.log
  compare ser dv 1e-6 1e-9
  GEMPIC_PIPELINE_ENGINE=1 OMP_NUM_THREADS=2 runx 3 "$(dirname "$EXE")/gempic_rz_gpuemu" "$DECK" $E \
    numerics.beam_field=mean numerics.beam_filter=2 output.dir=dvm > dvm.log
  compare mser dvm 1e-6 1e-9
  echo "  z_start on the device (GPU emulation, 3 ranks, also beam_field = mean) = serial"
fi
cd - > /dev/null
rm -rf "$WORK"
echo "mpi events test passed"
