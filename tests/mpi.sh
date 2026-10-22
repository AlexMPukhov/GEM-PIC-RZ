#!/bin/sh
# MPI pipeline along zeta (tests/regression.deck, 4 steps):
#   1. the pipelined step on one rank (GEMPIC_PIPELINE=1) is identical to the serial step
#   2. 2, 3 and 4 ranks agree with the serial run (to rounding: 1e-6 relative)
#   3. checkpoint on 3 ranks, restart on 2 ranks / serially; wall-time stop on 3 ranks
#   4. free-hand zeta and z grids (tests/flexgrid.deck): 3 ranks agree with the serial run
#   5. hybrid: 2 ranks x 2 OpenMP threads (backend = openmp) agree with the serial run, also with
#      field ionization (tests/ionization_N.deck) and a beam-dominated undulator run
# usage: sh tests/mpi.sh <gempic_rz> <compare_columns> <mpiexec> [extra mpiexec flags]
set -e
abspath() { (cd "$(dirname "$1")" && echo "$(pwd)/$(basename "$1")"); }
EXE=$(abspath "$1"); CMP=$(abspath "$2"); MPIEXEC=${3:-mpirun}; shift 3 2>/dev/null || true
FLAGS="$*"
if [ "$(id -u)" = 0 ]; then  # Open MPI refuses to run as root without these (other MPIs ignore them)
  export OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
fi
# Open MPI: allow more ranks than cores (the test uses up to 4 ranks)
if $MPIEXEC --version 2>&1 | grep -qiE "open mpi|openrte|open-mpi"; then FLAGS="$FLAGS --oversubscribe"; fi
run() { np=$1; shift; $MPIEXEC $FLAGS -n $np "$EXE" "$@"; }  # -n: mpirun and srun (srun reads -np as "-n p")
SRC=$(cd "$(dirname "$0")/.." && pwd)
WORK=$(mktemp -d)
cd "$WORK"
DECK="$SRC/tests/regression.deck"
C="numerics.backend=serial output.every=1 output.map_stride=8 grid.z_end=160k0^-1"
compare() {  # <dir> <tolerances>
  for f in axis_00002.txt axis_00004.txt map_00004.txt bunches.txt; do "$CMP" ser/$f $1/$f $2 $3; done
}
"$EXE" "$DECK" $C output.dir=ser > ser.log
GEMPIC_PIPELINE=1 "$EXE" "$DECK" $C output.dir=p1 > p1.log
compare p1 0 0
for np in 2 3 4; do
  run $np "$DECK" $C output.dir=mp$np > mp$np.log
  compare mp$np 1e-6 1e-9
done
run 3 "$DECK" $C output.dir=ck grid.z_end=80k0^-1 output.checkpoint_every=2 > ck1.log
run 2 "$DECK" $C output.dir=ck output.checkpoint_every=2 --restart latest > ck2.log
compare ck 1e-6 1e-9
cp -r ck ck_serial && "$EXE" "$DECK" $C output.dir=ck_serial --restart ck_serial/checkpoint_00004.h5 grid.z_end=200k0^-1 > cks.log
# free-hand zeta and z grids (tests/flexgrid.deck) on 3 ranks
F="numerics.backend=serial output.map_stride=8"
"$EXE" "$SRC/tests/flexgrid.deck" $F output.dir=fser > fser.log
run 3 "$SRC/tests/flexgrid.deck" $F output.dir=fmp3 > fmp3.log
for f in axis_00004.txt map_00004.txt bunches.txt; do "$CMP" fser/$f fmp3/$f 1e-6 1e-9; done
# hybrid MPI + OpenMP threads
H="numerics.backend=openmp output.every=1 output.map_stride=8 grid.z_end=160k0^-1"
OMP_NUM_THREADS=2 run 2 "$DECK" $H output.dir=hy > hy.log
grep -q "openmp on 2 MPI ranks, 2 threads per rank" hy.log
compare hy 1e-6 1e-9
N="output.every=1 grid.z_end=6k0^-1"
"$EXE" "$SRC/tests/ionization_N.deck" $N numerics.backend=serial output.dir=nser > nser.log
OMP_NUM_THREADS=2 run 2 "$SRC/tests/ionization_N.deck" $N numerics.backend=openmp output.dir=nhy > nhy.log
for f in ionization.txt axis_00006.txt; do "$CMP" nser/$f nhy/$f 1e-6 1e-12; done
U="grid.z_end=10mm bunch.macroparticles=20000"
"$EXE" "$SRC/tests/undulator.deck" $U numerics.backend=serial output.dir=user output.format=text > user.log
OMP_NUM_THREADS=2 run 2 "$SRC/tests/undulator.deck" $U numerics.backend=openmp output.dir=uhy output.format=text > uhy.log
for f in bunches.txt fel.txt; do "$CMP" user/$f uhy/$f 1e-6 1e-12; done
rc=0
run 3 "$DECK" $C output.dir=w output.walltime_hours=1e-5 > w1.log || rc=$?
if [ $rc -ne 3 ]; then echo "expected exit code 3 at the wall-time limit, got $rc"; exit 1; fi
run 3 "$DECK" $C output.dir=w --restart latest > w2.log
compare w 1e-6 1e-9
cd - > /dev/null
rm -rf "$WORK"
echo "mpi test passed"
