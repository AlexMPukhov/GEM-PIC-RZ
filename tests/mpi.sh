#!/bin/sh
# MPI pipeline along zeta (tests/regression.deck, 4 steps):
#   1. the pipelined step on one rank (GEMPIC_PIPELINE=1) is identical to the serial step
#   2. 2, 3 and 4 ranks agree with the serial run (to rounding: 1e-6 relative)
#   3. checkpoint on 3 ranks, restart on 2 ranks / serially; wall-time stop on 3 ranks
#   4. free-hand zeta and z grids (tests/flexgrid.deck): 3 ranks agree with the serial run
#   5. hybrid: 2 ranks x 2 OpenMP threads (backend = openmp) agree with the serial run, also with
#      field ionization (tests/ionization_N.deck) and a beam-dominated undulator run
#   6. (python3 with numpy and h5py) a laser from a LASY file fed through the absorber layer
#      (laser.inject = boundary) and the LASY export (output.lasy): 3 ranks = serial
#   7. mean-field beam push (numerics.beam_field = mean): the pipelined step on one rank is identical to the
#      serial step, 2-4 ranks agree with it, also across a checkpoint (regression deck: wake, trapping,
#      absorber; beam filter 2); the intense beam in a planar undulator at dz = lambda_u/16
#      (benchmarks/harmonics/fel.deck, the beam spans the 3 blocks): 3 ranks = serial
#   8. the march of the blocks on the device (Engine::marchBlock; GPU emulation build, or a CUDA build with a
#      GPU): 1-3 ranks, the mean-field push and the free-hand grids agree with the serial run
#   9. metal wall (tests/wall.deck): steep narrowing and widening tubes on 3 ranks (also on the device) = serial
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
runx() { np=$1; x=$2; shift 2; $MPIEXEC $FLAGS -n $np "$x" "$@"; }  # with another executable
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
# mean-field beam push: second push of step k at the start of step k+1 (meanPhase in pipeline.cpp)
MN="$C numerics.beam_field=mean numerics.beam_filter=2"
"$EXE" "$DECK" $MN output.dir=mser > mser.log
grep -q "beam field mean" mser.log
GEMPIC_PIPELINE=1 "$EXE" "$DECK" $MN output.dir=mp1 > mp1.log
for f in axis_00002.txt axis_00004.txt map_00004.txt bunches.txt; do "$CMP" mser/$f mp1/$f 0 0; done
for np in 2 3 4; do
  run $np "$DECK" $MN output.dir=mmp$np > mmp$np.log
  for f in axis_00002.txt axis_00004.txt map_00004.txt bunches.txt; do "$CMP" mser/$f mmp$np/$f 1e-6 1e-9; done
done
run 3 "$DECK" $MN output.dir=mck grid.z_end=80k0^-1 output.checkpoint_every=2 > mck1.log
run 2 "$DECK" $MN output.dir=mck output.checkpoint_every=2 --restart latest > mck2.log
for f in axis_00004.txt map_00004.txt bunches.txt; do "$CMP" mser/$f mck/$f 1e-6 1e-9; done
FD="$SRC/benchmarks/harmonics/fel.deck"
FO="numerics.backend=serial grid.z_end=1cm grid.dzeta=1.328nm bunch.macroparticles=24000 output.format=text output.every=100000"
"$EXE" "$FD" $FO output.dir=fel_ser > fel_ser.log
run 3 "$FD" $FO output.dir=fel_mp3 > fel_mp3.log
for f in fel.txt bunches.txt; do "$CMP" fel_ser/$f fel_mp3/$f 1e-6 1e-12; done
echo "  mean-field beam push: 1 rank pipelined = serial, 2-4 ranks and a 3 -> 2 rank restart agree"
# free-hand zeta and z grids (tests/flexgrid.deck) on 3 ranks
F="numerics.backend=serial output.map_stride=8"
"$EXE" "$SRC/tests/flexgrid.deck" $F output.dir=fser > fser.log
run 3 "$SRC/tests/flexgrid.deck" $F output.dir=fmp3 > fmp3.log
for f in axis_00004.txt map_00004.txt bunches.txt; do "$CMP" fser/$f fmp3/$f 1e-6 1e-9; done
# hybrid MPI + OpenMP threads (not in a CUDA build: there the parallel backend is cuda, which is not combined
# with MPI; the pure-MPI parts above cover the pipeline)
if "$EXE" "$DECK" numerics.backend=openmp --check > hycheck.log 2>&1; then
H="numerics.backend=openmp output.every=1 output.map_stride=8 grid.z_end=160k0^-1"
OMP_NUM_THREADS=2 run 2 "$DECK" $H output.dir=hy > hy.log
grep -q "openmp on 2 MPI ranks, 2 threads per rank" hy.log
compare hy 1e-6 1e-9
N="output.every=1 grid.z_end=6k0^-1"
"$EXE" "$SRC/tests/ionization_N.deck" $N numerics.backend=serial output.dir=nser > nser.log
OMP_NUM_THREADS=2 run 2 "$SRC/tests/ionization_N.deck" $N numerics.backend=openmp output.dir=nhy > nhy.log
for f in ionization.txt axis_00006.txt; do "$CMP" nser/$f nhy/$f 1e-6 1e-12; done
U="grid.z_end=10mm bunch.macroparticles=20000"  # numerics.beam_field = mean (the default with undulators)
"$EXE" "$SRC/tests/undulator.deck" $U numerics.backend=serial output.dir=user output.format=text > user.log
OMP_NUM_THREADS=2 run 2 "$SRC/tests/undulator.deck" $U numerics.backend=openmp output.dir=uhy output.format=text > uhy.log
for f in bunches.txt fel.txt; do "$CMP" user/$f uhy/$f 1e-6 1e-12; done
else
  echo "  hybrid MPI + OpenMP part skipped: no openmp backend in this build (CUDA build)"
fi
# the march of the blocks on the device (Engine::marchBlock): the GPU code paths emulated on CPU threads
# (gempic_rz_gpuemu of this build, GEMPIC_PIPELINE_ENGINE=1), or a CUDA build on its GPU(s)
DEV=""
if [ -x "$(dirname "$EXE")/gempic_rz_gpuemu" ]; then
  DEV="$(dirname "$EXE")/gempic_rz_gpuemu"; DB="numerics.backend=openmp"; export GEMPIC_PIPELINE_ENGINE=1
elif "$EXE" "$DECK" numerics.backend=cuda --check > devcheck.log 2>&1 && command -v nvidia-smi > /dev/null; then
  DEV="$EXE"; DB="numerics.backend=cuda"
fi
if [ -n "$DEV" ]; then
  E="output.every=1 output.map_stride=8 grid.z_end=160k0^-1 $DB"
  GEMPIC_PIPELINE=1 "$DEV" "$DECK" $E output.dir=d1 > d1.log
  compare d1 1e-6 1e-9
  for np in 2 3; do
    OMP_NUM_THREADS=2 runx $np "$DEV" "$DECK" $E output.dir=d$np > d$np.log
    compare d$np 1e-6 1e-9
  done
  OMP_NUM_THREADS=2 runx 3 "$DEV" "$DECK" $E numerics.beam_field=mean numerics.beam_filter=2 output.dir=dm > dm.log
  for f in axis_00002.txt axis_00004.txt map_00004.txt bunches.txt; do "$CMP" mser/$f dm/$f 1e-6 1e-9; done
  OMP_NUM_THREADS=2 runx 3 "$DEV" "$SRC/tests/flexgrid.deck" $DB output.map_stride=8 output.dir=dfx > dfx.log
  for f in axis_00004.txt map_00004.txt bunches.txt; do "$CMP" fser/$f dfx/$f 1e-6 1e-9; done
  unset GEMPIC_PIPELINE_ENGINE
  echo "  device march of the blocks ($DB): 1-3 ranks, mean-field push, free-hand grids = serial"
else
  echo "  device march of the blocks skipped: no GPU emulation build and no GPU"
fi
# metal wall (tests/wall.deck): steep narrowing and widening tubes (the wall crosses a face every ~4 steps), the
# sloped-wall E_z and the field continued into the opened annulus: 3 ranks = serial, also on the device
for wl in "12 9" "9 12"; do
  sed "s/^wall.*/wall { z_points = 0 0.3 mm; r_points = $wl um }/" "$SRC/tests/wall.deck" > wall.deck
  WL="grid.z_end=0.3mm laser.waist=3.5um output.every=60 numerics.backend=serial"
  "$EXE" wall.deck $WL output.dir=wser > wser.log
  run 3 wall.deck $WL output.dir=wmp > wmp.log
  "$CMP" wser/map_00060.txt wmp/map_00060.txt 1e-6 1e-9
  if [ -n "$DEV" ]; then
    GEMPIC_PIPELINE_ENGINE=1 OMP_NUM_THREADS=2 runx 3 "$DEV" wall.deck $WL $DB output.dir=wdv > wdv.log
    "$CMP" wser/map_00060.txt wdv/map_00060.txt 1e-6 1e-9
  fi
done
echo "  metal wall, narrowing and widening tube: 3 ranks = serial$([ -n "$DEV" ] && echo ", also on the device")"
rc=0
run 3 "$DECK" $C output.dir=w output.walltime_hours=1e-5 > w1.log || rc=$?
if [ $rc -ne 3 ]; then echo "expected exit code 3 at the wall-time limit, got $rc"; exit 1; fi
run 3 "$DECK" $C output.dir=w --restart latest > w2.log
compare w 1e-6 1e-9
if python3 -c "import numpy, h5py" 2>/dev/null; then
  python3 "$SRC/tests/lasy_files.py"
  sed -e 's|^flying_focus.*|laser { file = ffm.h5; zeta_center = 90 um; inject = boundary }|' \
    "$SRC/benchmarks/flying_focus/vacuum.deck" > ffl.deck
  L="grid.z_end=50um grid.absorber_strength=1 output.every=0 output.maps=false numerics.backend=serial output.lasy=true"
  "$EXE" ffl.deck $L output.dir=lser > lser.log
  run 3 ffl.deck $L output.dir=lmp > lmp.log
  "$CMP" lser/axis_00005.txt lmp/axis_00005.txt 1e-6 1e-9
  python3 -c "
import h5py, numpy as np
a, b = (h5py.File(d + '/lasy/laser_00005.h5')['/data/5/meshes/laserEnvelope'][()] for d in ('lser', 'lmp'))
assert np.abs(a - b).max() <= 1e-6 * np.abs(a).max()"
  echo "  LASY boundary injection and export: 3 ranks = serial"
fi
cd - > /dev/null
rm -rf "$WORK"
echo "mpi test passed"
