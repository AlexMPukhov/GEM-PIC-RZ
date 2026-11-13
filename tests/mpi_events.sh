#!/bin/sh
# Events of the MPI pipeline that act on the whole beam (synchronization points, pipeline.cpp: drain, syncEvent):
#   bunch.z_start: the witness of tests/regression.deck loaded at the end of step 2 (z = 80 k0^-1) in the plasma
#   wake. Reference: the serial run (Simulation::loadBunches). Checked: the pipelined step on one rank = serial
#   (to all digits), 2-4 ranks = serial (rounding: 1e-6), also with beam_field = mean, with a checkpoint before
#   the loading (3 ranks, restart on 2) and one after it (restart: the bunch is not loaded again), and on the
#   device (GPU emulation build) on 3 ranks.
#   chicane: the chirped flat-top bunch of tests/transport.sh through a thin chicane at z = 10 cm (R56 = 100 um,
#   delay 10 um: the beam moves 10 um back, across the block boundaries of 2-4 ranks, and is compressed by 2).
#   Reference: the serial run (Simulation::applyChicanes). Checked: 2-4 ranks = serial (bunches.txt every step,
#   axis fields before, at and after the chicane), also with beam_field = mean, a checkpoint before and one after
#   the chicane (3 -> 2 ranks), and on the device (GPU emulation) on 3 ranks.
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
# ---- chicane (tests/transport.sh: 80 slices, 20 per rank on 4 ranks)
sed -e 's/^bunch .*/bunch { name = b; charge = 1 pC; sigma_x = 100 um; emittance_x = 0 um; distribution_z = function; profile_z = "1"; zeta_min = 28 um; zeta_max = 32 um/'     -e 's/^        sigma_z = 2 um.*/        gamma_function = "1000 * (1 - 0.01 * (zeta - 30 um) \/ (2 um))"; macroparticles = 20000 }/'     "$SRC/tests/transport.deck" > c.deck
echo 'chicane { z = 10 cm; r56 = 100 um; delay = 10 um; energy = 510.49 MeV }' >> c.deck
grep -q gamma_function c.deck
CO="numerics.backend=serial output.every=1 output.axis=true grid.z_end=0.11m"
ccmp() {  # <reference dir> <dir>
  for f in axis_00019.txt axis_00020.txt axis_00022.txt bunches.txt; do "$CMP" $1/$f $2/$f 1e-6 1e-9; done
}
"$EXE" c.deck $CO output.dir=cser > cser.log
grep -q "beam moved by -10 um" cser.log
for np in 2 3 4; do
  run $np c.deck $CO output.dir=cmp$np > cmp$np.log
  grep -q "beam moved by -10 um" cmp$np.log
  ccmp cser cmp$np
done
"$EXE" c.deck $CO numerics.beam_field=mean output.dir=cmser > cmser.log
run 3 c.deck $CO numerics.beam_field=mean output.dir=cmmp > cmmp.log
ccmp cmser cmmp
echo "  chicane: 2-4 ranks = serial (the beam moves across the block boundaries), also beam_field = mean"
run 3 c.deck $CO output.dir=cc1 grid.z_end=0.095m output.checkpoint_every=19 > cc1a.log
run 2 c.deck $CO output.dir=cc1 --restart latest > cc1b.log
grep -q "beam moved by -10 um" cc1b.log
ccmp cser cc1
run 3 c.deck $CO output.dir=cc2 grid.z_end=0.105m output.checkpoint_every=21 > cc2a.log
run 2 c.deck $CO output.dir=cc2 --restart latest > cc2b.log
if grep -q "beam moved" cc2b.log; then echo "FAIL: the chicane was applied again after the restart"; exit 1; fi
"$CMP" cser/axis_00022.txt cc2/axis_00022.txt 1e-6 1e-9
echo "  chicane across checkpoints (3 -> 2 ranks): before and after the chicane = serial"
if [ -x "$(dirname "$EXE")/gempic_rz_gpuemu" ]; then
  CE="$(echo $CO | sed s/backend=serial/backend=openmp/)"
  GEMPIC_PIPELINE_ENGINE=1 OMP_NUM_THREADS=2 runx 3 "$(dirname "$EXE")/gempic_rz_gpuemu" c.deck $CE output.dir=cdv > cdv.log
  ccmp cser cdv
  echo "  chicane on the device (GPU emulation, 3 ranks) = serial"
fi
cd - > /dev/null
rm -rf "$WORK"
echo "mpi events test passed"
