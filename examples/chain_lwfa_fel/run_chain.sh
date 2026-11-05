#!/bin/sh
# The chained example (README.md): stage 1 (FBPIC, down-ramp injection) -> handoff.py -> stage 2 (GEM-PIC-RZ,
# acceleration) -> stage 3 (GEM-PIC-RZ, plasma lens, chicane, undulator). Prints the beam after each import and the
# power at the resonant wavelength along the undulator.
# usage: sh run_chain.sh                 (in a work directory; the decks are taken from this script's directory)
# environment: GEMPIC = gempic_rz executable (default: gempic_rz in PATH), PYTHON (default python3, with numpy and
#              h5py; stage 1 needs FBPIC), OMP_NUM_THREADS = threads (numerics.backend = openmp).
# Stage 1 is skipped if stage1/hdf5/data00006000.h5 exists (e.g. computed on a GPU machine and copied here).
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
GEMPIC=${GEMPIC:-gempic_rz}
PYTHON=${PYTHON:-python3}
export OMP_PROC_BIND=${OMP_PROC_BIND:-false}
RUN="$GEMPIC"
S1=stage1/hdf5/data00006000.h5
if [ ! -f "$S1" ]; then
  echo "== stage 1: FBPIC, density down-ramp injection (0.2 mm)"
  $PYTHON "$HERE/stage1_fbpic.py" 200 > stage1.log 2>&1 || { tail -20 stage1.log; echo "stage 1 failed (needs FBPIC; see README.md)"; exit 1; }
  tail -1 stage1.log
fi
echo "== hand-off: the injected beam (gamma > 20) behind the laser"
$PYTHON "$HERE/handoff.py" "$S1" stage1/laser_centre.txt 20 30.56 | tee handoff.log
ZC=$(tail -1 handoff.log)
echo "== stage 2: GEM-PIC-RZ, acceleration to 1.4 mm"
$RUN "$HERE/stage2_accel.deck" numerics.backend=openmp bunch.file="$S1" "$ZC" > stage2.log
grep -A1 "^bunch\[0\]" stage2.log
grep "^done" stage2.log
awk '$3 == "witness" { l = $0 } END { print "  at the exit:", l }' stage2/bunches.txt
echo "== stage 3: GEM-PIC-RZ, plasma lens, chicane, undulator (80 periods)"
$RUN "$HERE/stage3_fel.deck" numerics.backend=openmp > stage3.log
grep -A1 "^bunch\[0\]" stage3.log
grep "^chicane at z = 40 mm: beam\|^done" stage3.log
awk '!/^#/ && $1 % 160 == 0 { printf "  z = %5.2f m   P(lambda_r) = %.3g W   bunching %.3g\n", $2 / 1000, $7, $9 }' stage3/fel.txt
