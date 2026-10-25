#!/bin/bash
# GPU vs serial on ONE GPU without Slurm and MPI (a workstation): the comparisons of
# benchmarks/jupiter/check.sbatch for one rank, plus wall times GPU vs CPU threads.
#   sh benchmarks/gpu/check_1gpu.sh [build-gpu] [build-cpu] [work dir]
# defaults: build-gpu, build-cpu (CPU build with OpenMP) and runs-1gpu in the GEM-PIC-RZ directory.
# THREADS (default: all cores) = OpenMP threads of the CPU timing runs.
# 1. regression deck (laser wake, trapping, ions): the single-GPU engine (whole window) and the pipelined step
#    on one rank (GEMPIC_PIPELINE=1: the block march of the MPI code) against serial, 1e-6 relative
# 2. FEL deck (4 kA beam, mean-field push) and flying focus through the absorber layer: pipelined GPU vs serial
# 3. 4 kA beam in a planar undulator at lambda_u/16, pipelined on the GPU: eps_nx at 5 cm (serial 0.0067 um)
# 4. wall times: lwfa_heavy.deck (plasma) and the FEL deck (beam, 96 steps), GPU vs CPU with THREADS threads
set -u
SRC=$(cd "$(dirname "$0")/../.." && pwd)
BG=$(cd "${1:-$SRC/build-gpu}" && pwd) || exit 1
BC=$(cd "${2:-$SRC/build-cpu}" && pwd) || exit 1
W=${3:-$SRC/runs-1gpu}; mkdir -p "$W"; cd "$W" || exit 1
G=$BG/gempic_rz; C=$BC/gempic_rz; CMP=$BC/compare_columns
THREADS=${THREADS:-$(nproc)}
fails=0
cmpf() {  # <reference dir> <dir> <files...>
  ref=$1; d=$2; shift 2
  for f in "$@"; do
    out=$("$CMP" "$ref/$f" "$d/$f" 1e-6 1e-9 2>&1); rc=$?
    echo "   $f: $(echo "$out" | tail -1)"
    [ $rc = 0 ] || fails=$((fails + 1))
  done
}
donel() { grep '^done' "$1" | sed 's/, particles.*//'; }

echo "== GPU: $(nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv,noheader | head -1)"
echo "== 1. regression deck: GPU against serial"
D=$SRC/tests/regression.deck
O="output.every=1 output.map_stride=8 grid.z_end=160k0^-1"
"$C" "$D" $O numerics.backend=serial output.dir=reg_ser > reg_ser.log 2>&1
"$G" "$D" $O numerics.backend=cuda output.dir=reg_eng > reg_eng.log 2>&1
echo "-- single-GPU engine: $(donel reg_eng.log)"
cmpf reg_ser reg_eng axis_00002.txt axis_00004.txt map_00004.txt bunches.txt
GEMPIC_PIPELINE=1 "$G" "$D" $O numerics.backend=cuda output.dir=reg_pip > reg_pip.log 2>&1
echo "-- pipelined step (block march), 1 GPU: $(donel reg_pip.log)"
cmpf reg_ser reg_pip axis_00002.txt axis_00004.txt map_00004.txt bunches.txt

echo "== 2. FEL deck and flying focus, pipelined on the GPU, against serial"
FD=$SRC/benchmarks/harmonics/fel.deck
FO="grid.z_end=1cm grid.dzeta=1.328nm bunch.macroparticles=24000 output.format=text output.every=100000"
"$C" "$FD" $FO numerics.backend=serial output.dir=fel_ser > fel_ser.log 2>&1
GEMPIC_PIPELINE=1 "$G" "$FD" $FO numerics.backend=cuda output.dir=fel_gpu > fel_gpu.log 2>&1
echo "-- FEL deck (mean-field push): $(donel fel_gpu.log)"
cmpf fel_ser fel_gpu bunches.txt fel.txt
V=$SRC/benchmarks/flying_focus/vacuum.deck
VO="grid.z_end=50um grid.absorber_strength=1 output.every=0 output.maps=false"
"$C" "$V" $VO numerics.backend=serial output.dir=ff_ser > ff_ser.log 2>&1
GEMPIC_PIPELINE=1 "$G" "$V" $VO numerics.backend=cuda output.dir=ff_gpu > ff_gpu.log 2>&1
echo "-- flying focus: $(donel ff_gpu.log)"
cmpf ff_ser ff_gpu "$(cd ff_ser && ls axis_*.txt | tail -1)"

echo "== 3. planar undulator at lambda_u/16, pipelined on the GPU (serial: eps_nx = 0.0067 um at 5 cm)"
PO="grid.z_end=5cm grid.dzeta=1.328nm grid.dz=0.625mm numerics.beam_dz_max=0.625mm bunch.macroparticles=128000 output.format=text output.every=100000"
GEMPIC_PIPELINE=1 "$G" "$FD" $PO numerics.backend=cuda output.dir=planar_gpu > planar_gpu.log 2>&1
echo "-- $(donel planar_gpu.log)"
echo "-- eps_nx at 5 cm: $(awk '$3 == "b" { e = $10 } END { print e }' planar_gpu/bunches.txt) um"

echo "== 4. wall times, GPU against $THREADS CPU threads (OpenMP)"
L=$SRC/benchmarks/gpu/lwfa_heavy.deck
LO="output.every=0 output.maps=false"
F96="grid.z_end=30mm output.every=0"
"$G" "$L" $LO numerics.backend=cuda output.dir=t_lwfa_gpu > t_lwfa_gpu.log 2>&1
OMP_NUM_THREADS=$THREADS OMP_PROC_BIND=close "$C" "$L" $LO numerics.backend=openmp output.dir=t_lwfa_cpu \
  > t_lwfa_cpu.log 2>&1
echo "-- lwfa_heavy, GPU:          $(donel t_lwfa_gpu.log)"
echo "-- lwfa_heavy, $THREADS threads:   $(donel t_lwfa_cpu.log)"
"$G" "$FD" $F96 numerics.backend=cuda output.dir=t_fel_gpu > t_fel_gpu.log 2>&1
OMP_NUM_THREADS=$THREADS OMP_PROC_BIND=close "$C" "$FD" $F96 numerics.backend=openmp output.dir=t_fel_cpu \
  > t_fel_cpu.log 2>&1
echo "-- FEL 96 steps, GPU:        $(donel t_fel_gpu.log)"
echo "-- FEL 96 steps, $THREADS threads: $(donel t_fel_cpu.log)"
for d in t_fel_gpu t_fel_cpu; do
  awk -v d=$d '!/^#/ && $1 == 96 && $3 == 1 && $4 == 1 {print "   " d ": P1 in slice 14 at step 96:", $19, "W"}' \
    $d/fel_harmonics.txt
done

echo "== summary: $fails failed comparison(s); outputs in $W"
