#!/bin/bash
# Correctness of the GPU code on real GPUs (1 node, 4 GPUs, ~1 h). Submit from the GEM-PIC-RZ directory:
#   sbatch --account=<project> benchmarks/jupiter/check.sbatch
# 1. nvidia-smi; ctest of build-gpu: the single-GPU tests (regression_cuda, physics_cuda, restart_cuda, ...) and
#    mpi_pipeline, whose device section marches the MPI blocks on the GPUs (1-3 ranks, mean-field push,
#    free-hand grids) against the serial run
# 2. the MPI pipeline on 1, 2, 4 GPUs (numerics.backend = cuda) against the serial run (build-cpu):
#    regression deck (laser wake, trapping, ions), FEL deck (4 kA beam across the blocks, mean-field push),
#    flying focus through the absorber layer; tolerance 1e-6 relative, as in the tests
# 3. physics: the 4 kA beam in a planar undulator at lambda_u/16 on 4 GPUs: eps_nx at 5 cm 0.0067 um (serial)
#SBATCH --job-name=gempic-check
#SBATCH --partition=booster
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --gpus-per-node=4
#SBATCH --cpus-per-task=16
#SBATCH --time=01:30:00
#SBATCH --output=gempic-check.%j.log
set -u
. "${SLURM_SUBMIT_DIR:?submit from the GEM-PIC-RZ directory}/benchmarks/jupiter/env.sh"
load_modules || exit 1
export OMP_NUM_THREADS=$CPUS_PER_TASK OMP_PROC_BIND=close OMP_PLACES=cores
G=$BUILD_GPU/gempic_rz; C=$BUILD_CPU/gempic_rz; CMP=$BUILD_CPU/compare_columns
W=$RUNS/check.$SLURM_JOB_ID; mkdir -p "$W"; cd "$W" || exit 1
srun_n() { np=$1; shift; srun --ntasks=$np --gpus-per-task=1 --cpus-per-task=$CPUS_PER_TASK "$@"; }
fails=0
cmpf() {  # <reference dir> <dir> <files...>
  ref=$1; d=$2; shift 2
  for f in "$@"; do
    out=$("$CMP" "$ref/$f" "$d/$f" 1e-6 1e-9 2>&1); rc=$?
    echo "$out" | tail -1
    [ $rc = 0 ] || fails=$((fails + 1))
  done
}

echo "== 1. GPUs and ctest"
nvidia-smi --query-gpu=index,name,memory.total,driver_version --format=csv,noheader
(cd "$BUILD_GPU" && ctest --output-on-failure > "$W/ctest.log" 2>&1) || fails=$((fails + 1))
grep -E "tests passed|Failed|\*\*\*|device march" "$W/ctest.log" "$BUILD_GPU/Testing/Temporary/LastTest.log" | tail -20

echo "== 2. MPI pipeline on the GPUs against the serial run"
D=$SRC/tests/regression.deck
O="output.every=1 output.map_stride=8 grid.z_end=160k0^-1"
"$C" "$D" $O numerics.backend=serial output.dir=reg_ser > reg_ser.log 2>&1
for np in 1 2 4; do
  if [ $np = 1 ]; then export GEMPIC_PIPELINE=1; fi  # one rank: the pipelined step (the block march)
  srun_n $np "$G" "$D" $O numerics.backend=cuda output.dir=reg_gpu$np > reg_gpu$np.log 2>&1
  unset GEMPIC_PIPELINE
  echo "-- regression deck, $np GPU(s): $(grep '^done' reg_gpu$np.log)"
  cmpf reg_ser reg_gpu$np axis_00002.txt axis_00004.txt map_00004.txt bunches.txt
done
FD=$SRC/benchmarks/harmonics/fel.deck
FO="grid.z_end=1cm grid.dzeta=1.328nm bunch.macroparticles=24000 output.format=text output.every=100000"
"$C" "$FD" $FO numerics.backend=serial output.dir=fel_ser > fel_ser.log 2>&1
srun_n 4 "$G" "$FD" $FO numerics.backend=cuda output.dir=fel_gpu4 > fel_gpu4.log 2>&1
echo "-- FEL deck (mean-field push, beam across the blocks), 4 GPUs: $(grep '^done' fel_gpu4.log)"
cmpf fel_ser fel_gpu4 bunches.txt fel.txt
V=$SRC/benchmarks/flying_focus/vacuum.deck
VO="grid.z_end=50um grid.absorber_strength=1 output.every=0 output.maps=false"
"$C" "$V" $VO numerics.backend=serial output.dir=ff_ser > ff_ser.log 2>&1
srun_n 4 "$G" "$V" $VO numerics.backend=cuda output.dir=ff_gpu4 > ff_gpu4.log 2>&1
echo "-- flying focus, 4 GPUs: $(grep '^done' ff_gpu4.log)"
cmpf ff_ser ff_gpu4 "$(cd ff_ser && ls axis_*.txt | tail -1)"

echo "== 3. planar undulator at lambda_u/16, 4 GPUs (serial: eps_nx = 0.0067 um at 5 cm)"
PO="grid.z_end=5cm grid.dzeta=1.328nm grid.dz=0.625mm numerics.beam_dz_max=0.625mm bunch.macroparticles=128000 output.format=text output.every=100000"
srun_n 4 "$G" "$FD" $PO numerics.backend=cuda output.dir=planar_gpu4 > planar_gpu4.log 2>&1
echo "-- $(grep '^done' planar_gpu4.log)"
echo "-- eps_nx at 5 cm: $(awk '$3 == "b" { e = $10 } END { print e }' planar_gpu4/bunches.txt) um"

echo "== summary: $fails failed comparison(s); outputs in $W"
