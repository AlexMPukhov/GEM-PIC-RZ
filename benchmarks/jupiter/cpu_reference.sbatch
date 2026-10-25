#!/bin/bash
# CPU reference for the scaling table: the same two decks on ONE node without GPUs, MPI pipeline with OpenMP
# threads (numerics.backend = openmp), RANKS ranks x THREADS threads. Submitted by submit_scaling.sh.
# JUPITER node (4 Grace CPUs, 4 x 72 cores): 4 x 72; JUWELS Booster node (2 x 24 cores): RANKS=4 THREADS=12.
#SBATCH --job-name=gempic-cpu
#SBATCH --nodes=1
#SBATCH --output=gempic-cpu.%j.log
set -u
. "${SLURM_SUBMIT_DIR:?submit from the GEM-PIC-RZ directory}/benchmarks/jupiter/env.sh"
load_modules || exit 1
RANKS=${RANKS:-4}; THREADS=${THREADS:-72}
EXTRA=${EXTRA:-}  # more overrides for all runs (as in scaling.sbatch)
export OMP_NUM_THREADS=$THREADS OMP_PROC_BIND=close OMP_PLACES=cores
C=$BUILD_CPU/gempic_rz
W=$RUNS/cpu_${RANKS}x${THREADS}.$SLURM_JOB_ID; mkdir -p "$W"; cd "$W" || exit 1
run() {  # <label> <deck> <overrides...>
  label=$1; deck=$2; shift 2
  t0=$(date +%s)
  srun --ntasks=$RANKS --cpus-per-task=$THREADS --cpu-bind=cores "$C" "$deck" "$@" $EXTRA numerics.backend=openmp \
    output.dir=$label > $label.log 2>&1
  rc=$?
  done=$(grep '^done' $label.log | head -1)
  steps=$(echo "$done" | awk '{print $2}'); sec=$(echo "$done" | awk '{print $5}')
  echo "cpu ${label}_cpu ${RANKS}x${THREADS} ${steps:-NA} ${sec:-NA} $(( $(date +%s) - t0 )) rc=$rc node" \
    | tee -a "$RUNS/scaling.txt"
}
run lwfa "$SRC/benchmarks/gpu/lwfa_heavy.deck" output.every=0 output.maps=false
run fel "$SRC/benchmarks/harmonics/fel.deck" grid.z_end=25cm output.every=0
