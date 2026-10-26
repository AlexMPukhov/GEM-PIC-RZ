#!/bin/bash
# Strong scaling of the MPI pipeline on GPUs (one rank per GPU). Submitted by submit_scaling.sh with
# NGPU = 1, 4, 16, 64 (and the matching --nodes); every job appends its timings to $RUNS/scaling.txt.
#   lwfa: benchmarks/gpu/lwfa_heavy.deck (LWFA, ~4e4 plasma particles per slice, 1200 slices, 200 steps)
#   fel:  benchmarks/harmonics/fel.deck (4 kA beam, 640 000 particles, 3 modes, 18 372 slices; 400 steps)
# With NGPU = 1 also the single-GPU engine without the pipeline (the reference of the speed-up).
#SBATCH --job-name=gempic-scale
#SBATCH --output=gempic-scale.%j.log
set -u
. "${SLURM_SUBMIT_DIR:?submit from the GEM-PIC-RZ directory}/benchmarks/jupiter/env.sh"
load_modules || exit 1
NGPU=${NGPU:?set NGPU (submit_scaling.sh does)}
EXTRA=${EXTRA:-}  # more overrides for all runs, e.g. EXTRA=grid.z_end=20um for a short trial
export OMP_NUM_THREADS=$CPUS_PER_TASK OMP_PROC_BIND=close OMP_PLACES=cores GEMPIC_PROFILE=1
G=$BUILD_GPU/gempic_rz
W=$RUNS/scale_${NGPU}gpu.$SLURM_JOB_ID; mkdir -p "$W"; cd "$W" || exit 1
nvidia-smi --query-gpu=name --format=csv,noheader | head -1 > gpu.txt
run() {  # <label> <ranks> <deck> <overrides...>
  label=$1; np=$2; deck=$3; shift 3
  t0=$(date +%s)
  srun --ntasks=$np --gpus-per-task=1 --cpus-per-task=$CPUS_PER_TASK "$G" "$deck" "$@" $EXTRA numerics.backend=cuda \
    output.dir=$label > $label.log 2>&1
  rc=$?
  done=$(grep '^done' $label.log | head -1)
  # "done: <steps> steps in <seconds> s, ..."
  steps=$(echo "$done" | awk '{print $2}'); sec=$(echo "$done" | awk '{print $5}')
  per=$(grep 'per step after the pipeline filled' $label.log | awk '{print $2}')  # (GEMPIC_PROFILE=1, rank 0)
  echo "$NGPU $label $np ${steps:-NA} ${sec:-NA} $(( $(date +%s) - t0 )) rc=$rc $(cat gpu.txt | tr ' ' '_') ${per:-NA}" \
    | tee -a "$RUNS/scaling.txt"
}
L=$SRC/benchmarks/gpu/lwfa_heavy.deck
LO="output.every=0 output.maps=false"
F=$SRC/benchmarks/harmonics/fel.deck
FO="grid.z_end=25cm output.every=0"
if [ "$NGPU" = 1 ]; then
  run lwfa_engine 1 "$L" $LO          # the single-GPU engine (Engine::step, whole window)
  run fel_engine 1 "$F" $FO
  export GEMPIC_PIPELINE=1            # the pipelined step on one rank (the block march)
fi
run lwfa $NGPU "$L" $LO
run fel $NGPU "$F" $FO
