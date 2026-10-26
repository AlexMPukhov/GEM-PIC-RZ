#!/bin/bash
# Weak scaling of the MPI pipeline on GPUs (one rank per GPU): the LWFA deck with a window of 1200 slices per
# GPU (benchmarks/gpu/lwfa_heavy.deck, nzeta = 1200 NGPU, the laser kept at the same distance from the front of
# the window), 200 steps. Submitted by submit_weak.sh with NGPU = 1, 4, 16, 64; every job appends to
# $RUNS/weak.txt: GPUs, run, ranks, steps, wall time, the time per step after the pipeline has filled (rank 0,
# GEMPIC_PROFILE=1) and the exit code. Ideal: the same time per step for all GPU counts.
#SBATCH --job-name=gempic-weak
#SBATCH --output=gempic-weak.%j.log
set -u
. "${SLURM_SUBMIT_DIR:?submit from the GEM-PIC-RZ directory}/benchmarks/jupiter/env.sh"
load_modules || exit 1
NGPU=${NGPU:?set NGPU (submit_weak.sh does)}
EXTRA=${EXTRA:-}
export OMP_NUM_THREADS=$CPUS_PER_TASK OMP_PROC_BIND=close OMP_PLACES=cores GEMPIC_PROFILE=1
[ "$NGPU" = 1 ] && export GEMPIC_PIPELINE=1  # the pipelined step also on one rank
G=$BUILD_GPU/gempic_rz
W=$RUNS/weak_${NGPU}gpu.$SLURM_JOB_ID; mkdir -p "$W"; cd "$W" || exit 1
nvidia-smi --query-gpu=name --format=csv,noheader | head -1 > gpu.txt
NZ=$(( 1200 * NGPU ))
ZC=$(( 240 + 300 * (NGPU - 1) ))  # k0^-1: 60 k0^-1 behind the front, as in the deck
srun --ntasks=$NGPU --gpus-per-task=1 --cpus-per-task=$CPUS_PER_TASK "$G" "$SRC/benchmarks/gpu/lwfa_heavy.deck" \
  grid.nzeta=$NZ "laser.zeta_center=${ZC}k0^-1" output.every=0 output.maps=false $EXTRA numerics.backend=cuda \
  output.dir=lwfa > lwfa.log 2>&1
rc=$?
done=$(grep '^done' lwfa.log | head -1)
steps=$(echo "$done" | awk '{print $2}'); sec=$(echo "$done" | awk '{print $5}')
per=$(grep 'per step after the pipeline filled' lwfa.log | awk '{print $2}')
echo "$NGPU lwfa_weak $NGPU ${steps:-NA} ${sec:-NA} ${per:-NA} rc=$rc $(tr ' ' '_' < gpu.txt)" | tee -a "$RUNS/weak.txt"
