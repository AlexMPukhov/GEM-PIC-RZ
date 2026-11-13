#!/bin/sh
# Submit the weak-scaling series (weak.sbatch): 1, 4, 16, 64 GPUs, 1200 slices per GPU.
# Run from the GEM-PIC-RZ directory after build.sh:   sh benchmarks/jupiter/submit_weak.sh [GPU counts]
. "$(dirname "$0")/env.sh"
[ "$ACCOUNT" = CHANGE_ME ] && { echo "set ACCOUNT (env.sh or ACCOUNT=<project> sh $0)"; exit 1; }
cd "$SRC" || exit 1
mkdir -p "$RUNS"
for n in ${*:-1 4 16 64}; do
  nodes=$(( (n + GPUS_PER_NODE - 1) / GPUS_PER_NODE ))
  tpn=$(( n < GPUS_PER_NODE ? n : GPUS_PER_NODE ))
  sbatch --account="$ACCOUNT" --partition="$PARTITION" --nodes=$nodes --ntasks-per-node=$tpn \
    --gpus-per-node="$GPUS_PER_NODE" --cpus-per-task="$CPUS_PER_TASK" --time=01:00:00 \
    --export=ALL,NGPU=$n benchmarks/jupiter/weak.sbatch
done
echo "results accumulate in $RUNS/weak.txt; table: sh benchmarks/jupiter/collect.sh"
