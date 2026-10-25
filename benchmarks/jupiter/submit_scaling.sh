#!/bin/sh
# Submit the scaling series: 1, 4, 16, 64 GPUs (one rank per GPU) and the CPU reference on one node.
# Run from the GEM-PIC-RZ directory after build.sh:   sh benchmarks/jupiter/submit_scaling.sh [GPU counts]
# e.g.  sh benchmarks/jupiter/submit_scaling.sh 1 4      (only the first two)
. "$(dirname "$0")/env.sh"
[ "$ACCOUNT" = CHANGE_ME ] && { echo "set ACCOUNT (env.sh or ACCOUNT=<project> sh $0)"; exit 1; }
cd "$SRC" || exit 1
mkdir -p "$RUNS"
for n in ${*:-1 4 16 64}; do
  nodes=$(( (n + GPUS_PER_NODE - 1) / GPUS_PER_NODE ))
  tpn=$(( n < GPUS_PER_NODE ? n : GPUS_PER_NODE ))
  # time: the FEL deck on one GPU dominates (400 steps of 18 372 slices); fewer GPUs, more time
  case $n in 1) t=04:00:00 ;; 4) t=02:00:00 ;; *) t=01:00:00 ;; esac
  sbatch --account="$ACCOUNT" --partition="$PARTITION" --nodes=$nodes --ntasks-per-node=$tpn \
    --gpus-per-node="$GPUS_PER_NODE" --cpus-per-task="$CPUS_PER_TASK" --time=$t \
    --export=ALL,NGPU=$n benchmarks/jupiter/scaling.sbatch
done
sbatch --account="$ACCOUNT" --partition="$CPU_PARTITION" --time=04:00:00 benchmarks/jupiter/cpu_reference.sbatch
echo "results accumulate in $RUNS/scaling.txt; table: sh benchmarks/jupiter/collect.sh"
