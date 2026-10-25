# Settings of the JUPITER / JUWELS Booster runs (sourced by the other scripts). Edit ACCOUNT first.
# Values in the environment take precedence: e.g.  ACCOUNT=myproj sh build.sh
ACCOUNT=${ACCOUNT:-CHANGE_ME}            # compute project (sbatch --account)
PARTITION=${PARTITION:-booster}          # GPU partition (JUPITER and JUWELS Booster: booster; tests: develbooster)
CPU_PARTITION=${CPU_PARTITION:-$PARTITION}  # partition of the CPU reference (the same nodes, GPUs unused)
GPUS_PER_NODE=${GPUS_PER_NODE:-4}        # 4 GH200 per JUPITER node, 4 A100 per JUWELS Booster node
CPUS_PER_TASK=${CPUS_PER_TASK:-16}       # host cores per rank (the beam bookkeeping, I/O, MPI)
CUDA_ARCH=${CUDA_ARCH:-90}               # 90: GH200 (JUPITER), 80: A100 (JUWELS Booster)
MODULES=${MODULES:-"GCC OpenMPI CUDA HDF5 CMake"}  # module names as on the JSC systems (module spider if absent)
# the GEM-PIC-RZ checkout: batch jobs are submitted from it (Slurm runs a copy of the job script elsewhere)
if [ -z "${SRC:-}" ]; then
  if [ -n "${SLURM_SUBMIT_DIR:-}" ] && [ -f "$SLURM_SUBMIT_DIR/src/engine.cpp" ]; then SRC=$SLURM_SUBMIT_DIR
  else SRC=$(cd "$(dirname "$0")/../.." && pwd); fi
fi
BUILD_GPU=${BUILD_GPU:-$SRC/build-gpu}   # CUDA build (numerics.backend = cuda)
BUILD_CPU=${BUILD_CPU:-$SRC/build-cpu}   # CPU build (serial references, CPU timings)
RUNS=${RUNS:-$SRC/benchmarks/jupiter/runs}  # output of all jobs

load_modules() {
  if command -v module > /dev/null 2>&1; then
    # shellcheck disable=SC2086
    module load $MODULES || { echo "module load $MODULES failed: set MODULES (see 'module spider CUDA')"; return 1; }
  fi
}
