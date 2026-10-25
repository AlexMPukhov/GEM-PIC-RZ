#!/bin/sh
# Build GEM-PIC-RZ for the JUPITER / JUWELS Booster runs, on a login node of the same system (JUPITER: aarch64
# Grace + GH200; JUWELS Booster: x86 + A100):
#   build-gpu  CUDA build (numerics.backend = cuda), MPI, HDF5; ctest runs the MPI tests with srun
#   build-cpu  CPU build (serial references, CPU timing on one node)
# usage: sh benchmarks/jupiter/build.sh           (settings in env.sh)
set -e
. "$(dirname "$0")/env.sh"
load_modules
command -v nvcc > /dev/null || { echo "nvcc not found: load the CUDA module (MODULES in env.sh)"; exit 1; }
SRUN=$(command -v srun || true)
echo "nvcc: $(nvcc --version | tail -1)"
echo "MPI launcher for ctest: ${SRUN:-mpiexec}"
# -DGEMPIC_NATIVE=OFF: the login nodes need not have the CPU of the compute nodes
cmake -S "$SRC" -B "$BUILD_GPU" -DCMAKE_BUILD_TYPE=Release -DGEMPIC_NATIVE=OFF \
  -DGEMPIC_CUDA=ON -DCMAKE_CUDA_ARCHITECTURES="$CUDA_ARCH" \
  ${SRUN:+-DMPIEXEC_EXECUTABLE=$SRUN -DMPIEXEC_PREFLAGS=--cpus-per-task=$CPUS_PER_TASK}
cmake --build "$BUILD_GPU" -j 16
cmake -S "$SRC" -B "$BUILD_CPU" -DCMAKE_BUILD_TYPE=Release -DGEMPIC_NATIVE=OFF \
  ${SRUN:+-DMPIEXEC_EXECUTABLE=$SRUN -DMPIEXEC_PREFLAGS=--cpus-per-task=$CPUS_PER_TASK}
cmake --build "$BUILD_CPU" -j 16
echo
echo "built: $BUILD_GPU/gempic_rz (GPU), $BUILD_CPU/gempic_rz (CPU)"
echo "next:  sbatch --account=\$ACCOUNT benchmarks/jupiter/check.sbatch   (see benchmarks/jupiter/README.md)"
