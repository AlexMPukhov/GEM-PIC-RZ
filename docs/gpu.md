# Running GEM-PIC-RZ on an NVIDIA GPU

GEM-PIC-RZ does **not** use Kokkos. The GPU kernels are written against a small
portability layer of its own (`src/portable.hpp`). The only compiler needed is
`nvcc` from the standard NVIDIA CUDA Toolkit (version 12 or 13). No
`nvcc_wrapper`, Kokkos build or other special compiler is involved, and the
rest of the code is compiled by your usual g++.

## 1. Check the machine

```
nvidia-smi            # GPU model, driver, and "CUDA Version: X.Y" = newest CUDA the driver supports
which nvcc            # is a CUDA toolkit already installed?
module avail cuda     # on clusters: CUDA toolkits as modules
```

The toolkit must not be newer than the "CUDA Version" that `nvidia-smi` shows.

| GPU | compute capability (`CMAKE_CUDA_ARCHITECTURES`) |
|---|---|
| V100 | 70 (CUDA 12 only; CUDA 13 dropped it) |
| T4, RTX 20xx | 75 |
| A100, A30 | 80 |
| RTX 30xx, A40, A10 | 86 |
| L40, RTX 40xx | 89 |
| H100, H200 | 90 |

## 2. Get `nvcc` (no root rights needed)

Choose one of the following.

**a) Cluster module.** Load the CUDA module, e.g. `module load cuda` (or `CUDA/12.4`, whatever `module avail` lists).

**b) Driver ≥ 580 (nvidia-smi shows CUDA 13.x): pip.** Run

```
sh tools/install_cuda_pip.sh $HOME/cuda-pip
```

This installs nvcc, the CUDA runtime and the headers of CUDA 13.4 from NVIDIA's packages on PyPI into `$HOME/cuda-pip`, in about 15 s, and prints the build command. It was tested with GEM-PIC-RZ: a clean build succeeds.

**c) Older driver (nvidia-smi shows CUDA 12.x): CUDA 12 toolkit in your home directory.** Download the Linux **runfile** of a CUDA 12 version not newer than your driver from https://developer.nvidia.com/cuda-toolkit-archive, then run

```
sh cuda_12.*_linux.run --silent --toolkit --toolkitpath=$HOME/cuda-12
```

`--toolkit` installs only the compiler and libraries, not a driver. Alternatively, with conda:

```
conda create -p $HOME/cuda-12 -c nvidia cuda-nvcc=12.6 cuda-cudart-dev=12.6
```

(Routes a, c and the conda variant could not be tried here; route b was.)

## 3. Build and test

```
cmake -S . -B build-gpu -DGEMPIC_CUDA=ON -DCMAKE_CUDA_ARCHITECTURES=native \
      -DCMAKE_CUDA_COMPILER=$HOME/cuda-pip/cuda/bin/nvcc
cmake --build build-gpu -j
cd build-gpu && ctest                                               # all tests, including the GPU ones
```

- **Paths:** use the path of your toolkit's `nvcc`. With a module or a system-wide toolkit, drop `-DCMAKE_CUDA_COMPILER`.
- **CUDA runtime:** it is linked statically, so the program runs without `LD_LIBRARY_PATH`. Only the NVIDIA driver must be installed.
- **Old CMake:** `native` needs CMake ≥ 3.24. With an older CMake, give the number from the table, e.g. `-DCMAKE_CUDA_ARCHITECTURES=80`.
- **Makefile:** `make CUDA=1 NVCC=/path/to/nvcc CUDA_ARCH="-arch=sm_80"`.

The GPU tests compare the GPU with the serial reference:

| test | checks |
|---|---|
| `regression_cuda` | the regression deck against the serial reference (to 1e-6) |
| `physics_cuda` | bunch field (Gauss), ion collapse, field ionization of H and N |
| `restart_cuda` | checkpoint and restart on the GPU |
| `mpi_pipeline` (MPI builds) | the ζ-march of the blocks on the GPU (`numerics.backend = cuda` with 1–3 ranks, mean-field push, free-hand grids) against the serial run |

Several GPUs (MPI, one rank per GPU): `mpirun -np 4 build-gpu/gempic_rz run.deck numerics.backend=cuda`. Every rank
takes the GPU of its node-local rank (Open MPI, Slurm, MPICH, MVAPICH variables), or the single GPU it sees with
`srun --gpus-per-task=1`. Each rank marches its block of slices on its GPU; the beam push and the exchanges go
through the host (see the MPI section of the README).

Then a timing run:

```
GEMPIC_PROFILE=1 build-gpu/gempic_rz examples/ionization_injection.deck grid.z_end=0.1mm
build/gempic_rz examples/ionization_injection.deck grid.z_end=0.1mm numerics.backend=serial   # CPU for comparison
```

The profile line at the end shows where the GPU time goes (push, field solve, launch overhead).

## 4. Problems when configuring

- **`Could NOT find Threads` (from `FindCUDAToolkit.cmake`):** this happens with GEM-PIC-RZ versions before 2026-09-30. The build no longer calls `find_package(CUDAToolkit)`, so this check is gone.
- **`nvcc fatal: Value 'gnu++17' is not defined for option 'std'` while "Building CXX object":** nvcc was used as the C++ compiler, typically because `CXX=nvcc` is set in the environment (as in Kokkos setups). Since 2026-09-30 the build detects this: g++ compiles the `.cpp` files, and that nvcc becomes the CUDA compiler. Delete the old build directory first, because it has nvcc cached as the C++ compiler: `rm -rf build-gpu`, then configure again with the full command.
- **Other configure or link errors:** send `build-gpu/CMakeFiles/CMakeConfigureLog.yaml` (CMake ≥ 3.26) or `CMakeFiles/CMakeError.log`, plus the full error text.

## 5. Messages at the start of a GPU run

| message | meaning / remedy |
|---|---|
| `no NVIDIA driver found` | not a GPU machine, or no driver: use `backend = serial` / `openmp` there |
| `the NVIDIA driver supports CUDA up to X, but this binary was built with CUDA Y` | toolkit newer than the driver: build with a toolkit ≤ X (route c), or update the driver |
| `no usable CUDA device` | driver present but no GPU visible (check `CUDA_VISIBLE_DEVICES`, the job's GPU request) |
| `the binary has no code for this GPU` | rebuild with `-DCMAKE_CUDA_ARCHITECTURES=native` or the number from the table |
| `error while loading shared libraries: libcudart.so...` | only if you set `-DCMAKE_CUDA_RUNTIME_LIBRARY=Shared`: add the toolkit's `lib` directory to `LD_LIBRARY_PATH` |

## 6. Without a GPU

Everything except the GPU run itself can be checked on any machine:

- **`backend = openmp`:** the same kernels on CPU threads.
- **`gempic_rz_gpuemu`:** built in CPU builds. It runs the GPU-only code paths (thread teams with barriers, atomic deposition, unordered appends) on CPU threads. The ctest `regression_gpu_emulation` checks it against the serial reference.
