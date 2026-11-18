---
name: gempic-rz-setup
description: Install, build and test the plasma simulation code GEM-PIC-RZ (Galilean electromagnetic particle-in-cell, quasi-cylindrical RZ geometry, C++17) and choose how to run it - serial, OpenMP threads, MPI ranks or an NVIDIA GPU. Use when gempic_rz is not installed yet, when a build or ctest fails, when the user asks how to run GEM-PIC-RZ on several cores, a cluster or a GPU, or where its output files go. To write or debug an input deck use the skill gempic-rz-simulation.
license: BSD-3-Clause
compatibility: Needs git, cmake >= 3.14 and a C++17 compiler (Linux, macOS, or WSL on Windows). Optional - MPI, the HDF5 C library, CUDA toolkit 12 or 13, python3 with numpy and h5py.
metadata:
  repository: https://github.com/AlexMPukhov/GEM-PIC-RZ
---

# GEM-PIC-RZ: build, test, run

GEM-PIC-RZ is one executable, `gempic_rz`, that reads one ASCII input deck. There is nothing to install besides the
build. Tell the user what you are going to do before you clone or build.

## 1. Is it there already?

```
gempic_rz --version            # or <repo>/build/gempic_rz --version
```

If this prints `GEM-PIC-RZ <version>`, skip to section 4 (use the path of that executable wherever this skill writes
`build/gempic_rz`). If the user has a clone, build it in place (section 2,
without the clone).

## 2. Clone and build

```
git clone https://github.com/AlexMPukhov/GEM-PIC-RZ.git
cd GEM-PIC-RZ
cmake -S . -B build && cmake --build build -j
build/gempic_rz --version
```

- **Required:** a C++17 compiler (g++ >= 9 or clang) and cmake >= 3.14. Without cmake: `make` (plain Makefile, same
  result in `build/`).
- **Found automatically, optional:** OpenMP (threads), MPI (any implementation, `mpicxx`), the HDF5 C library
  (`output.format = hdf5 / openpmd`, checkpoints and restart, LASY laser files, bunches from files). On Debian/Ubuntu:
  `apt install g++ cmake libhdf5-dev libopenmpi-dev`. Read the cmake output to see what was found; a missing HDF5 is
  not an error, but then every HDF5 feature of a deck is refused at run time.
- **Windows:** build and run inside WSL (Ubuntu); the commands are the same.
- **Portable binary:** `-DGEMPIC_NATIVE=OFF` builds without `-march=native` (needed when the build machine and the
  compute nodes have different CPUs).
- **Python (only for some tests, the GUI and plotting):** python3 with numpy, h5py, matplotlib.

## 3. Test

```
cd build && ctest --output-on-failure          # about 35 tests, ~20 min on 2 cores
ctest -R "deck_parser|regression|agent_skill"   # quick check, about 1 min
```

- A test that reports `Skipped` lacks an optional tool (python3 with numpy/h5py, MPI, HDF5); that is not a failure.
- In containers and WSL set `OMP_PROC_BIND=false`. MPI tests as root need
  `OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1`. Under Slurm run the MPI tests with `--ntasks >= 4`.
- If a test fails on an unmodified clone, report the test name and its output to the user; do not edit the tests.

## 4. Run

```
build/gempic_rz my.deck --check                   # read the deck, print summary and warnings, exit (always do this first)
build/gempic_rz my.deck                           # run, one thread (backend serial)
OMP_NUM_THREADS=8 build/gempic_rz my.deck numerics.backend=openmp   # CPU threads
mpirun -np 16 build/gempic_rz my.deck             # MPI: 16 blocks of slices along zeta
build/gempic_rz my.deck grid.z_end=2mm output.dir=scan_a   # any deck value can be overridden: block.key=value
build/gempic_rz my.deck --restart latest          # continue from the newest checkpoint (needs output.checkpoint_every)
build/gempic_rz --template > my.deck              # annotated starting deck
```

| backend | set by | use it for |
|---|---|---|
| `serial` | default of a CPU build | reference; small tests |
| `openmp` | `numerics.backend=openmp`, threads from `OMP_NUM_THREADS` | a workstation: parallel over the particles and the radial grid of a slice |
| MPI | `mpirun -np P` (with `serial` or `openmp` on every rank) | long windows: P blocks of slices, pipelined over the steps; needs >= 4 slices per rank and many more steps than ranks |
| `cuda` | a GPU build, default there | slices with many particles (dense loading, high modes); one GPU per MPI rank |

GPU build (NVIDIA; only `nvcc` is needed, `tools/install_cuda_pip.sh` installs it without root; details in the
repository's `docs/gpu.md`):

```
cmake -S . -B build-gpu -DGEMPIC_CUDA=ON -DCMAKE_CUDA_ARCHITECTURES=native && cmake --build build-gpu -j
```

All backends run the same algorithm and read the same deck; `openmp` and MPI reproduce `serial` (to rounding with
`-march=native`), `cuda` to rounding.

## 5. Where the output goes

Everything is written to `output.dir` (default `diag`), relative to the directory where `gempic_rz` is started:
`summary.txt`, `deck_used.txt`, `bunches.txt` (beam moments, one line per step), and every `output.every` steps
`axis_NNNNN.txt`, `map_NNNNN.txt`, `beam_NNNNN.txt` (text) or `data_NNNNN.h5` / `openpmd/openpmd_<step>.h5`. The
first output step prints `output volume: ... GB per output step ...`; read it before a long run, and do not put
`output.dir` on tmpfs (`/dev/shm`). The last console line gives the peak memory.

Other tools in the build: `build/h5_extract <file.h5> axis|lineout|energy|spectrum ...` (tables from the HDF5 files),
`python3 tools/gempic_gui.py my.deck` (deck editor, run control and live plots in the browser),
`gnuplot -e "step=50; dir='diag'" tools/plot_maps.gp` (maps of a text output step as PDF).

## 6. Cluster jobs

Never run simulations on a login node. Build there, then submit; with a wall-time limit let the run stop itself and
continue:

```
# job.sh -- deck with  output { checkpoint_every = 100; walltime_hours = 23.5 }
R=""; if ls diag/checkpoint_*.h5 > /dev/null 2>&1; then R="--restart latest"; fi
srun build/gempic_rz run.deck $R
if [ $? -eq 3 ]; then sbatch job.sh; fi       # exit code 3: time limit reached, checkpoint written
```

## Citing

A. Pukhov, *GEM-PIC-RZ*, arXiv:2610.31622 (2026); method: A. Pukhov, N. Elkina, T. Wilson, J. Comput. Phys. 567,
115331 (2026); code: doi:10.5281/zenodo.23600178. Do not describe the code as "quasi-static": it keeps the slow
derivative of the forward field.
