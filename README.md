# GEM-PIC-RZ

Galilean electromagnetic particle-in-cell code in quasi-cylindrical geometry
(azimuthal modes), with free-hand non-uniform grids along r, ζ and z: the
step is given at reference points and interpolated linearly in between. The
scheme is described in a manuscript in preparation. This project grew out of
a single-file prototype (not included). Two discretization errors of the
prototype have been corrected (see "Corrections"), so results with particle
bunches differ from the prototype.

## Build

```
cmake -S . -B build && cmake --build build -j     # or simply: make
cd build && ctest                                 # or: make test
```

C++17 is required. MPI (any implementation) and the HDF5 C library are optional (HDF5 output; e.g. `apt install libhdf5-dev`); there are no other external libraries. `-DGEMPIC_NATIVE=OFF` builds without `-march=native`.

GPU build (NVIDIA, CUDA toolkit 12 or 13; only plain `nvcc` is needed, no Kokkos. See
[`docs/gpu.md`](docs/gpu.md) for installing nvcc without root rights, e.g.
`sh tools/install_cuda_pip.sh`):

```
cmake -S . -B build-gpu -DGEMPIC_CUDA=ON -DCMAKE_CUDA_ARCHITECTURES=native   # or 80 (A100), 90 (H100), ...
cmake --build build-gpu -j
cd build-gpu && ctest                         # includes regression_cuda and physics_cuda on the GPU
make CUDA=1 CUDA_ARCH="-arch=sm_80"           # the same with the plain Makefile
```

## Run

```
build/gempic_rz examples/lwfa_witness.deck
build/gempic_rz examples/lwfa_witness.deck --check              # check the deck, print the summary
build/gempic_rz examples/lwfa_witness.deck laser.a0=2.5 output.dir=a0_2.5
build/gempic_rz --template > my.deck                            # annotated starting point
OMP_NUM_THREADS=8 build/gempic_rz my.deck numerics.backend=openmp  # CPU threads
build-gpu/gempic_rz my.deck                                     # GPU (default backend of a GPU build)
mpirun -np 16 build/gempic_rz my.deck                           # MPI: 16 ranks along zeta
gnuplot -e "step=50; dir='lwfa_witness'" tools/plot_maps.gp     # maps_00050.pdf
build/gempic_rz my.deck output.format=hdf5                      # results in data_NNNNN.h5
build/h5_extract diag/data_00050.h5 lineout 2.0 90              # E, B at r = 2 um, phi = 90 deg
build/h5_extract diag/data_00050.h5 energy                      # field energy and its zeta-Nyquist part
build/gempic_rz my.deck output.checkpoint_every=100             # checkpoints for restarts
build/gempic_rz my.deck --restart latest                        # continue from the newest one
```

## Input deck

The input is a structured ASCII deck with one block per kind of object:

```
constants{ ... }    named values for expressions       (optional)
grid     { ... }    gridding, modes, window, steps     (once)
laser    { ... }    one block per laser pulse          (0 ... N)
plasma   { ... }    background plasma profile          (0 or 1)
ion      { ... }    one block per mobile ion species   (0 ... N)
bunch    { ... }    one block per particle bunch       (0 ... N)
undulator{ ... }    one block per undulator section    (0 ... N)
output   { ... }    diagnostics                        (optional)
numerics { ... }    algorithm parameters               (optional)
```

- **Units:** values are given in physical units (`um`, `fs`, `cm^-3`, `pC`, `MeV`, …) or in normalized code units (`k0^-1`, `n_c`). A dimensional value without a unit is an error.
- **Error messages:** unknown keys, a unit of the wrong kind, and missing or conflicting keys are reported with file and line number.
- **Functions:** profiles can be given as expressions, as in HiPACE++ and QUARZ: the plasma density
  `density_function(x, y, z, r)`, the laser `profile_function(r)` and `envelope_function(zeta)`, the bunch current
  `profile_z(zeta)` and energy chirp `gamma_function(zeta)`, and an undulator taper `K_function(z)`. Expressions accept
  numbers with units (`r / (20 um)`), the usual functions and operators, and the named values of the block
  `constants`; any single value of the deck may be an expression of them (`waist = "2 * w0"`). See
  [Functions and expressions](docs/input_deck.md#functions-and-expressions).
- **Reference:** all keys, units and defaults are listed in [`docs/input_deck.md`](docs/input_deck.md).

| example | shows |
|---|---|
| `examples/lwfa_witness.deck` | LWFA with an external witness bunch on a free-hand radial grid (the prototype's default case) |
| `examples/fel_sase.deck` | compact high-gain FEL (SASE from the macroparticle noise): gain length against Ming Xie |
| `examples/undulator_emission.deck` | coherent spontaneous emission of a short bunch in a planar undulator (λ_r resolved on the ζ grid) |
| `examples/lwfa_flexgrid.deck` | the same on free-hand grids along ζ and z (fewer slices and steps, small steps in the up-ramp) |
| `examples/function_profiles.deck` | profiles as functions: density up-ramp and channel, super-Gaussian laser with a sin² envelope, ramped and chirped witness |
| `examples/pwfa_driver_witness.deck` | beam-driven wake with drive and witness bunches, flat-top witness, m = 0 only |
| `examples/two_pulses_channel.deck` | two laser pulses (focused drive pulse, circularly polarized 2nd harmonic), parabolic channel, trapezoidal profile |
| `examples/ion_collapse.deck` | mobile hydrogen ions collapsing inside a dense electron bunch (Rosenzweig et al. 2005) |
| `examples/ionization_injection.deck` | ionization injection: N⁵⁺ dopant field-ionized to N⁶⁺/N⁷⁺ by an a0 = 2 pulse (ADK, NIST energies) |

## Source layout

| file | content |
|---|---|
| `src/deck.*` | generic parser: blocks, `key = value`, lists, strings, comments, line numbers, command-line overrides |
| `src/units.*` | unit table and conversion to code units |
| `src/beam_import.*` | reading particle beams from openPMD and native HDF5 files |
| `src/expr.*` | function parser: expressions of the deck (profiles, constants) compiled to a stack program |
| `src/config.*` | typed structures `GridConfig`, `LaserConfig`, `PlasmaConfig`, `BunchConfig`, `OutputConfig`, `NumericsConfig`; reading, validation, summary |
| `src/atomic_data.*` | NIST ionization energies of all charge states of 21 elements, standard atomic masses |
| `src/ionization.*` | ADK tunnel-ionization rates per level (SI reference function and code-unit tables) |
| `src/mapping.*` | free-hand grid mapping: step given at reference points, linear in between (analytic, used for r, ζ and z) |
| `src/grid.*` | staggered non-uniform radial grid, interpolation weights |
| `src/external.hpp` | external fields acting on the particles: undulators, quadrupoles, plasma lenses, solenoids |
| `src/eigen.*` | symmetric tridiagonal eigenproblem (initial fields of bunches with v < c) |
| `src/zgrid.hpp` | longitudinal grid of the window: slices, half-cells, location of particles, interpolation weights of the R equation |
| `src/fields.hpp` | mode fields and currents |
| `src/solver.*` | tridiagonal field update per slice and mode (matrices per distinct half-cell length and step) |
| `src/particles.hpp` | particle record, Boris pusher |
| `src/simulation.*` | laser and bunch initialization, plasma loading, pushers, deposition, time step (serial reference) |
| `src/engine.*` | parallel engine (GPU / OpenMP): the time step as kernels, device-resident state |
| `src/pipeline.*` | MPI: domain decomposition along ζ, pipelined in the steps; diagnostics and output on rank 0 |
| `src/comm.*` | thin MPI layer (packed messages; a single-rank stub without MPI) |
| `src/portable.hpp`, `src/hd.hpp` | portability layer: memory, kernel launchers, teams, atomics, counter-based RNG, complex numbers |
| `src/diagnostics.cpp` | console output, `bunches.txt`, axis/map/particle files |
| `src/output_hdf5.cpp` | HDF5 output (`data_NNNNN.h5`), openPMD output, checkpoints (`checkpoint_NNNNN.h5`) and restart |
| `src/main.cpp` | command-line driver |
| `tests/test_deck.cpp` | about 230 checks (the number depends on the build options): parser, units, configuration, backends, output format, free-hand mapping and ζ/z grids, ion and ionization blocks, error messages, deposition shape, NIST table, ADK = Landau rate; all example decks must be valid |
| `tests/regression.*` | runs `tests/regression.deck` (laser, plasma, mobile H ions, witness) and compares with `tests/reference/`, for every backend |
| `tests/physics.*`, `tests/check_physics.cpp` | bunch self-field in vacuum against Gauss's law (stretched and free-hand radial grid, free-hand ζ grid); a γ = 100 bunch keeps its energy in its own field; ion collapse against a test-ion model; field ionization of H and N⁵⁺→N⁷⁺ against the ADK rate equations |
| `tools/h5_extract.cpp` | reads the HDF5 files: tables of axis, map, beam, field line-outs (also an example reader) |
| `tests/hdf5.sh` | HDF5 files read back against the text output (double and single precision) |
| `tests/flexgrid.sh`, `tests/flexgrid.deck` | free-hand ζ and z grids: constant-step tables = uniform grid; parallel backend = serial; restart, also with a new z grid |
| `tests/undulator.sh`, `tests/undulator.deck` | undulator radiation: energy, spectrum and slippage against the far-field radiation integral (`tools/undulator_radint.py`) |
| `tests/cherenkov.sh` | numerical Cherenkov: ζ-Nyquist field energy of a relativistic bunch in vacuum and in plasma, with and without the beam filter |
| `tests/mpi.sh` | MPI: one-rank pipeline = serial; 2–4 ranks against serial; restarts across rank counts; wall-time stop |
| `tests/restart.sh` | interrupted + restarted run against the uninterrupted run (checkpoint file, wall-time stop) |
| `tools/install_cuda_pip.sh` | installs nvcc + CUDA runtime from PyPI into a user directory (no root) |
| `tools/compare_columns.cpp` | numeric comparison of ASCII tables (used by the regression test) |
| `tools/plot_maps.gp` | gnuplot → PDF: E_z and density maps, on-axis lineout |

## MPI: domain decomposition along ζ

`mpirun -np P gempic_rz run.deck` splits the ζ-window into P blocks of
consecutive slices. Rank P−1 is at the front, rank 0 at the rear. MPI is
found automatically by CMake and by the Makefile (`mpicxx`).

**Pipeline in the steps.** Within one step s⁰ → s¹ the ζ-march is sequential,
from the front of the window to the rear. A ζ-decomposition therefore runs as
a pipeline: rank q computes step k while rank q+1 already computes step k+1.

- **Efficiency:** about N_steps / (N_steps + P − 1) for equal blocks. The streaming plasma loads every slice equally.
- **Block size:** at least 4 slices per rank.

**Exchanges per step and rank (block [a, b)):**

| direction | when | content |
|---|---|---|
| front → rear | end of the step | slice a of the new level, the current prediction used at slice a (the R-correction needs exactly this value), beam currents deposited into the half-cells below a, the streaming particles (plasma, freed electrons, waiting particles), beam particles drifting into the next block (mid-push if necessary) |
| rear → front | after the first 3 slices | corrected R of slice b, the final slices b−2, b−1 and currents of half-cell b−1 (ghosts for the next step), currents of particles trapped at slice b−1 into half-cell b, copies of the beam particles in half-cell b−1 |

- **Halo copies:** the receiving rank pushes these copies of the rear rank's beam particles itself, a small redundant computation, so that its prediction of the beam current below its block is complete.
- **Beam ownership:** a beam particle belongs to the rank of the upper half-cell of its deposition.

**Speed.** The ionization-injection example (20 steps) runs 1.90× faster on 2 ranks / 2 cores than serially: 51.2 s vs 97.2 s, as predicted by 2 × 20/21. The results are identical.

**Hybrid MPI + OpenMP.** With `numerics.backend = openmp` every rank pushes its streaming particles (plasma,
ions, freed electrons) slice by slice and its beam particles with OpenMP threads: contiguous chunks of the particle
list per thread, per-thread current buffers of the half-cell, freed electrons and the beam bookkeeping in the serial
order. One rank with 2 threads is as fast as the OpenMP engine (regression deck: 14.8 s vs 14.9 s; serial 25.7 s).
Run e.g. `OMP_NUM_THREADS=8 mpirun -np 4 --bind-to none gempic_rz run.deck numerics.backend=openmp` on 32 cores.
`tests/mpi.sh` checks 2 ranks × 2 threads against the serial run (regression deck, field ionization, undulator).

**Output.** Rank 0 finishes every step last. It collects the diagnostics of all
ranks and writes the console line, `bunches.txt`, `ionization.txt`, the
text/HDF5 files and the checkpoints, in the same format as a serial run. A
checkpoint can be continued with any number of ranks, or serially, and vice
versa. The wall-time limit is decided by the front rank and passed down the
pipeline. All ranks then stop at the same step and write a checkpoint (exit
code 3).

**Validation (`tests/mpi.sh`, ctest `mpi_pipeline`):**

- **One rank:** the pipelined step (`GEMPIC_PIPELINE=1`) is identical to the serial step.
- **2, 3, 4 ranks:** agree with the serial run to rounding (1e-6 relative). In practice they agree to all printed digits on the regression deck, including a bunch centred on a block boundary.
- **Restarts:** checkpoint on 3 ranks → restart on 2 ranks → continuation serially; wall-time stop and `--restart latest` on 3 ranks.
- **Hard case:** the ionization-injection example at a0 = 2.5 on 3 ranks (13 steps: trapping, re-streaming of 1936 particles, field ionization). `bunches.txt` and `ionization.txt` are identical to the serial run, and the axis fields agree to rounding.

**Current limits:**

- **Backend:** `numerics.backend = serial` runs the serial step on every rank; `numerics.backend = openmp` adds OpenMP threads on every rank (hybrid MPI + threads, `OMP_NUM_THREADS` threads per rank). The GPU engine is not yet combined with MPI.
- **Memory:** every rank stores the field and current arrays of the whole window (the particles are distributed). This keeps the indexing of the serial code and is small for typical grids, but it limits very large windows.

## GPU and multi-core backends

`numerics.backend` selects the implementation of the time step:

| backend | where | build |
|---|---|---|
| `serial` | one CPU thread: the reference implementation (`simulation.cpp`) | always |
| `openmp` | CPU threads (`OMP_NUM_THREADS`) | CPU build (default) |
| `cuda` | NVIDIA GPU | `-DGEMPIC_CUDA=ON` |

The parallel engine (`src/engine.cpp`) keeps fields, currents and particles in
device memory for the whole run. Per step:

- the beam list is pushed with one thread per particle;
- the ζ-march runs slice by slice, each slice with the predictor-corrector passes;
- per pass there are three kernels:
  - streaming push, deposition and ionization, with one thread per particle, atomic current deposition, and freed electrons appended through an atomic counter;
  - field solve, with one thread block per azimuthal mode: the radial operators are applied in parallel, the tridiagonal solves run in one thread;
  - slice preparation: the predictor and the |E| bound for ionization.

There are no host-device transfers inside the march. After each step, the fields,
densities and beam list are copied back for the diagnostics.

The kernels are written once, as lambdas, behind a thin portability layer (`src/portable.hpp`). The same source compiles with nvcc for the GPU and with a C++ compiler for OpenMP threads.

**Validation:**

- **`openmp`:** reproduces `serial` bit for bit, provided both are compiled without FMA contraction (`-DGEMPIC_NATIVE=OFF`). This was checked on the ionization-injection example at a0 = 2.5 over 20 steps: field ionization, self-trapping, and 10 750 particles returned from the beam list to streaming. With `-march=native` the compiler contracts the complex arithmetic of the two codes differently, so last-bit differences appear. They stay at rounding level on the regression deck (1200 axis values equal to 7 digits, on 1 and 2 threads) but grow in chaotic runs with trapping.
- **Emulated GPU paths:** the CPU build also makes `gempic_rz_gpuemu` (ctest `regression_gpu_emulation`). It runs the GPU code paths on CPU threads: teams of 8 threads with real barriers, atomic deposition, and unordered appends to the particle lists. It agrees with `serial` to 1e-9 on the regression deck and passes all physics tests.
- **Ionization:** the `openmp` and emulated runs give identical results, so the counter-based random numbers do not depend on the thread schedule.
- **Real GPU:** run `ctest` in the GPU build. `regression_cuda` compares the GPU run with the serial reference (1e-6 relative); `physics_cuda` runs the bunch-field, ion-collapse and ionization benchmarks on the GPU.

**Performance notes:**

- **Where the time goes:** the march is sequential in ζ, so the GPU gains from the particles within a slice. A slice with a few thousand particles, as in the examples, is launch-latency bound (about 10 kernel launches per slice).
- **Profiling:** `GEMPIC_PROFILE=1` prints the time per phase.
- **OpenMP:** it scales with the threads (1.8× on 2 cores for the regression deck). Each thread deposits into a private copy of the current half-cell, so no atomics are needed there.

## HDF5 output

`output { format = hdf5 }` (or `both`) writes one file `data_NNNNN.h5` per
output step. Each file holds:

- every azimuthal mode of E and B on the native staggered grid;
- the radial and ζ coordinates;
- the electron and ion densities;
- the beam list: bunches and trapped particles, with species ids;
- the input deck.

All data are in code units, with a `unitSI` attribute on every dataset. The
layout is listed in `docs/input_deck.md`.

- **Size:** `hdf5_precision = single` and `hdf5_compression = 1..9` reduce the size (26.5 MB → 7.6 MB for the regression deck with single precision and level 6).
- **Build:** HDF5 is optional: CMake and the Makefile enable it when the HDF5 C library is found. Without it, `format = hdf5` is rejected when the deck is read.
- **Reading, C++:** `tools/h5_extract.cpp` is an example reader. It prints `info`, `axis`, `map`, `beam`, a `lineout` at any (r, φ), the field `energy` of the window (optionally of a ζ range, per mode group) together with its ζ-Nyquist part (a monitor of numerical Cherenkov radiation), or the `spectrum` dW/dk of the radiation (modes m ≥ 1). It reconstructs the fields from the modes exactly as the particle gather does, including the axis.
- **Reading, Python:** h5py reads the complex fields directly, e.g. `h5py.File('data_00050.h5')['fields/Ez'][1]` for mode 1.
- **Test:** `hdf5_output` checks that the files, read back with `h5_extract`, reproduce the text output value by value (1200 axis values, 60 000 map values, 20 000 beam particles).

## openPMD output

`output { format = openpmd }` (combinable: `format = text, openpmd`) writes files that follow the openPMD standard
1.1.0, `<dir>/openpmd/openpmd_<step>.h5`, readable with openPMD-viewer, openPMD-api, VisualPIC or ParaView. The fields
are written in `thetaMode` on a uniform r–ζ grid (`openpmd_dr`, `openpmd_dzeta`; interpolated as in the particle
gather), the lab position is z = ζ + s and the time s/c, and each bunch is a particle species with SI units. The
files pass the openPMD validator without errors or warnings; details in
[`docs/input_deck.md`](docs/input_deck.md#openpmd-files).

## Checkpoints and restart

`output { checkpoint_every = 100 }` writes `checkpoint_NNNNN.h5` every 100
steps and at the end of the run; `checkpoint_keep` (default 2) limits how many
are kept. `walltime_hours = 23.5` stops a cluster job with a checkpoint before
its time limit (exit code 3).

`gempic_rz run.deck --restart latest` (or `--restart <file>`) continues the
run. The continuation is exact: the checkpoint holds the complete state in
double precision, including the random-generator state. The deck may change
`z_end`, `dz` (or the z table), the output settings and the backend, so a GPU
run can be continued on a CPU. The r and ζ grids must stay the same, and the
restart checks them. See
`docs/input_deck.md`.

- **Tests:** `restart_serial` and `restart_openmp` (`restart_cuda` in a GPU build) compare an interrupted and restarted run with the uninterrupted one. They cover a restart from a checkpoint file and a stop at the wall-time limit followed by `--restart latest`.
- **Hard case:** the ionization-injection example at a0 = 2.5 was also restarted at step 10, after which 4256 more particles are re-streamed. All tables, 6915 beam particles and the HDF5 files at step 14 are identical to the uninterrupted run.

## Free-hand grids (r, ζ, z)

All three axes take their step from reference points, with linear
interpolation in between and a constant step outside, the same principle as
in QUARZ:

```
radial_grid  = table                     # radial step dr(r)
r_points     = 0,  1.5, 4,   24 um
dr_points    = 40, 40,  400, 400 nm
zeta_grid    = table                     # slice spacing dzeta(zeta), zeta = 0 at the rear
zeta_points  = 0,   8,   12,    20,    22.5, 25,    38.2 um
dzeta_points = 127, 127, 31.83, 31.83, 63.7, 31.83, 31.83 nm
z_grid       = table                     # step dz(z) along the propagation
z_points     = 0,   51,  150 um
dz_points    = 2.5, 2.5, 10.2 um
```

- **Mapping:** one mapping (`src/mapping.*`) for all axes. The node positions follow from dx/dk = h(x), with k the
  cell or step index. Within a segment the steps grow geometrically, and the nodes are analytic.
  The summary prints the number of cells and steps, the step range and the largest change from one cell (or step)
  to the next. A warning is issued above 25 %; keep it below ~20 %. Uniform grids (`dr`, `dzeta`, `dz`) remain the
  default and run bit for bit as before.
- **ζ:** the two-slice equations of the scheme (L, B_z, E_z) are trapezoidal rules over one half-cell and hold for
  any half-cell length h_n. The equation of R lives on the slice; its current is the linear interpolation of the
  two neighbouring half-cells, `(h_n j_{n−1/2} + h_{n−1} j_{n+1/2})/(h_{n−1} + h_n)`, i.e. the plain average on a
  uniform grid. The tridiagonal matrices `16 − h_n(2Δ + h_n) Lap` are factorized once per distinct h_n.
  Deposition divides by the local h_n, and beam particles deposit linearly between half-cell centres.
  The beam filter acts in index space, so it removes the odd–even (Nyquist) mode at any local step.
  The window is set by `window_length` (rounded up to a whole cell) or `nzeta`. Laser pulses need about 20 points
  per wavelength over their full length; a warning is issued if the step at a pulse is too large.
- **z:** the step from s is `X(ξ(s) + 1) − s`. Every step refactorizes the (cheap) tridiagonal matrices when Δ
  changes. A restart may continue with a different `dz` or `z_grid`; the ζ grid must be the one of the checkpoint.

**Checks:**
- Constant-step tables reproduce the uniform run exactly (`tests/flexgrid.sh`).
- On a non-uniform ζ and z grid, the OpenMP, GPU-emulation and MPI results equal the serial code to rounding.
- Restarts are exact, also when the z grid changes at the restart.
- The bunch self-field on a ζ grid with steps of 2.5–8 k0⁻¹ across the bunch satisfies Gauss's law to 0.2%,
  which is the shot-noise level (physics test).

| ζ grid: beam-driven wake, n = 0.002 n_c (σ_z = 1.5 k0⁻¹ driver, witness 130 k0⁻¹ behind) | slices | Ez error at the driver / witness | witness ΔE error | witness dE/E |
|---|---|---|---|---|
| uniform, h = 0.25 k0⁻¹ (reference) | 2800 | — | — | 1.932 × 10⁻⁵ |
| uniform, h = 2 k0⁻¹ | 350 | 16% / 6.6% | 1.3% | 1.267 × 10⁻⁵ |
| free-hand: 0.25 at the driver and witness, 2 elsewhere | 734 | 0.11% / 0.84% | 0.06% | 1.933 × 10⁻⁵ |

The E_z errors are relative to the largest |E_z|; ΔE is the energy change of the witness over 50 steps. The free-hand run is 3.6 times faster than the fine uniform run.

| LWFA, `examples/lwfa_witness.deck` to z = 0.2 mm | slices × steps | witness energy at 0.2 mm | run time (2 cores, OpenMP) |
|---|---|---|---|
| uniform ζ (31.8 nm), dz = 5.09 µm | 1200 × 40 | 529.854 MeV | 54 s |
| free-hand ζ (31.8–127 nm, as in `examples/lwfa_flexgrid.deck`) | 897 × 40 | 529.853 MeV | 42 s |
| uniform, dz = 2.5 µm | 1200 × 80 | 530.133 MeV | 142 s |
| free-hand z: 2.5 µm in the up-ramp → 10.2 µm | 1200 × 44 | 530.215 MeV | 68 s |

- **Free-hand ζ grid:** ahead of the back of the bubble, the fields agree with the uniform run to 10⁻⁴–10⁻³ of
  their maximum, also where the step is doubled. Behind it, the sheath crossing is chaotic; the differences there
  are the same as those between two uniform runs with different random seeds.
- **Free-hand z grid:** the energy gain is set mainly by the steps in the up-ramp. With 44 steps, the free-hand
  z grid gives an energy gain within 0.4% of the 80-step run. Uniform 5.09 µm steps are 1.4% off.

## Undulators (FEL)

`undulator` blocks add static planar or helical undulator fields that act on the particles (`src/external.hpp`,
`docs/input_deck.md`). No other change to the scheme is needed. GT2 keeps the full carrier in ζ and is paraxial in
s, so a forward wave stays at fixed ζ, and an electron slips back by λ_r = λ_u(1 + K²/2)/(2γ²) per period. The
transverse current j_⊥ = −e n v_⊥ drives R, and the radiation acts back on the particles through the gather. So
spontaneous and coherent emission, and the FEL interaction, are contained in the existing equations. The ζ grid
resolves the radiation (λ_r/20 or finer). The steps along z resolve the wiggling (λ_u/16 or finer); the summary
warns otherwise.

**Test: coherent spontaneous emission** (`examples/undulator_emission.deck`). The setup:
- planar undulator, 50 periods of 14.22 mm, K = 0.5;
- bunch of 1 fC at γ = 100, σ_z = 0.1 µm, σ_r = 100 µm, zero emittance;
- λ_r = 0.8 µm, resolved by 32 slices.

The reference is a far-field radiation integral (Jackson 14.65, `tools/undulator_radint.py`) over the same trajectory (with the ramps), times the
bunch form factor exp(−k²σ_z² − k²θ²σ_r²), integrated over 0.6 < k/k_r < 1.3 and θ < 8 mrad.

| quantity | code | radiation integral |
|---|---|---|
| radiated energy (modes m ≥ 1 ahead of the bunch), code units | 5.79 × 10⁻⁸ (dζ = λ_r/32), 5.86 × 10⁻⁸ (λ_r/64) | 6.00 × 10⁻⁸ (1.77 × 10⁻¹⁴ J) |
| peak of the angle-integrated spectrum, k/k_0 | 0.9922 | 0.9920 |
| FWHM of the spectrum, k/k_0 | 0.9775 … 1.0029 | 0.9780 … 1.0030 |
| on-axis wavelength (FFT of E_x) | 0.7995 µm | 0.8 µm |
| radiation pulse | 50 cycles over 5 … 45 µm, the slippage | N_u λ_r = 40 µm |

- **Convergence:** the energy converges with dζ at second order, to −2% of the integral. Halving dz, or refining the
  radial grid threefold, changes it by less than 0.2%.
- **Energy balance:** the additional kinetic-energy loss of the bunch in the undulator (K = 0.5 minus K = 0) is
  5.68 × 10⁻⁸, against 5.79 × 10⁻⁸ radiated.
- **Spectral density:** `h5_extract spectrum` agrees with the integral in absolute value to 2–5% across the line.
- **Backends:** serial, OpenMP, the emulated GPU and 2 MPI ranks give the same results.
- **Helical undulator:** the on-axis polarization is circular (equal |E_x|, |E_y|). The energy is about twice the
  planar value, and λ_r = λ_u(1 + K²)/(2γ²).
- **ctest `undulator`:** 10 periods, 15 s. It checks the energy (±5%), the spectral peak (±0.3%) and the slippage of
  the bunch (8.262 µm vs 8.260 µm from the trajectory).

**Test: high-gain FEL (SASE)** (`examples/fel_sase.deck`, 7.7 min on 2 cores). The setup:
- helical undulator, period 1 cm, K = 1, 60 periods;
- bunch at γ = 1000, 4 kA, σ_r = 10 µm, zero emittance, flat-top over 1.3 µm = 130 λ_r;
- λ_r = 10 nm, resolved by 16 slices;
- radial boundary at 1.5 mm, just beyond γL = 1.3 mm.

The FEL parameter is ρ = 0.009. The 1D power gain length is 5.1 periods, and Ming Xie's formula gives 6.4
(diffraction η_d = 0.40). 2πL_c = 55 λ_r. The amplification starts from the macroparticle noise. A copy with
1000 times less charge has no gain and serves as the reference.

| z (m) | 0.20 | 0.28 | 0.36 | 0.44 | 0.52 | 0.60 |
|---|---|---|---|---|---|---|
| radiated energy / Q², relative to the low-charge copy | 0.95 | 0.93 | 1.00 | 1.22 | 1.76 | 3.77 |

- **Gain length:** the energy above the spontaneous emission (ratio − 1) grows by 3.5–3.6 per 8 periods from
  z = 0.44 m on: a power gain length of 6.3 periods, as Ming Xie's 6.4 (with corrections 7 and 8; before them,
  the ratio reached 3.2 instead of 3.77).
- **After that:** the growth slows down. The bunch is only 2.4 × 2πL_c long, and the energy spread from its
  longitudinal space charge reaches 0.33 ρ.
- **Beam size:** the beam grows from 10 to 12.2 µm by space charge, as the envelope equation predicts.
- **Two corrections were needed** for this test (see "Corrections", items 5 and 6): bunch self-fields for v < c,
  and the Vay pusher. Before them, the beam blew up from 10 to 30 µm and there was no gain.
- **Short bunches:** a bunch of 10 λ_r (≪ 2πL_c) shows no exponential gain. The radiation slips out after
  about 2 gain lengths. The bunching peaks at 0.047, and the radiated energy is at most 3 times the
  coherent emission.

**Limitations:**
- **Helical field:** it is the on-axis field, without natural focusing. External focusing comes from the
  transport elements (quadrupoles, solenoids, plasma lenses), which may overlap the undulator.
- **Shot noise:** in this example the noise comes from the macroparticles. `bunch.shot_noise = fawley` gives
  the true shot noise of the electron number (beamlets).

## Mobile ions

`ion` blocks make ion species of the plasma mobile (element or mass, charge state,
share of the neutralizing charge). Ions are loaded at the positions of the plasma
electrons (neutral pairs) and pushed through the window like them, with
q/m = Z/M. The rest of the positive charge stays an immobile background. The
mobile-ion density is written as the column `Zni/n0` of the axis and map files;
the console shows its on-axis maximum. See `docs/input_deck.md`.

## Field ionization

`ion { element = N; charge_state = 5; ionization = adk; density = 3.5e16 cm^-3 }`
makes a species field-ionizable.

- **Rate:** the ADK tunnel rate of every level, evaluated with the instantaneous |E| at the ion (the laser carrier is resolved).
- **Energies:** ionization energies of all charge states of H, D, He, Li, Be, B, C, N, O, F, Ne, Na, Mg, Al, Si, Ar, Cu, Kr, Ag, Xe and Au, taken from the NIST Atomic Spectra Database (ver. 5.12, retrieved 2026-09-27). Provenance and NIST uncertainty flags are in `src/atomic_data.cpp`. `ionization_potentials` can replace them.
- **Algorithm:** a Monte Carlo per macroparticle, with sequential ionization of several levels within one step. The random numbers are counter-based (`src/rng.hpp`): a hash of (seed, particle, step, slice, pass). All backends therefore draw the same numbers, independent of the thread schedule. The rate is averaged over the two slices of each step. The rate of the averaged field would underestimate strongly nonlinear ionization, by 6% for H near threshold at 25 points per wavelength.
- **Freed electrons:** created at rest relative to the ion and tagged by their parent species. Ionization-injected electrons appear as `inj_<name>` in `bunches.txt`.
- **Charge-state output:** `ionization.txt` gives the charge-state distribution of the ions leaving the window.
- **Neutral gas:** `plasma.free_electrons = false` starts from a neutral gas.

## Numerical Cherenkov radiation

**Test.** A short ultra-relativistic bunch (1 pC, γ = 1000, σ_r = 2 k0⁻¹,
σ_z = 4 k0⁻¹, i.e. 16 cells) propagates in vacuum for 1000 steps (dz = 20 k0⁻¹,
2.5 mm at 0.8 µm). Without a remedy a field pattern grows behind the bunch. Its
signature is a sign change from slice to slice (kh = π), with E_x = B_y, i.e.
pure R. The field energy at the ζ-Nyquist wavenumber reaches 9% of the total.

**Cause.** The GT2 dispersion relation of a single Fourier mode is

  tan(ωΔ/2) = Δ h k⊥² sin θ / (16 sin²(θ/2) − h² k⊥² cos²(θ/2)),  θ = k_ζ h.

ω = 0 for θ = π and every k⊥. That mode is stationary in the window, like the
bunch, so the bunch's macroparticle noise drives it resonantly. Its group
velocity dζ/ds = −h²k⊥²/8 carries it backwards, slowly. The field does not act
back on the bunch coherently, so this is radiation, not an instability. It
scales with the shot noise (∝ 1/N_macro) and with the propagation distance.
Smaller dz does not help; it is a property of the ζ-discretization.

**Remedy.** `numerics.beam_filter` (default 1) applies a (1,2,1)/4 filter along
ζ to the current of the beam-list particles before the field solve. The
filter's response cos²(θ/2) is zero at θ = π, so the resonant source is gone.
The plasma current is not filtered, so the plasma response and the wake are
unchanged. The filter is the same in the serial code, the OpenMP/GPU engine
and the MPI pipeline: the results are identical up to rounding, also with a
bunch straddling a rank boundary.

**Rear end of the window.** Testing in plasma revealed a second, independent
error. At the last slice (ζ = 0) the current sum of the scheme needs the
half-cell behind the window, which was taken as zero. That halved the source of
R there, and a spike grew at the rear boundary. Now the missing half-cell
equals the last one (zero-gradient closure).

| case (W_nyquist / W_total at the end) | filter 0 | filter 1 |
|---|---|---|
| vacuum, 16-cell bunch, 1000 steps | 8.6 × 10⁻² | 1.3 × 10⁻⁴ |
| plasma wake (0.002 n_c, 0.2 pC driver, 600 slices), 200 → 400 steps | 1.3 → 1.9 × 10⁻², growing | 9 → 8 × 10⁻⁵, constant |
| `tests/cherenkov.sh`: vacuum, 300 steps | 4.3 × 10⁻² | 1.3 × 10⁻⁴ |
| `tests/cherenkov.sh`: plasma, 100 steps | 8.6 × 10⁻³ | 1.3 × 10⁻⁴ |
| `examples/lwfa_witness.deck`, 30 steps | 2.36 × 10⁻⁴ | 2.36 × 10⁻⁴ |

With the filter the ratio sits at the floor of the diagnostic: the estimator
(2F_n − F_{n−1} − F_{n+1})/4 also picks up the smooth field, by (kh)⁴/16. The
Nyquist energy itself drops by a factor of 750 in the vacuum case. In the plasma case the driver and witness parameters agree to 4 digits with and
without the filter, and the wake differs by 0.35%. In the LWFA example the
Nyquist part comes from the laser and the plasma, not from the bunch; the
witness energy changes by 2 × 10⁻⁴ MeV. Use `h5_extract <file> energy` to
monitor W_nyquist in your own runs.

**Limits.** A bunch with structure on the scale of 2–3 cells loses part of it:
one pass damps kh = 0.2 (31 cells per wavelength) by 1%, and kh = π/2 by 50%.
Such a bunch is under-resolved anyway. `beam_filter = 0` reproduces the
earlier versions except in the rearmost slice.

## Corrections of earlier versions

1. **Bunch self-fields** (found with the ion-collapse benchmark). Bunches now start with their relativistic vacuum
   field: R = 2E_perp from the discrete Poisson problem -Lap_perp psi = j_z, L = 0.
   Before, bunches started with R = 0. For v -> c the source of R
   (d_r E_z - j_r) vanishes, so the correct field never built up, and a bunch
   carried the field of a backward-moving charge (E_r = -B_phi): wrong sign and
   a spurious E_z. The Poisson problem uses the operator of the field solver, so
   a bunch in vacuum is an exact discrete equilibrium (E_z = B_z = L = 0).
2. **Deposition on the stretched grid.** Plain area weighting deposits volumes
   V_i that differ from the measures A_i of the solver's discrete Gauss law by
   O(R''/r) on a stretched grid (1.4% error in the field of a bunch in the
   transition zone). The deposition now adds a quadratic correction per interval,
   w = x + kappa x(1-x), with x linear in r^2. kappa is chosen so that
   V_i = A_i exactly: a uniform density is still deposited exactly, and the
   solver sees exactly the deposited charge.

3. **Slow particles that reach the next s-level** (found with the ionization-injection
   example, a0 = 2.5). A plasma electron in the laser moves forward by about
   (a0²/4) L sqrt(pi/2) while it crosses the pulse (γ − u_z = 1): 49 k0^-1 for
   a0 = 2.5, more than a step of 40 k0^-1. Such electrons reached z = s^1, were
   put on the beam list and then deleted there for v_z < 0.5. This broke charge
   conservation and made the fields blow up (E_z = 3 after 7 steps). Now a
   particle that reaches s^1 joins the beam list only if v_z >= vz_beam (0.9).
   Slower ones continue streaming in the next step, entering the zeta-march at
   their zeta; the first half-cell is pushed over the partial path. Beam particles
   slower than vz_min (0.8) return to streaming the same way, and no particle is
   deleted any more. Runs without such particles are unchanged bit for bit
   (regression test).

4. **Rear end of the window** (found with the numerical Cherenkov test in plasma).
   The current sum of the last slice used zero for the half-cell behind the
   window, which halved the source of R there. Now a zero-gradient closure is used.
   See "Numerical Cherenkov radiation" above.

5. **Self-fields of bunches with v < c** (found with the undulator tests at
   γ = 100). Bunches used to start with the v = c field, R = 2E_⊥ from
   −Lap_⊥ψ = j_z. A bunch with v < c slips backwards in the window, and the
   part of that field that is not its co-moving field stayed behind and
   exchanged energy with the bunch. A 16 µm, 107 pC bunch at γ = 100 gained
   0.13 % in energy in 8 cm, while its field energy fell by 19 %; this was
   independent of dζ, dz and dr. Now every bunch starts with the field of its
   own β: −(Lap_⊥ + γ⁻²∂_ζ²)ψ = j_z/β, R = (1 + β)E_⊥, L = (1 − β)E_⊥,
   E_z = −γ⁻²∂_ζψ. This is solved with the radial eigenmodes of the
   discrete Laplacian and a tridiagonal solve along ζ (`src/eigen.*`).
   - **Box size:** the box must hold the field of such a bunch, and the
     summary prints the smallest box that does (`docs/input_deck.md`, bunch
     block).
     - Radial: r_max ≥ γL; a closer wall shields the space charge like a beam
       pipe, by 1.4 % at r_max = γL for a γ = 100 bunch.
     - Window: gaps of min(10 √(a·max(a, γL))/γ, ln(100) r_max/(2.405 γ)) at
       both ends; behind the bunch after its slippage. The scheme sets E_z = 0
       at the front of the window, and a field cut off there acts as a uniform
       E_z on the whole bunch.
     - A larger box changes nothing: the FEL example gives the same bunching
       and radiation (3 digits) with r_max = 1.5 mm as with 2 cm, and runs
       16 % faster.
     - A capacitor estimate of the wall field, E_z ≈ Q/(ε0 π r_max²), does not
       apply: the field lines end radially on the wall. At r_max = 1 mm the
       capacitor estimate predicts ΔE/E = 3 × 10⁻³ for the bunch below; the
       code gives 1.6 × 10⁻⁵.
   - **Test (same bunch, no undulator, 8 cm):** mean energy 50.5889 →
     50.5899 MeV (2 × 10⁻⁵); field energy within 2 %. The gap was 70 µm, and
     the radial boundary was at 10 cm with a free-hand grid of 5 µm → 5 mm (the same with 2 mm ≥ γL).
     With a 2 µm gap the bunch loses 1.2 × 10⁻³ instead.
   - **Bunches with γ > 10⁵:** they keep the v = c field, unchanged.

6. **Pusher of the beam list** (found with the FEL tests). The Boris push does not keep the cancellation
   E + v × B ≈ 0 of a beam in its own field. Its error is a spurious force ~ (qE Δt)², and the beam sub-steps
   are large (`beam_dz_max`). At γ = 1000 (4 kA, σ_r = 10 µm, sub-step 0.625 mm) the radial force on the
   particles was 17 times the field force E_r/γ². The beam blew up from 10 to 30 µm in 0.6 m instead of
   growing by 2 µm. The beam list is now pushed with the Vay pusher (`numerics.beam_pusher = vay`), which is
   exact for E + v × B = 0 at any step: the force agrees with E_r/γ² to 1% at that sub-step. Streaming
   particles keep the Boris push; their steps are small.
7. **Undulator kick per step** (found with the GENESIS benchmark, `benchmarks/genesis_seeded`). The push kicks
   the momentum with the field at one point of the step. For the undulator field, a sinusoid in z, this
   overestimates the kick by 1/sinc(k_u Δz/2): the orbit amplitude was K/γ·(1 + (k_u Δz)²/24), 0.6 % too large
   at Δz = λ_u/16, and the resonance was shifted by 0.6 %·K²/(1 + K²). The undulator field is now multiplied by
   sinc(k_u Δz/2), the exact mean over the step.
8. **Time of the beam current in an undulator** (same benchmark). The beam push is a leapfrog: after a sub-step
   the position is at its end, the momentum half a sub-step earlier. In an undulator u⊥ turns by k_u Δz per
   sub-step, so the deposited current lagged by k_u Δz/2 in phase, a phase error of the FEL coupling of first
   order in Δz. The current is now deposited with the momentum at the time of the position (the half step
   completed with the external field, `velocityAtPosition` in `src/external.hpp`).
   In the 10 nm benchmark the power gain length went from 9.5 cm (before 7 and 8) to 8.7 cm (7) and 8.3 cm (7 and
   8); halving Δz or Δζ then changes it by less than 1 %.

## What is new compared with the prototype

- **Input deck:** the deck, with physical units, validation and overrides, replaces the `key=value` command-line parameters.
- **Several laser pulses:** each pulse has its own wavelength, focal position (initialized with wavefront curvature and Gouy phase), polarization (linear at any angle, or circular) and CEP.
- **Plasma profiles:** a piecewise-linear longitudinal profile, and a uniform, parabolic-channel or tabulated radial profile.
- **Several bunches:** electrons or positrons (the charge-to-mass ratio is stored per particle), Gaussian or flat-top shape, offsets, and a waist position.
- **Ions:** mobile ion species and field ionization (see above).
- **Free-hand grids** along r, ζ and z, and **undulators** acting on the particles (see above).
- **Beams from files:** a bunch can be read from an openPMD file of any code or from a native data file (selection, placement, resampling), and **Fawley shot noise** gives the bunching noise of the real electrons for FEL start-up.
- **Dopant profiles:** an ion block can have its own density profile (`profile = own`), e.g. a nitrogen-doped section for localized ionization injection.
- **FEL diagnostics:** `fel.txt` gives the radiation energy, peak power and bunching factor at every step, so gain curves and saturation are read directly.
- **Transport elements:** quadrupoles, active plasma lenses, solenoids and thin chicanes (R56, delay) for matching and compressing a beam between stages.
- **Function parser:** profiles of plasmas, lasers, bunches and undulator tapers as expressions in the deck.
- **GPU and multi-core:** the parallel engine runs the time step on an NVIDIA GPU (`cuda`) or on CPU threads (`openmp`); see above.
- **Diagnostics:** `bunches.txt` holds the statistics of every bunch and of the self-trapped charge at every step. The output is in µm and pC. Each run saves `deck_used.txt` and `summary.txt`.

## Benchmarks against other codes

`benchmarks/` compares GEM-PIC-RZ with established codes on the same cases (details, scripts and data in
[`benchmarks/README.md`](benchmarks/README.md)):

- **Seeded FEL vs GENESIS 1.3 v4** at 10 nm and 0.8 nm. The power gain length agrees to 0.6 % at 0.8 nm and to
  4 % at 10 nm with GENESIS runs that include the physics GEM-PIC-RZ has beyond GENESIS: longitudinal space
  charge, and the widening of the beam by its transverse space charge. The benchmark found two numerical errors
  of the undulator coupling (corrections 7 and 8 below).
- **LWFA vs FBPIC** (a0 = 2, n0 = 3.5·10¹⁸ cm⁻³): the on-axis wake of the first period agrees to 1.6 % rms after
  100 µm, at 1/50 of the computing time.

## Checks

- **Prototype equivalence:** before the corrections, the deck version reproduced the prototype bit by bit (built without `-march=native`). The regression reference was regenerated after each correction and after the beam filter was introduced.
- **Bunch field in vacuum:** an ultra-relativistic flat-top bunch has E_r equal to Gauss's law to 0.05% outside the core, on stretched, free-hand (kinks inside the bunch) and uniform grids. The rms E_z (from macroparticle noise) is below 0.2% of E_r.
- **Ion collapse:** the benchmark uses n_b = 100 n0 and hydrogen. In the compression phase, the on-axis ion density agrees with a test-ion model (analytic bunch field, full blowout) to 1-3%. Deuterium collapses later by a factor of 1.418 (√2 = 1.414).
- **Immobile limit:** an ion species with a mass of 1e15 u gives fields identical to a run without ions.
- **ADK rate:** it equals the Landau rate for hydrogen and hydrogen-like ions to 1e-12, and the code-unit tables equal the SI formula.
- **Field ionization, hydrogen:** one level, a0 = 0.01 laser, fixed atoms. The ionized fraction over the disk agrees with 1 − exp(−∫w dt) in the same field: 0.02294 against 0.02288 (seed 12345; the Monte Carlo scatter is about ±1%).
- **Field ionization, nitrogen:** N⁵⁺ → N⁶⁺ → N⁷⁺ at a0 = 2.2, mobile ions. The charge-state fractions agree with the rate equations: 0.93568 / 0.05995 / 0.00437 (scatter about ± 0.001 / 0.001 / 0.0003) against 0.93667 / 0.05927 / 0.00406.
- **Laser focusing in vacuum:** a pulse started 0.2 mm before focus (w0 = 5 µm, zR = 98 µm) follows a0/√(1+((z−z_f)/zR)²) to within 3% over 0.4 mm. At focus it reaches 0.098 instead of 0.100.
- **Bunch initialization:** the initial rms sizes, emittances, energy spread and flat-top length agree with the deck.
- **Function profiles:** a deck with the plasma ramp, the Gaussian laser and a flat-top bunch written as functions (`tests/functions.sh`) reproduces the built-in profiles to rounding. In a tapered undulator, K = 0.5 → 0.8, the slippage is 8.83 µm against 8.78 µm from ∫(1 + K²/2)/(2γ²)dz.
- **Sign checks:** in the PWFA example the electron witness is accelerated and focused. The same bunch as positrons is decelerated and defocused.

## Known limitations

- **Performance:** the ζ-march is sequential within a step. Parallelism comes from the particles and the radial grid within a slice (OpenMP/GPU), or from the pipeline over steps (MPI, serial step per rank). MPI combines with OpenMP threads (`backend = openmp`), not yet with the GPU engine.
- **Missing physics:** the energy that ionization takes from the laser (ionization current) and collisional ionization are not modelled. Ion species follow the plasma profile or their own (`profile = own`).
- **Laser envelope:** the laser is initialized as a paraxial Gaussian beam (or, with `profile_function`, as the given field at z = 0 with flat phase fronts). The temporal envelope is not corrected for focusing (spatio-temporal coupling), which is fine for pulses much longer than a cycle.
- **Radial absorber:** it is a simple damping of R and has not been tested for reflection.
