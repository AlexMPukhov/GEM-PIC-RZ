# GEM-PIC-RZ

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23600178.svg)](https://doi.org/10.5281/zenodo.23600178)

Galilean electromagnetic particle-in-cell code in quasi-cylindrical geometry
(azimuthal modes), with free-hand non-uniform grids along r, ζ and z: the
step is given at reference points and interpolated linearly in between. The
scheme is described in A. Pukhov, [arXiv:2610.31622](https://arxiv.org/abs/2610.31622)
(see "Citing GEM-PIC-RZ" below). This project grew out of
a single-file prototype (not included). Two discretization errors of the
prototype have been corrected (see "Corrections"), so results with particle
bunches differ from the prototype.

## Citing GEM-PIC-RZ

If you use GEM-PIC-RZ, please cite the paper

> A. Pukhov, *GEM-PIC-RZ: a Galilean electromagnetic particle-in-cell code in quasi-cylindrical geometry with
> nonuniform grids*, [arXiv:2610.31622](https://arxiv.org/abs/2610.31622) (2026),

and the method paper: A. Pukhov, N. Elkina, and T. Wilson, *Galilean electromagnetic particle-in-cell code*,
[J. Comput. Phys. **567**, 115331 (2026)](https://doi.org/10.1016/j.jcp.2026.115331).
For the code itself: A. Pukhov, *GEM-PIC-RZ*, Zenodo, [doi:10.5281/zenodo.23600178](https://doi.org/10.5281/zenodo.23600178) (all versions;
version 1.0.0: [doi:10.5281/zenodo.23600179](https://doi.org/10.5281/zenodo.23600179)).
`CITATION.cff` has the machine-readable form (GitHub: "Cite this repository").

The code is released under the BSD 3-Clause License (`LICENSE`). Version: `gempic_rz --version`.

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
python3 tools/gempic_gui.py my.deck                             # GUI in the web browser (see below)
```

### GUI

`tools/gempic_gui.py` is an optional graphical front end for desktop and notebook runs, like the one of
QUARZ. It is a small local web server (Python standard library, numpy, and h5py for the HDF5 frames) and
a page without external dependencies, so it works offline:

```
python3 tools/gempic_gui.py examples/lwfa_witness.deck            # opens the browser
python3 tools/gempic_gui.py my.deck --no-browser --port 8766       # remote: ssh -L 8766:127.0.0.1:8766 host
```

- **Deck:** editor with the example decks, Save, Check (`gempic_rz deck --check`), extra command-line
  overrides that are not saved in the deck.
- **Run:** threads (`numerics.backend = openmp` for more than one), MPI ranks, Stop, the log and the progress.
  The deck is checked first; a deck with errors is refused with the message.
- **Live frames:** the GUI adds `hdf5` to `output.format` (and sets `output.every` if "frames every" is
  given) on the command line; the deck file is not changed. Every output step then writes
  `data_NNNNN.h5` (written as `.tmp` and renamed, so a frame is never read half-written).
- **Plots:** 2D map in the plane y = 0 (x = +r at φ = 0, x = −r at φ = π) of E_x, E_y, E_z, B_x, B_y, B_z
  (code units m c ω₀/e, m ω₀/e) and of n_e/n₀, Z n_i/n₀, with a robust colour range; click to set the
  cuts along ζ (at x) and across (at ζ); frame slider and playback; the bunch table `bunches.txt`
  (energy, spread, charge, sizes, emittances against z). Finished runs can be browsed the same way.
- **Fields:** built from the modes as the particle gather does (staggered faces and centres, regular
  values on the axis). The test `gui` compares them with the code's own axis and map output to the
  print precision of those files, runs a 2-rank MPI run through the GUI, and checks the token and the
  error handling. Without h5py, or with a build without HDF5, the GUI shows the text maps
  `map_NNNNN.txt` (E_x, E_z, densities).
- **Security:** the server listens on 127.0.0.1 only, and every request needs the random token in the
  URL it prints.

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
  `profile_z(zeta)` and energy chirp `gamma_function(zeta)`, its transverse profile `profile_r(x, y, r, zeta)` and mean
  momenta `ux_function`, `uy_function` (rings, helical orbits), and an undulator taper `K_function(z)`. Expressions accept
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
| `src/laser_file.*` | lasers from LASY / openPMD envelope files: reader, azimuthal modes, exact propagation (Hankel + FFT), envelope → R |
| `src/h5read.hpp` | small read-only HDF5 helper shared by the importers |
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
| `src/solver.*` | tridiagonal field update per slice and mode (matrices per distinct half-cell length and step); pentadiagonal with the per-mode implicitness `numerics.theta_stiff` |
| `src/smoothing.*` | radial smoothing of the plasma sources (`plasma.smooth_length`): vector and scalar filters (1 − a²∇²)⁻¹ that commute with the discrete divergence |
| `src/particles.hpp` | particle record, Boris pusher |
| `src/simulation.*` | laser and bunch initialization, plasma loading, pushers, deposition, time step (serial reference) |
| `src/engine.*` | parallel engine (GPU / OpenMP): the time step as kernels, device-resident state |
| `src/pipeline.*` | MPI: domain decomposition along ζ, pipelined in the steps; diagnostics and output on rank 0 |
| `src/comm.*` | thin MPI layer (packed messages; a single-rank stub without MPI) |
| `src/portable.hpp`, `src/hd.hpp` | portability layer: memory, kernel launchers, teams, atomics, counter-based RNG, complex numbers |
| `src/diagnostics.cpp` | console output, `bunches.txt`, axis/map/particle files |
| `src/output_hdf5.cpp` | HDF5 output (`data_NNNNN.h5`), openPMD output, LASY envelope files, checkpoints (`checkpoint_NNNNN.h5`) and restart |
| `src/main.cpp` | command-line driver |
| `tests/test_deck.cpp` | about 230 checks (the number depends on the build options): parser, units, configuration, backends, output format, free-hand mapping and ζ/z grids, ion and ionization blocks, error messages, deposition shape, NIST table, ADK = Landau rate; all example decks must be valid |
| `tests/regression.*` | runs `tests/regression.deck` (laser, plasma, mobile H ions, witness) and compares with `tests/reference/`, for every backend |
| `tests/physics.*`, `tests/check_physics.cpp` | bunch self-field in vacuum against Gauss's law (stretched and free-hand radial grid, free-hand ζ grid); a γ = 100 bunch keeps its energy in its own field; ion collapse against a test-ion model; field ionization of H and N⁵⁺→N⁷⁺ against the ADK rate equations |
| `tools/h5_extract.cpp` | reads the HDF5 files and the native records of the openPMD files: tables of axis, map, beam, field line-outs (also an example reader) |
| `tools/gempic_gui.py`, `tools/gempic_gui.html` | GUI: deck editor, run control and live plots in the browser (see "GUI"); test `tests/gui.sh`, `tests/gui_check.py` |
| `tools/openpmd_rectilinear.py` | Python reader of openPMD meshes on regular and rectilinear (non-uniform) grids, GEM-PIC-RZ and QUARZ |
| `tests/hdf5.sh` | HDF5 files read back against the text output (double and single precision) |
| `tests/flexgrid.sh`, `tests/flexgrid.deck` | free-hand ζ and z grids: constant-step tables = uniform grid; parallel backend = serial; restart, also with a new z grid |
| `tests/undulator.sh`, `tests/undulator.deck` | undulator radiation: energy, spectrum and slippage against the far-field radiation integral (`tools/undulator_radint.py`) |
| `tests/smoothing.sh` | `plasma.smooth_length`: a → 0 gives the reference, the axis density spike is capped, the fields ahead of the bubble back hardly change, parallel backend = serial |
| `tests/lasy.sh`, `tests/lasy_files.py` | LASY files: import against the built-in Gaussian (thetaMode, cartesian, circular), exact propagation, export and re-import (read by LASY if installed), boundary injection against the built-in flying focus |
| `tests/energy_spread.sh`, `tests/energy_spread.py` | dE/E of `bunches.txt` against the particle data (two-pass, extended precision): γ = 1000, dE/E from 0 to 2 × 10⁻⁵, the printed value must be the rounded reference |
| `tests/cherenkov.sh` | numerical Cherenkov: ζ-Nyquist field energy of a relativistic bunch in vacuum and in plasma, with and without the beam filter |
| `tests/thin_ring.sh`, `tests/thin_ring.py`, `tests/thin_ring.deck` | thin, intense, low-γ ring (`beam_consistent`, `theta_stiff`): static and prescribed-expansion rings at γ = 10⁶ against the exact field relations; a 10 kA, γ = 50 ring at rest and on matched orbits in a plasma lens against a long-beam shell model, per backend |
| `tests/mpi.sh` | MPI: one-rank pipeline = serial; 2–4 ranks against serial (also `beam_field = mean`); restarts across rank counts; wall-time stop |
| `tests/restart.sh` | interrupted + restarted run against the uninterrupted run (checkpoint file, wall-time stop) |
| `tools/install_cuda_pip.sh` | installs nvcc + CUDA runtime from PyPI into a user directory (no root) |
| `tools/compare_columns.cpp` | numeric comparison of ASCII tables (used by the regression test) |
| `tools/plot_maps.gp` | gnuplot → PDF: E_z and density maps, on-axis lineout |

## Lasers from and to LASY files

GEM-PIC-RZ reads and writes laser pulses as [LASY](https://github.com/LASY-org/lasy) envelope files (openPMD,
mesh `laserEnvelope`), so a pulse can be prepared in LASY (measured data, axiparabola or flying-focus optics,
chirps, aberrations), passed between codes, and handed from one stage to the next. Details and all keys:
[`docs/input_deck.md`](docs/input_deck.md#lasers-from-lasy--openpmd-envelope-files-file).

```
laser  { file = axiparabola_00000.h5; zeta_center = 60 um; file_z = 2 mm; inject = boundary }
output { lasy = true }                    # diag/lasy/laser_NNNNN.h5 at z_end, readable by LASY
```

- **No conversion needed:** in GT2 a level s is the plane z = s and its slices are the times t = s − ζ at which the
  fields pass it. That is exactly LASY's representation (the envelope at a plane versus time), so a file is a level
  of the code, with all its spatio-temporal structure and its absolute carrier phase.
- **Import:** at the file's plane the envelope is interpolated; at another plane (`file_z`) it is propagated
  exactly, per azimuthal mode (Hankel transform, exact for the interpolant) and frequency (FFT), in vacuum or a
  uniform plasma. Thetamode (any number of modes) and cartesian files.
- **Boundary injection** (`inject = boundary`): the field is also fed through the radial absorber layer at every
  step (total field / scattered field, as the `flying_focus` block), so converging waves of real optics come in from
  the side over the whole length. This is open task 6 of the flying-focus handover.
- **Export** (`output.lasy`): the forward field at the plane z = s, demodulated at ω_L and projected on its
  dominant polarization, in LASY's layout.
- **Checks** (ctest `lasy`): built-in Gaussian vs the same pulse from a file (1.2·10⁻⁴ of the peak), exact
  propagation vs the paraxial Gaussian (0.1 % in the peak), export → re-import (7·10⁻⁴), and the luminal Bessel
  pulse of the `flying_focus` block given as a file with `inject = boundary` vs the built-in flying focus after
  0.3 mm (peak 0.25 %, rms 1.3·10⁻³). LASY reads the exported files. All backends give the same fields.

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
- **Mean-field beam push** (`numerics.beam_field = mean`, the default with undulators): the second push of the beam
  in (F⁰ + F¹)/2 needs F¹ below the block, which the rank behind computes later. Each rank therefore does the second
  push of step k at the start of step k+1, after the exchanges above have brought F¹ below its block (rear → front)
  and the second-push results of the rank in front (front → rear, with the next end-of-step message). It runs from
  the front to the rear like the first push, with the same halo copies and continuation on the rank behind; a tag
  on every beam particle pairs its two pushes across ranks, so that the rule of the serial code holds (a particle
  the second push would lose keeps its first-push state). The output of step k follows its second push; the last
  step takes one extra exchange.

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
- **Mean-field beam push:** the pipelined step on one rank is identical to the serial step; 2, 3, 4 ranks and a
  checkpoint on 3 ranks continued on 2 agree with the serial run to all printed digits (regression deck with
  `beam_field = mean` and beam filter 2: wake, trapping, absorber). The 4 kA beam in a planar undulator at λ_u/16
  (`tests/planar_step.sh` settings, 5 cm, the beam spans the 3 blocks): ε_nx = 0.0067 µm on 3 ranks as serially
  (`level`: 0.118 µm); `fel.txt` and `bunches.txt` agree to all printed digits.
- **Hard case:** the ionization-injection example at a0 = 2.5 on 3 ranks (13 steps: trapping, re-streaming of 1936 particles, field ionization). `bunches.txt` and `ionization.txt` are identical to the serial run, and the axis fields agree to rounding.

**Current limits:**

- **Backend:** `numerics.backend = serial` runs the serial step on every rank; `numerics.backend = openmp` adds OpenMP threads on every rank (hybrid MPI + threads, `OMP_NUM_THREADS` threads per rank). `numerics.backend = cuda` runs the ζ-march of every block on a GPU (one GPU per rank, chosen by the node-local rank of the launcher, or the one GPU a task sees with `srun --gpus-per-task=1`): the beam push (`Engine::pushBeams`, also of the halo copies), the beam deposits and the filter across the block edges, then the streaming plasma, the field solves and the R corrections of the own slices as kernels (`Engine::marchBlock`), staged through the host for the messages. The mean-field second push (`numerics.beam_field = mean`) runs on the device as well: the pushes in (F⁰ + F¹)/2 (F¹ of the absorber layer before the absorber, as serially), the deposits and the filter. The host keeps the beam list and its bookkeeping (which particle stays, goes to the rank behind or returns to streaming, which push of a particle counts). Tested in the GPU emulation build (`GEMPIC_PIPELINE_ENGINE=1 gempic_rz_gpuemu ... numerics.backend=openmp`: the GPU code paths on CPU threads): 1–4 ranks, a restart from 3 to 2 ranks, the mean-field push, free-hand grids, field ionization, the flying focus and the FEL beam across 3 blocks agree with the serial run to all printed digits (`tests/mpi.sh`). On a real GPU (RTX 5090) one rank with the pipelined step (`GEMPIC_PIPELINE=1`) equals the serial run to all printed digits on the regression, FEL and flying-focus decks (`benchmarks/gpu/check_1gpu.sh`); several GPUs (MPI) have not run yet.
- **Memory:** every rank builds and stores the fields and currents only of its own block of slices and ghost slices
  (2 × beam-filter passes + 6 on each side); the particles are distributed as well. Nothing is built for the whole
  window on any rank at the start: lasers (also from LASY files), the flying focus and the bunch self-fields are
  computed per block (for bunches with v < c the solve along ζ runs across the ranks, with the arithmetic of the
  serial code), and a restart reads only the block of each rank from the checkpoint. The memory per rank therefore
  falls as 1/P, which makes long windows (XFEL, multi-period wakes) possible. The results are identical to those with
  the whole window on every rank (`GEMPIC_FULL_STORAGE=1`, testing). Example: `tests/regression.deck` with a 300 µm
  window (9400 slices) on 4 ranks: 0.43 GB per rank during the steps; the whole window on every rank needed more
  than 1 GB per rank and did not fit on a 4 GB machine.
- **Output without the whole window:** every output file (text, HDF5, openPMD, LASY, checkpoints) is written block
  by block: the datasets are created at their full size and filled with hyperslabs. Rank 0 receives the blocks of
  the ranks one after the other and never holds more than two of them (the next block supplies the one slice and
  half-cell that the interpolations need beyond a block). The files are the same as those of a serial run. In the
  example above, with HDF5, openPMD, text and checkpoint output, the largest rank needs 0.60 GB. Exception: the
  LASY export keeps R at its output radii for all time samples (its demodulation needs the whole time series of a
  radius); the beam list (bunches, trapped particles) is gathered on rank 0 for the files.

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
  - field solve, with one thread block per azimuthal mode: the radial operators are applied in parallel, and the tridiagonal solves use the Thomas/PCR partition method on all threads of the block (`numerics.tridiag`; Thomas in one thread on CPU threads);
  - slice preparation: the predictor and the |E| bound for ionization.

There are no host-device transfers inside the march. After each step, the fields,
densities and beam list are copied back for the diagnostics.

The kernels are written once, as lambdas, behind a thin portability layer (`src/portable.hpp`). The same source compiles with nvcc for the GPU and with a C++ compiler for OpenMP threads.

**Validation:**

- **`openmp`:** reproduces `serial` bit for bit, provided both are compiled without FMA contraction (`-DGEMPIC_NATIVE=OFF`). This was checked on the ionization-injection example at a0 = 2.5 over 20 steps: field ionization, self-trapping, and 10 750 particles returned from the beam list to streaming. With `-march=native` the compiler contracts the complex arithmetic of the two codes differently, so last-bit differences appear. They stay at rounding level on the regression deck (1200 axis values equal to 7 digits, on 1 and 2 threads) but grow in chaotic runs with trapping.
- **Emulated GPU paths:** the CPU build also makes `gempic_rz_gpuemu` (ctest `regression_gpu_emulation`). It runs the GPU code paths on CPU threads: teams of 8 threads with real barriers, atomic deposition, and unordered appends to the particle lists. It agrees with `serial` to 1e-9 on the regression deck and passes all physics tests.
- **Ionization:** the `openmp` and emulated runs give identical results, so the counter-based random numbers do not depend on the thread schedule.
- **Real GPU:** run `ctest` in the GPU build. `regression_cuda` compares the GPU run with the serial reference (1e-6 relative); `physics_cuda` runs the bunch-field, ion-collapse and ionization benchmarks on the GPU. First run on a real GPU (RTX 5090, 2026-10-25): all GPU tests pass, GPU = serial to all printed digits; numbers and timings in `benchmarks/gpu/README.md`.

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

openPMD 1.1 has no non-uniform grids. Besides the interpolated fields the files therefore hold the exact data
on the code's own nodes (`output.openpmd_grid = both`, the default; or `native`): records `native_*` whose
node positions are data sets referenced by the mesh attribute `gridCoordinates`. This is the form we propose for
openPMD 2.0 ([`docs/openpmd_rectilinear.md`](docs/openpmd_rectilinear.md)); QUARZ writes an interim variant.
`h5_extract` and `tools/openpmd_rectilinear.py` read both.

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
warns otherwise. With undulators the beam is pushed in the mean of the fields of the two levels of a step
(`numerics.beam_field = mean`, "Corrections", item 9); this keeps intense beams in planar undulators stable at λ_u/16.

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
  the true shot noise of the electron number (beamlets); `shot_noise = quiet` gives the same beamlets without
  noise (quiet start, for seeded amplifiers compared with codes that load a quiet beam).

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

## Bubble closure: regularization (`plasma.smooth_length`)

**Problem.** Behind the first bubble of a laser wake the on-axis E_z can show oscillations on the scale of the
laser wavelength that FBPIC does not have (`benchmarks/lwfa_fbpic`, ζ < 13 µm at 100 µm). They come from the back
of the bubble, where the electrons that close it converge on the axis: for a cold plasma and a perfectly
axisymmetric driver this is a caustic, with a density that grows without bound as the cells are refined (124 n0
on the axis in `tests/regression.deck`). Two numerical effects of the scheme act there.
- **Fast spurious modes.** The reduced field equations have a resonance at k_ζ = ω_p/2 where ω → ∞ (the
  dispersion relation ω = 2k_ζω_p²/(4k_ζ² − ω_p²)). The trapezoid in s maps such modes onto undamped
  modes that change sign from step to step. With the local ω_p of the density spike they fall on resolved
  scales. With Δ = 5 µm instead of 2.5 µm, they grow at the bubble back: E_z reaches ±0.5 instead of 0.04, and
  the axis density alternates between 30 and 1800 n0 from step to step.
- **Switching levels.** Sheath and laser-driven electrons move forward by more than Δ while they cross the
  window. When their z reaches the next level they are taken out of the current level and continue in the next
  one from the same ζ; every such switch is a small jump in the current. Their number grows as Δ decreases,
  and so does the noise: with Δ = 1.25 µm the rms deviation from FBPIC behind the bubble is 5 times larger than
  with 2.5 µm (0.0098 vs 0.0018 m c ω0/e).

**Remedy**, as in QUARZ: the sources of the streaming plasma (electrons, ions, freed electrons) are smoothed
radially after the deposit of every half-cell,

    j⊥ ← (1 − a² ∇²_vec)⁻¹ j⊥,     (ρ, j_z) ← (1 − a² ∇²_m)⁻¹ (ρ, j_z),

per azimuthal mode. Every ring of charge gets a finite width a, which caps the spike. ∇²_vec = grad div − rot curl
is built from the operators of the field solver on the staggered grid. div ∇²_vec = ∇²_m div then holds exactly, so
the filter commutes with the divergence and the smoothed currents conserve the smoothed charge (unit test in
`test_deck`). The vector system is pentadiagonal (j_r and j_φ interleaved), the scalar one tridiagonal. Both are
factorized once, so the filter costs less than 1 % of the run time. The beam list and the fields acting on the
particles are not smoothed. MPI, OpenMP and the GPU engine give the same results as the serial code.

| LWFA vs FBPIC at 100 µm (`benchmarks/lwfa_fbpic`), rms of E_z − E_z(FBPIC) [m c ω0/e] | Δ = 2.5 µm | Δ = 1.25 µm | Δ = 5 µm |
|---|---|---|---|
| behind the bubble (2 … 12.5 µm), cold plasma | 0.0018 | 0.0098 | 0.053 (spike ±0.5) |
| behind the bubble, a = 1 k0⁻¹ | 0.0010 | 0.0017 | 0.030 |
| behind the bubble, a = 2 k0⁻¹ (0.25 µm, 0.09 c/ω_p) | **0.0007** | **0.0007** | 0.0079 (stable) |
| first wake period (14 … 38 µm), cold / a = 2 k0⁻¹ | 0.0004 / 0.0005 | 0.0005 / 0.0005 | |

- **Convergence:** with a = 2 k0⁻¹ the wake behind the bubble is the same for Δ = 2.5 and 1.25 µm, and matches
  FBPIC.
- **Cost of the model change:** the smoothing shifts the bubble back by about a. The first period changes by
  0.0001 rms (1.6 % → 1.9 % rms of the peak against FBPIC), and the peak E_z at the back changes from 0.047
  (spike) to 0.040 (FBPIC 0.042).
- **Choice of a:** at least the axis cell and the slice spacing. Off by default.
- **Dense plasmas** (0.01 n_c, `benchmarks/laser_evolution`): the fast modes also grow in the nonlinear wake several
  periods behind the laser, with E_z up to 20 times the FBPIC value. They need a ≈ 0.3 c/ω_p for the wake (mode
  m = 0), while the laser-driven currents (m = 1) need only ≈ 0.1 c/ω_p and are changed by more:
  `plasma.smooth_length_m0` sets a for m = 0 separately (here 3 and 1 k0⁻¹).

**Other switches** (both off by default, kept for studies):
- `numerics.trap_advance = k`: an electron leaves its level only after advancing k steps. With k = 4 to 8
  (5 … 10 µm) the noise from switching levels disappears at Δ = 2.5 and 1.25 µm (rms 0.0015 … 0.0017).
- `numerics.trap_gamma = g`: a particle that advances beyond its level is moved to the beam list (trapped) only if
  its γ ≥ g. Electrons in a laser of a0 = 3 quiver with γ ≈ 1 + a0²/2 and v_z ≈ 0.9 c; they advance far beyond Δ
  while crossing the pulse without being trapped, and with g = 0 (default) are counted as trapped.
- `numerics.theta > 0.5`: an off-centred update of R in s that damps the fast spurious modes by (1 − θ)/θ per
  step. θ = 0.6 gives rms 0.0003 for Δ = 1.25, 2.5 and 5 µm. It is first order in Δ for any oscillating source,
  so it must not be used with undulators.

## Thin, intense beams: consistent beam sources and per-mode implicitness

A thin, intense, low-γ beam structure (here a ring of radius 0.75 µm and rms width 20 nm, half a radial cell,
10 kA, γ = 50, on matched circular orbits in a plasma lens) needs its radial self-force, which is only 1/γ² = 4 × 10⁻⁴
of E_r, to that accuracy. Three errors of the standard scheme spoil it (numbers from `tests/thin_ring.sh` and the
private development runs):

1. **Static mismatch of j_z.** The discrete Ampère law of the E_z equation relates the mean of div B at the slices
   ζ_n, ζ_{n+1} to j_z of the half-cell between them, while the initial bunch field and every R update attach to a
   slice the mean of the two neighbouring half-cells. The difference, −(h²/4) ∂²_ζ j_z, drives a spurious E_z along
   ζ; at the edges of a flat-top it is O(j_z) for any ζ resolution, and d_r E_z ~ 1/w at a sheet of width w enters
   L = E_r − B_φ. A static ring at γ = 10⁶ (an exact steady state at s = 0) got E_z = 4.5 × 10⁻³ and L = 6 × 10⁻⁴
   of max E_r after one step.
2. **Non-conserving radial current.** j_r (faces) and j_z (centres) were deposited with separate linear weights, so a
   sheet thinner than a cell left a Gauss-law error ("ghost" charge) behind wherever it moved.
3. **Undamped odd-even mode of the stiff E_z constraint.** For a radially moving β → 1 charge, E_z jumps across the
   sheet such that d_r E_z = −j_r cancels j_r in the L equation (the field of each charge is a tilted pancake). At
   large Δ the scheme enforces θ G(E_z¹) + (1 − θ) G(E_z⁰) ≈ −j_r; with θ = 1/2 that is E_z¹ = −2 j_r − E_z⁰, an
   undamped odd-even mode in s excited by every change of j_r, and L sees E_z¹, not the mean.

`numerics.beam_consistent = true` fixes 1 and 2: the half-cell equations (E_z, L_r, L_φ) get the (¼, ½, ¼) mean of the
beam current the R update uses, and the R update gets the charge-conserving radial flux of the step (per particle
the change of its centre weights between the two positions, with the j_φ of the R_φ update for m ≥ 1; with
`beam_field = mean` the flux starts from the first-push positions of the last step, with which R agrees).
`numerics.theta_stiff = 1` fixes 3: θ becomes an operator, θ_k = (θ + μ_k)/(1 + μ_k) with μ_k = c h Δ |λ_k| per
radial eigenmode of E_z and B_z, so θ_k = θ for every mode the step resolves (all radiation: c h Δ |λ| ≪ 1, damping
negligible) and θ_k → 1 (L-stable) for the stiff near-field modes. The implicit matrices become pentadiagonal.
With `beam_consistent`, bunches also start with the co-moving field of their azimuthal current (B_z, E_φ of a
rotating beam; without it a rotating ring gets a large torque in the first step), discretized as the steady state of
the scheme (E_φ from the difference of A_φ across half-cells).

| ring test (0.5 mm, dz = 20 µm) | old scheme | `beam_consistent` + `theta_stiff = 1` | reference |
|---|---|---|---|
| static ring, γ = 10⁶: max L / max E_r after a step | 6 × 10⁻⁴ | 5 × 10⁻¹³ | 0 |
| prescribed expansion v_r = 3 × 10⁻⁴ c: max L / max E_r | 1.7 × 10⁻³ | 2 × 10⁻⁷ … 9 × 10⁻⁶ | O(v_r²) |
| 10 kA ring at rest, no lens: r 1/50/99 % (µm) | 0.69 / 2.49 / 4.48 | 0.74 / 1.37 / 1.95 | 0.73 / 1.44 / 1.98 (shell model) |
| same ring on matched orbits in a lens: u_r rms | 0.32 | 0.0025 | 0.0022 (shell model) |
| angular momentum ⟨r u_φ⟩ (start 0.5365) | 1.27 | 0.535 | 0.537 (shell model) |

The reference is a long-beam shell model (every electron a charged, current-carrying cylindrical shell; Gauss,
Ampère, B_z of the rotating shells, the lens; `tests/thin_ring.py`). With one azimuthal mode only (m_max = 0) the old
breathing instability is gone completely; m = 1 needed the charge-conserving flux with j_φ. Beyond ~0.5 mm the
longitudinal space charge of the sharp-edged flat-top (±10 % in γ per mm at head and tail, physical) mismatches head
and tail to the lens; the shell model, which has no E_z, does not describe that. Both switches are off by default;
`beam_consistent` is not yet available in the MPI pipeline.

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
9. **Beam push in the fields of one level** (found with the planar harmonics benchmark, `benchmarks/harmonics`).
   The beam list was pushed from s⁰ to s¹ in the fields of s⁰, first order in the step. In an undulator the
   large self-field of the beam turns with the wiggle angle, so the push used the field of a beam tilted by the
   angle of s⁰, and E + v × B no longer cancelled. A planar undulator with 4 kA (γ = 1000, K = 1.5, σ_r = 10 µm,
   zero emittance) at Δz = λ_u/16 blew the beam up from the tail forward: ε_n 0 → 2.8 µm and σ_x 10 → 21 µm in
   25 periods (λ_u/32: 0.04 µm), and the gain length came out 1.8 times too long. Now, with undulators, the beam
   is pushed again after the march of s¹, from s⁰ in the mean (F⁰ + F¹)/2, and its current of s¹ is deposited
   anew (`numerics.beam_field = mean`). At λ_u/16: ε_n 0.04 µm, σ_x 10.41 µm in 25 periods, the same as a full
   predictor–corrector with a second solve of s¹ (10.413 µm) and as the converged λ_u/64 run. Stable at λ_u/8.
   Extrapolating the field from s⁻¹ and s⁰ instead made it worse. Without undulators nothing changes
   (`auto` = `level`), bit for bit. Test: `tests/planar_step.sh`.

## What is new compared with the prototype

- **Input deck:** the deck, with physical units, validation and overrides, replaces the `key=value` command-line parameters.
- **Several laser pulses:** each pulse has its own wavelength, focal position (initialized with wavefront curvature and Gouy phase), polarization (linear at any angle, or circular) and CEP.
- **Plasma profiles:** a piecewise-linear longitudinal profile, and a uniform, parabolic-channel or tabulated radial profile.
- **Several bunches:** electrons or positrons (the charge-to-mass ratio is stored per particle), Gaussian or flat-top shape, offsets, and a waist position.
- **Ions:** mobile ion species and field ionization (see above).
- **Free-hand grids** along r, ζ and z, and **undulators** acting on the particles (see above).
- **Lasers from and to LASY files:** import (exact propagation, boundary injection) and export of envelope files.
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
  100 µm, at 1/50 of the computing time. Behind the bubble the cold-plasma run has spurious oscillations that
  `plasma.smooth_length` removes (see "Bubble closure").
- **Ionization injection vs FBPIC** (N5+ dopant, a0 = 2.3, 200 µm): the injected charge above 5 MeV agrees to 10 %
  and its mean energy, energy spread and spectrum to 2–5 % at the same lab position of the bunch, at 1/50 of the
  computing time.
- **Laser evolution vs FBPIC and Wake-T** (6.8 fs, a0 = 1.5, 0.01 n_c, 300 µm): laser energy (depletion) to 1.3 %,
  red-shift to 0.002 k0 and peak field to 5 % against FBPIC, at 1/60 of the computing time; the envelope code
  Wake-T underestimates the depletion fivefold for this 2.5-cycle pulse.

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
- **Laser envelope:** the built-in laser is initialized as a paraxial Gaussian beam (or, with `profile_function`, as the given field at z = 0 with flat phase fronts). The temporal envelope is not corrected for focusing (spatio-temporal coupling), which is fine for pulses much longer than a cycle. A pulse from a LASY file (`laser.file`) has no such limitation: it is propagated exactly.
- **Radial absorber:** it is a simple damping of R and has not been tested for reflection.
- **Density down-ramp injection (open, future work):** not reproduced (2–5x less trapped charge than FBPIC, step
  dependent; `benchmarks/downramp_fbpic`). Attributed to the GT2 frame of the particles; GT1 should describe it but
  produces stripes in the fields where plasma particles reach the next slow-time level inside the window. No solution
  yet; to be discussed before work starts.
