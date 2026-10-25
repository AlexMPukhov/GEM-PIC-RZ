# GEM-PIC-RZ on a GPU: what to run

Build (see `docs/gpu.md`):
```
cmake -S . -B build-gpu -DGEMPIC_CUDA=ON -DCMAKE_CUDA_ARCHITECTURES=native
cmake --build build-gpu -j
```
Then `cd benchmarks/gpu && sh gpu_check.sh`. It runs four steps and prints the numbers to report back.

| step | deck | what it checks | CPU reference (2 cores, OpenMP) |
|---|---|---|---|
| 1 | `ctest` in `build-gpu` | GPU = serial reference to 1e-6 (`regression_cuda`), bunch fields, ion collapse, ionization (`physics_cuda`), restart (`restart_cuda`) | all pass on CPU |
| 2 | `../lwfa_fbpic/gem.deck` (LWFA, a0 = 2, 3.5e18 cm^-3, 100 µm, 40 steps) | physics on the GPU against FBPIC: `compare.py gem_out_gpu` prints the same differences as for the CPU run `gem_out` (on-axis E_z of the first period 1.6 % rms). The GPU run is not bit-identical to the CPU (atomic deposition): expect differences at the 1e-6 level, more after self-trapping | 1.6 s per step |
| 3 | `lwfa_heavy.deck` (the same LWFA, 8 × 32 particles per cell, ~4·10⁴ per slice, 500 µm, 200 steps) | speed, plasma-dominated; `GEMPIC_PROFILE=1` prints the time per phase | 10.6 s per step |
| 4 | `../harmonics/fel.deck` (seeded planar FEL, 640 000 beam particles, 3 modes × 122 radial × 18 372 slices; needs ~6 GB), first 30 mm (48 steps at dz = λ_u/16) | speed, beam- and field-dominated; P₁ in slice 14 at step 48 should equal the CPU value 9.96529·10⁶ W to ~1e-4 | 1.6 s per step (20 threads) |

What to expect / look for:
- The ζ-march is sequential: each slice launches about 10 kernels. In the 100 µm deck (~1600 plasma particles
  per slice) the GPU is launch-latency bound and may hardly beat 2 CPU cores; the heavy deck shows the real
  speed-up. If step 3 is also slow, the `GEMPIC_PROFILE` output tells whether the field solve (one block per mode,
  tridiagonal solves in one thread) dominates; that would be the place to bring in a parallel tridiagonal
  solver, as was done in QUARZ.
- This check uses one GPU. Several GPUs (the MPI pipeline, one rank per GPU) and the scaling series:
  `benchmarks/jupiter/README.md`.
- Please send back: the `nvidia-smi` line, the ctest summary, the `compare.py` lines of step 2, the `done:` lines
  and the profile of steps 3 and 4.

## Without Slurm and MPI: `check_1gpu.sh`

`sh benchmarks/gpu/check_1gpu.sh [build-gpu] [build-cpu] [work dir]` (from the GEM-PIC-RZ directory; `build-cpu` is a
CPU build with OpenMP). It runs the comparisons of `benchmarks/jupiter/check.sbatch` on one GPU, for one rank,
and adds wall times against `THREADS` CPU threads (default: all cores):

1. the regression deck: the single-GPU engine (whole window) and the pipelined step on one rank
   (`GEMPIC_PIPELINE=1`, the block march of the MPI code) against serial, to 1e-6 relative;
2. the FEL deck (4 kA beam, mean-field push) and the flying focus through the absorber layer, pipelined, against serial;
3. the 4 kA beam in a planar undulator at λ_u/16, pipelined: ε_nx at 5 cm (serial 0.0067 µm);
4. wall times of `lwfa_heavy.deck` and of the FEL deck to 30 mm, GPU against the CPU threads.

The last line is `summary: 0 failed comparison(s)` when all is well.

## First real-GPU run: RTX 5090 (2026-10-25)

A workstation of the group: NVIDIA GeForce RTX 5090 (compute capability 12.0, 170 SMs, 32 GB), driver 580.95,
CUDA 13.0, gcc 14.1, 20 CPU cores, no MPI (so one rank only). Build:
`cmake -S . -B build-gpu -DCMAKE_BUILD_TYPE=Release -DGEMPIC_CUDA=ON -DCMAKE_CUDA_ARCHITECTURES=120`, HDF5 1.14.6.
This was the first time the GPU code ran on a real GPU.

**ctest of the GPU build: 29 of 30 pass, 1 skipped** (`gui`: needs python ≥ 3.7, the machine has 3.6). On the GPU
(`numerics.backend = cuda`, the default of a CUDA build) ran `physics_cuda`, `regression_cuda`, `restart_cuda`,
`flexgrid_cuda`, `functions_cuda`, `dopant_cuda`, `smoothing_cuda`, `transport_cuda`, `planar_step_cuda`, `lasy`,
`flying_focus`, `undulator`, `cherenkov`, `axis_refine` and `bunch_functions`. Numbers from the log:

| test | GPU result |
|---|---|
| `regression_cuda` | axis 1200 rows and bunches against serial: worst \|diff\| / tolerance 0.129 (tolerance 1e-6 relative) |
| `physics_cuda` | bunch field ⟨E_r⟩/Gauss 1.0002–1.0017 (tol. 0.5 %); dE/E = 7.9·10⁻⁶ over 4 cm at γ = 100; ion density on axis 2.84 / 3.88 / 5.76 / 9.64 / 19.65 (model 2.88 / 3.89 / 5.76 / 9.79 / 20.59); H ionized 0.02294 (ADK 0.02288), N⁵⁺ 0.93568 (rate equations 0.93667), N⁶⁺ 0.05995 (0.05927), N⁷⁺ 0.00437 (0.00406) |
| `restart_cuda`, `flexgrid_cuda`, `functions_cuda`, `dopant_cuda`, `smoothing_cuda` | every compared file agrees (worst \|diff\| / tolerance ≤ 0.242) |
| `planar_step_cuda` | ε_nx at 5 cm: mean-field push 0.0067 µm, level push 0.1179 µm (as serial) |
| `lasy` | file = built-in Gaussian to 1.2·10⁻⁴ of the peak (thetaMode), 9.8·10⁻⁴ (cartesian); flying focus from a file: rms difference 1.28·10⁻³ |
| `axis_refine`, `bunch_functions` | identical to the CPU build: 200 000 → 463 680 macroparticles, σ_x 100.188 → 100.184 µm; ring r = 99.979 µm, width 1.00046 µm after 1 cm |

**`check_1gpu.sh`: 0 failed comparisons.** GPU and serial agree to all printed digits (worst |diff| / tolerance = 0)
for the regression deck (engine and pipelined step: axis, map, bunches), the FEL deck (bunches.txt, fel.txt) and the
flying focus (1e-18). Planar undulator at λ_u/16 on the GPU: ε_nx at 5 cm = 0.0067 µm (serial 0.0067 µm).

**Wall times, GPU against 20 CPU threads** (OpenMP, the same machine):

| case | GPU | 20 threads | GPU speed-up |
|---|---|---|---|
| `lwfa_heavy.deck` (plasma-dominated, ~4·10⁴ plasma particles per slice, 1200 slices, 200 steps) | 75.9 s (0.38 s per step) | 484.2 s (2.42 s per step) | 6.4 |
| FEL deck to 30 mm (beam- and field-dominated, 1.4·10⁶ beam particles, 3 modes, 18 372 slices, 48 steps) | 84.2 s (1.75 s per step) | 75.2 s (1.57 s per step) | 0.9 |

The diagnostics of the timing runs are identical on GPU and CPU (`lwfa_heavy` step 200: |E_x|max = 1.083,
E_z in [−0.0279, 0.0258]; FEL step 48: P₁ in slice 14 = 9.96529·10⁶ W).

`GEMPIC_PROFILE=1` (seconds; the profile synchronizes, so the totals are ~10 % above the plain runs):

| case | streaming push + deposit | field solve | slice preparation | R correction | backup/restore | trapped/freed | wall |
|---|---|---|---|---|---|---|---|
| lwfa, engine | 47.8 | 23.5 | 4.1 | 2.0 | 4.6 | 1.5 | 84.9 |
| lwfa, pipelined (1 rank) | 47.8 | 23.6 | 5.1 | 2.0 | 4.7 | 1.5 | 88.6 |
| FEL, engine | 4.9 | 43.5 | 14.2 | 7.6 | 0.2 | 4.9 | 99.8 |
| FEL, pipelined (1 rank) | 4.9 | 43.5 | 14.3 | 7.7 | 0.1 | 5.0 | 127.3 |

**These timings are not representative of GH200 / A100.** A GeForce GPU computes FP64 at 1/64 of its FP32 rate
(RTX 5090: ~1.6 TFLOP/s FP64, against ~34 for a GH200 and ~10 for an A100, without tensor cores), and the code
computes in double precision. The field solve (tridiagonal solves per mode, FP64) dominates the FEL case and is
where the GeForce is weakest; the plasma push gains even here. In the FEL case the pipelined step needs 27 s more
than the engine, all outside the engine phases (equal in both): the copies of the whole block between host and GPU
every step and the host-side beam bookkeeping (milestone 5e). Real numbers for the scaling series must come from
JUPITER (`benchmarks/jupiter/`).
