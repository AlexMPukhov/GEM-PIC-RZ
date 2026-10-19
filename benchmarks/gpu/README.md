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
| 4 | `../harmonics/fel.deck` (seeded planar FEL, 640 000 beam particles, 3 modes × 122 radial × 18 372 slices; needs ~6 GB), first 96 steps | speed, beam- and field-dominated; P₁ in slice 14 at step 96 should equal the CPU value 1.0006·10⁷ W to ~1e-4 | 2.1 s per step |

What to expect / look for:
- The ζ-march is sequential: each slice launches about 10 kernels. In the 100 µm deck (~1600 plasma particles
  per slice) the GPU is launch-latency bound and may hardly beat 2 CPU cores; the heavy deck shows the real
  speed-up. If step 3 is also slow, the `GEMPIC_PROFILE` output tells whether the field solve (one block per mode,
  tridiagonal solves in one thread) dominates; that would be the place to bring in a parallel tridiagonal
  solver, as was done in QUARZ.
- The GPU engine is not combined with MPI yet: one GPU per run.
- Please send back: the `nvidia-smi` line, the ctest summary, the `compare.py` lines of step 2, the `done:` lines
  and the profile of steps 3 and 4.
