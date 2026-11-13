# GEM-PIC-RZ on JUPITER / JUWELS Booster: first real-GPU runs and the scaling series

The MPI pipeline with `numerics.backend = cuda` gives every rank one GPU. The rank pushes its beam, marches its
block of slices and does the mean-field second push there (README, MPI section). So far this was tested only in
the GPU emulation build (the GPU code paths on CPU threads) and compiled with nvcc 13.4 for sm_90. **On one GPU
it has now run on a real GPU** (RTX 5090, one rank, no MPI, 2026-10-25): the ctest of the GPU build and the
one-rank comparisons of `check.sbatch` (`benchmarks/gpu/check_1gpu.sh`) agree with the serial code, see
`benchmarks/gpu/README.md`. **Several GPUs (MPI) have never run, and there are no representative timings yet**
(GeForce FP64 is 1/64 of FP32). These scripts do that, in two parts:

| part | script | resources | what it shows |
|---|---|---|---|
| 1 | `check.sbatch` | 1 node, 4 GPUs, ≤ 1.5 h | the GPU results equal the serial code (1e-6), on 1, 2 and 4 GPUs |
| 2 | `submit_scaling.sh` | 1, 4, 16, 64 GPUs (1, 1, 4, 16 nodes) + 1 CPU node | time per step, speed-up, parallel efficiency |

Run part 2 only after part 1 has passed.

## 0. Once: get the code and build it

On a **login node of the system you run on**. JUPITER is aarch64 (Grace + GH200); JUWELS Booster is x86 + A100.

```
git clone https://github.com/AlexMPukhov/GEM-PIC-RZ.git && cd GEM-PIC-RZ
# edit benchmarks/jupiter/env.sh: ACCOUNT (your compute project); on JUWELS Booster also CUDA_ARCH=80
sh benchmarks/jupiter/build.sh
```

`build.sh` loads `MODULES` (default `GCC OpenMPI CUDA HDF5 CMake`, as for QUARZ on JUPITER) and builds two copies:

- **`build-gpu`:** CUDA, MPI, HDF5. Its ctest uses `srun`.
- **`build-cpu`:** serial references and the CPU timing.

If a module has another name, set `MODULES` (find the names with `module spider CUDA` and `module spider HDF5`).

## 1. Correctness on the GPUs

From the GEM-PIC-RZ directory (the jobs find the code through the submit directory):

```
sbatch --account=<project> benchmarks/jupiter/check.sbatch      # --partition=develbooster for a quicker start
```

The log `gempic-check.<job>.log` contains:

1. **`nvidia-smi` and ctest of `build-gpu`.** All tests should pass, including:
   - `regression_cuda`, `physics_cuda`, `restart_cuda` (one GPU);
   - `mpi_pipeline`, whose line "device march of the blocks (numerics.backend=cuda)" means 1–3 ranks on the GPUs equal the serial run.
2. **The MPI pipeline on 1, 2 and 4 GPUs against the serial run** (`build-cpu`), one `OK` line per file:
   - the regression deck (laser wake, trapping, ions);
   - the FEL deck (4 kA beam across the blocks, mean-field push);
   - the flying focus through the absorber layer.
3. **The 4 kA beam in a planar undulator at λ_u/16 on 4 GPUs:** ε_nx at 5 cm should be 0.0067 µm (serial: 0.0067 µm; without the mean-field push it would be 0.118 µm).

The last line is `summary: 0 failed comparison(s)` when all is well.

## 2. Scaling

```
sh benchmarks/jupiter/submit_scaling.sh          # 1, 4, 16, 64 GPUs and the CPU reference: 5 jobs
sh benchmarks/jupiter/collect.sh                 # the table, when the jobs are done
```

**Two test cases, each with `GEMPIC_PROFILE=1`** (time per phase of the GPU engine, printed per rank):

- **lwfa:** `benchmarks/gpu/lwfa_heavy.deck`. LWFA in a plasma, about 4·10⁴ plasma particles per slice, 1200 slices, 200 steps; dominated by the plasma.
- **fel:** `benchmarks/harmonics/fel.deck` to 25 cm. 4 kA beam, 640 000 particles, 3 modes, 18 372 slices, 400 steps; dominated by the beam and the field solves.

**With 1 GPU, two runs per case:**

- **`*_engine`:** the single-GPU engine (whole window, no pipeline).
- **plain:** the pipelined step on one rank, which is the baseline of the speed-up.

**The CPU reference (`*_cpu`):** the same decks on one node without GPUs, 4 ranks × 72 OpenMP threads. On JUWELS Booster set `RANKS=4 THREADS=12`.

`collect.sh` prints for every run:

- time per step;
- speed-up and efficiency against 1 GPU;
- "pipeline": the ideal efficiency of the pipeline, N/(N+P−1) for N steps on P ranks. The measured efficiency cannot be higher.

**Short trial before the full series:**
`EXTRA=grid.z_end=20um sh benchmarks/jupiter/submit_scaling.sh 1 4` (8 steps of lwfa, 1 step of fel).

**Budget, upper limits** (set by `--time` in `submit_scaling.sh`):

- 1 GPU: 4 h
- 4 GPUs: 2 h
- 16 GPUs: 1 h
- 64 GPUs: 1 h

That is at most about 100 GPU-hours, plus 4 node-hours for the CPU reference. The real times should be much shorter. The `done:` lines give them.

## Results on JUPITER (Thomas Wilson, 2026-10-26; `results_2026-10-26/`)

**Correctness.** All 35 ctests of the GPU build pass. The regression, FEL (mean-field push, beam across the blocks)
and flying-focus decks on 1, 2 and 4 GPUs equal the serial code (0 failed comparisons). The planar undulator at
λ_u/16 on 4 GPUs gives ε_nx = 0.0067 µm, as serially.

**Strong scaling** (s per step; speed-up and efficiency against the pipelined step on 1 GH200; "pipeline" is the
ideal N/(N + P − 1)):

| GPUs | LWFA (`lwfa_heavy`, 200 steps) | speed-up | eff. | pipeline | FEL (400 steps) | speed-up | eff. | pipeline |
|---|---|---|---|---|---|---|---|---|
| 1 (engine, no pipeline) | 0.366 | 1.02 | | | 1.894 | 1.14 | | |
| 1 | 0.371 | 1 | 100 % | 100 % | 2.151 | 1 | 100 % | 100 % |
| 4 | 0.105 | 3.5 | 88 % | 99 % | 0.842 | 2.6 | 64 % | 99 % |
| 16 | 0.035 | 10.6 | 66 % | 93 % | 0.284 | 7.6 | 47 % | 96 % |
| 64 | 0.020 | 18.1 | 28 % | 76 % | 0.112 | 19.2 | 30 % | 86 % |
| CPU node, 4 × 72 cores | 0.273 | 1.36 | | | 1.166 | 1.84 | | |

Reading:
- One GH200 node (4 GPUs) is 2.6× (LWFA) and 1.4× (FEL) faster than one CPU node.
- **The GPU part scales ideally.** The `engine profile` of rank 0 (the march on the GPU) is 609 / 150 / 38 / 9.4 s
  for FEL and 73 / 18.3 / 4.7 / 1.18 s for LWFA on 1 / 4 / 16 / 64 GPUs, i.e. exactly 1/P. All the losses are
  outside the march: wall minus engine time is 251 / 187 / 75 / 35 s (FEL) and 1 / 2.7 / 2.3 / 2.9 s (LWFA).
- **Why the host part hurts more than its share.** Rank r can start step k only after rank r+1 has finished it, and
  rank r+1 can start step k+1 only after rank r has marched the first slices of step k (the ghosts). In the steady
  state a step therefore takes the time of a whole step **plus** everything a rank does before its first slices:
  receiving the messages, the mean-field second push, the per-step diagnostics, the beam push, the uploads of the
  whole block. On one GPU this prefix counts once, with several GPUs twice. It is largest for the FEL deck (640 000
  beam particles on the ranks that hold the beam, the per-step FEL diagnostics of 18 372 slices on the host).
- Rank 0 receives the FEL slice tables of every rank every step (5 MB per rank, mostly zeros outside its block).
- At 64 GPUs the LWFA window (1200 slices) leaves about 19 slices per rank: strong scaling ends there. A weak-scaling
  series (window or step count growing with the GPU count) would be the fairer test of large runs.

## Round 2 (2026-10-26, branch `device-resident`): what changed, what to run

Changes against the code of the first series (results bit-identical in the emulation build, 1–4 ranks):
- **Level s⁰ and its currents stay on the GPU** between the steps. The host uploads only the ghost slices it
  receives from the rank behind; the mean-field push swaps the pre-absorber values of the absorber layer on the
  GPU instead of uploading the whole block twice. Per step and rank this removes three whole-block uploads (fields
  of s⁰, beam and plasma currents of s⁰), two more with the mean-field push.
- **The per-step diagnostics run on a host thread** while the GPU marches the next step (they read the host copy of
  the previous level, which the step does not change), and rank 0 writes the tables of step k after its step k+1.
  Output and checkpoint steps stay synchronous. The diagnostics are no longer in the prefix described above.
- **The FEL slice tables go to rank 0 as the non-zero range of each rank** (the slices of its block).
- **`GEMPIC_PROFILE=1` prints a `pipeline profile`** on rank 0: wait (messages), beam push+current, upload, march,
  download, messages/swap, diagnostics+output, finite check, mean phase; rank 0 and the maximum over the ranks.
  The line after it gives the **time per step after the pipeline has filled** (rank 0), the number that matters
  for long runs; `scaling.sbatch` stores it as column 9 of `scaling.txt`, `collect.sh` lists it.
- **Weak scaling:** `submit_weak.sh` (`weak.sbatch`): the LWFA deck with 1200 slices per GPU (window × P, laser
  at the same distance from the front), 200 steps, 1 / 4 / 16 / 64 GPUs, each ≤ 10 min. `collect.sh` prints the
  table (efficiency = time per step after the fill on 1 GPU / on P GPUs).

To run (after `git pull` and `sh benchmarks/jupiter/build.sh`):

```
sbatch --account=<project> benchmarks/jupiter/check.sbatch      # correctness first, as before
NO_CPU=1 sh benchmarks/jupiter/submit_scaling.sh                # strong scaling again (without the CPU reference)
sh benchmarks/jupiter/submit_weak.sh                             # weak scaling, 1 / 4 / 16 / 64 GPUs
sh benchmarks/jupiter/collect.sh
```

Send back the check log, `collect.sh` output, `runs/scaling.txt`, `runs/weak.txt`, and from every
`runs/scale_*/*.log` and `runs/weak_*/lwfa.log` the last lines (`engine profile`, `pipeline profile`, `pipeline: ...
per step`): `tail -n 5 runs/*/*.log`.

## What to send back

- `gempic-check.<job>.log`, and if a comparison failed, the directory `benchmarks/jupiter/runs/check.<job>`.
- The output of `collect.sh` and the file `benchmarks/jupiter/runs/scaling.txt`.
- One `engine profile` line per run from `benchmarks/jupiter/runs/scale_*/{lwfa,fel}.log` (rank 0 is enough). It tells where the GPU time goes.
- Any error at the start (module, `no usable CUDA device`, MPI). `docs/gpu.md`, section 5, explains the messages.

## Known limits of this version

- **Copies between CPU and GPU every step, for every rank** (since round 2 only downwards for the fields):
  - level s¹ and the currents come down after the march (diagnostics, messages, output);
  - the beam list goes up and down for the push.

  The `pipeline profile` shows their cost (upload / download / beam push).
- **MPI messages go through the host.** CUDA-aware MPI is not used yet.
- **One GPU per rank.** The rank picks its GPU from `SLURM_LOCALID`, or uses the single GPU it sees with `--gpus-per-task=1`.
