# GEM-PIC-RZ: notes for Claude sessions

Read this first. README.md is the user documentation; docs/ holds the details.

- Owner: Alexander Pukhov (HHU Duesseldorf). Repository is PUBLIC: everything pushed is visible.
- Never call GEM-PIC(-RZ) "quasi-static" or a "quasi-static approximation": it keeps the slow derivative of the
  forward field (quasi-static = dropping ALL slow-time derivatives).
- Build/test: `cmake -B build && cmake --build build -j && (cd build && ctest)`; for many threads use
  `numerics.backend=openmp` and `OMP_PROC_BIND=false` in containers.
- Unsubmitted papers that are not directly about the code go to the private repository CLAUDE-notes, never here.
- Commit at milestones with descriptive messages. Larger features go on a branch and are merged into main by a pull
  request after the tests pass; the owner prefers that Claude pushes and opens the pull requests itself.
- Every new feature gets a test (ctest) against an independent reference and an entry in docs/ and README.md.
- Testing aids: `-DCMAKE_CXX_FLAGS=-DGEMPIC_RANGE_CHECK` (every field/current access checked against the stored
  slices), an AddressSanitizer build, and `tests/mpi.sh <exe> <compare_columns> mpirun` for 2-4 ranks.
- **Setting up runs with a flexible radial grid (`radial_grid = table`): the grid must not reflect radiation back**
  (Alexander). (1) A cell larger than the transverse wavelength cannot carry the wave (discrete cutoff
  k_perp = 2/dr): it is thrown back. (2) Below the cutoff the numerical dispersion, k_perp,eff = (2/dr) sin(k_perp dr/2),
  falls as dr grows outward and acts like a lower refractive index: numerical total internal reflection of grazing
  waves. Estimate the largest k_perp = k theta that must leave the beam region (radiation cones of all relevant
  harmonics, also off-axis harmonics landing at the diagnosed frequency, e.g. for a planar undulator harmonic h' seen
  at h: gamma theta = sqrt((1 + K^2/2)(h'/h - 1)); diffraction of the smallest features), keep k_perp dr <~ 1 wherever
  it must propagate freely, let dr grow only outside that region or behind an absorber that removes it first, and
  stretch slowly. Symptom of a violation: fringes near the axis with a fixed radial wavenumber that disappear when the
  fine region is extended.
- Compile-check the CUDA build after every engine change, also without a GPU (nvcc catches what g++ does not,
  e.g. device lambdas in private member functions): `-DGEMPIC_CUDA=ON -DCMAKE_CUDA_ARCHITECTURES=90
  -DCMAKE_CUDA_COMPILER=<nvcc>` (nvcc from tools/install_cuda_pip.sh).
- Ongoing work:
  - **FEL** (separate paper): draft, notes and handoff are in the PRIVATE repository AlexMPukhov/CLAUDE-notes
    (`gem-pic-rz/`); read them there if that repository is in the session. Benchmarks here:
    `benchmarks/harmonics/` (emission test, FEL deck), `benchmarks/genesis_seeded/`; the FEL result data are private too.
  - Code paper: arXiv:2610.31622 (v1 source in `paper/arxiv_v1/`); journal version for J. Comput. Phys. in
    `paper/gempic_rz.tex` (elsarticle, highlights, Elsevier end statements). Code DOI 10.5281/zenodo.23600178.
  - GPU check: `benchmarks/gpu/` (gpu_check.sh, check_1gpu.sh; first real-GPU results on an RTX 5090 in the README there).
  - **Flying focus** (luminal Bessel pulse fed through the radial boundary, `flying_focus` block): handover in
    `docs/flying_focus_handoff.md`; tests and production decks in `benchmarks/flying_focus/` (README there).
  - **LASY** laser files (merged, PR #1): `laser.file` (import, exact propagation, `inject = boundary`) and
    `output.lasy` (export), `src/laser_file.*`, ctest `lasy`, docs/input_deck.md. Motivation: flying focus /
    structured light prepared in LASY. ctest lasy passes on a real GPU (RTX 5090); MPI checked in tests/mpi.sh (3 ranks = serial).
  - **Multi-GPU plan** (MPI pipeline + GPU engine), branch `distfields`:
    1. done: distributed storage in the MPI pipeline: `Fields`/`Currents` store only [n0, n0 + ns) with
       global-index accessors; every rank keeps its block and 2 pf + 6 ghost slices (`pipelineBlock()` in
       pipeline.cpp, `StorageBlock` in simulation.hpp).
    2. done: every rank BUILDS only its block (lasers, LASY lasers, flying focus, bunch self-fields: for v < c the
       zeta solve is a Thomas sweep across the ranks, bit-identical to one rank); restart reads the block from the
       checkpoint (hyperslabs); the initial output goes through the pipeline. Results bit-identical to
       `GEMPIC_FULL_STORAGE=1`; tests/mpi.sh passes with range checks.
    3. done: block-wise output (Simulation::beginOutput / outputBlock / endOutput with one sink per format, HDF5
       hyperslabs); the pipeline streams the blocks of the ranks through rank 0 (at most two blocks there). Files
       of 3 ranks = serial files (h5diff: only date / input_deck differ).
    4. done (branch `mpi-mean-field`): mean-field beam push (`beam_field = mean`) in MPI: `PipelineRun::meanPhase`
       does the second push of step k at the start of step k+1 (F^1 below the block from the R message, the front
       rank's results in the F message); tags in Particle::id pair the two pushes across ranks; F^1 of the absorber
       layer is used before the absorber, as serially. tests/mpi.sh: 1 rank pipelined = serial (exact), 2-4 ranks
       and a 3 -> 2 restart = serial; planar undulator at lambda_u/16 on 3 ranks: eps_nx 0.0067 um as serial.
       Merged (PR #2). dE/E in bunches.txt is now accumulated centred (Moments::gmean, gM2; branch
       centred-energy-spread): the old sum gamma^2 - (sum gamma)^2 / Q was 27 % low at dE/E = 4e-6, gamma = 1000,
       and differed between rank counts; ctest energy_spread checks it against the particle data.
    5. done (PRs #4-#6, merged 2026-10-24/25): the GPU engine on the block of an MPI rank, backend = cuda (or
       GEMPIC_PIPELINE_ENGINE=1 with backend = openmp, the GPU emulation build, for testing); one GPU per rank
       from the node-local rank. 5a device arrays in the host layout with slice offsets (FieldsV n0/ns, CurrV
       k0/nk; whole window bit-identical). 5b the zeta march of a block (Engine::marchBlock, Engine::march shared
       with step(), absorbSlices per final slice with pre-absorber values; PipelineRun::deviceMarch). 5c the beam
       push of the own list, the continued pushes and the halo copies (Engine::pushBeams; host bookkeeping on the
       returned states), deposits into c1b / H and the filter across the block edges (deviceBeamCurrent with
       Engine::combine/filter/zero/addCurrents, AUX_H / AUX_RW). 5d the mean phase (pushBeams with wMean = 0.5
       after uploading F^1 with the pre-absorber values; deviceBeamCurrent(C0B, ..., meanPhase = true)).
       Emulation = serial on all tests/mpi.sh decks; CUDA build compiles (nvcc 13.4, sm_90). Open, 5e: transfers
       per step are still whole blocks (fold, c0b, c0s up; fnew, c1b, c1s down; the beam list up and down) and
       the MPI messages go through the host (CUDA-aware MPI later) - decide after the first JUPITER timings.
    6. real-GPU tests and the scaling series (1, 4, 16, 64 GPUs) on JUPITER / JUWELS Booster, run by Thomas
       Wilson: scripts in benchmarks/jupiter/ (README there: build.sh, check.sbatch = GPU vs serial on 1-4 GPUs,
       submit_scaling.sh + scaling.sbatch + cpu_reference.sbatch, collect.sh = table). Dry-run in WSL with
       stand-ins for srun/sbatch/nvidia-smi and the emulation binary. Results 2026-10-26 in
       benchmarks/jupiter/README.md (35/35 ctests, 1-4 GPUs = serial; 4/16/64 GPUs: LWFA 88/66/28 %, FEL 64/47/30 %
       efficiency; 1 GH200 node = 2.6x (LWFA) / 1.4x (FEL) one CPU node). Next: 5e (device-resident blocks,
       CUDA-aware MPI) and a weak-scaling series.
       6a done (2026-10-25, branch gpu-rtx5090): first real-GPU run, one RTX 5090 (cc 12.0, CUDA 13.0, no MPI, so
       1 rank): ctest of the GPU build 29/30 pass + gui skipped (python 3.6); benchmarks/gpu/check_1gpu.sh (the
       1-rank part of check.sbatch, no Slurm) 0 failed comparisons, GPU = serial to all printed digits (engine and
       pipelined block march, FEL mean-field push, flying focus), planar eps_nx 0.0067 um. Fixed: axis_refine /
       bunch_functions ran backend = openmp, rejected by a CUDA build (now the parallel backend of the build).
       Timings vs 20 threads: lwfa_heavy 6.4x faster, FEL 0.9x (FP64-bound field solve: GeForce FP64 = 1/64,
       not representative of GH200); FEL pipelined +27 s over the engine outside the engine phases (5e transfers).
       Still open: several GPUs (MPI) and representative timings on JUPITER.
  - **Thin, intense, low-gamma beams** (branch `beam-self-field-2`; README "Thin, intense beams"; ctest thin_ring):
    `numerics.beam_consistent` ((1/4,1/2,1/4) beam sources in the half-cell equations, charge-conserving radial
    flux in the R update, all modes) and `numerics.theta_stiff` (theta per radial eigenmode of E_z, B_z:
    pentadiagonal implicit matrices, FieldSolver::St; engine: leader-sequential solves in solveSliceTeam). Off by
    default (bit-identical when off). Open: beam_consistent in the MPI pipeline (needs c1b at half-cell b and the
    flux across block edges); team-parallel band solves on the GPU. Private context: CLAUDE-notes
    gem-pic-rz/self_field/HANDOVER.md.
  - **Metal wall** (branch `metal-wall`; block `wall`, README "Metal walls", ctest wall_*): perfectly conducting
    wall r = r_w(z) inside the grid (FieldSolver::setWall: last face K, ghost centre of E_z at r_w, zero B_z flux,
    fields outside zero; absorber off). Sloped wall: E_z = -r_w' E_r at the wall (kept in e[K]); widening wall: the
    field continues into the opened annulus (Simulation/Engine::extendWall). Serial, engine, pipeline and MPI.
    Against the modal theory of the tube (tests/wall.py, the scheme's dispersion relation): 0.8 % straight, 1.0 %
    12 -> 9 um taper. Open: varying radius with theta_stiff; reflections (backward waves) are outside the model.
    Asked for by Alexander (waveguides of variable radius); dielectric structures were discussed, not started.
  - Open physics issue: density down-ramp injection vs FBPIC (attributed to the GT2 choice; GT1 future work).
