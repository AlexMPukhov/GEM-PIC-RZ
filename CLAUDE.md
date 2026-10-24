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
- Ongoing work:
  - **FEL** (separate paper): draft, notes and handoff are in the PRIVATE repository AlexMPukhov/CLAUDE-notes
    (`gem-pic-rz/`); read them there if that repository is in the session. Benchmarks here:
    `benchmarks/harmonics/` (emission test, FEL deck), `benchmarks/genesis_seeded/`; the FEL result data are private too.
  - Code paper: arXiv:2610.31622 (v1 source in `paper/arxiv_v1/`); journal version for J. Comput. Phys. in
    `paper/gempic_rz.tex` (elsarticle, highlights, Elsevier end statements). Code DOI 10.5281/zenodo.23600178.
  - GPU check: `benchmarks/gpu/` (not yet run on a GPU).
  - **Flying focus** (luminal Bessel pulse fed through the radial boundary, `flying_focus` block): handover in
    `docs/flying_focus_handoff.md`; tests and production decks in `benchmarks/flying_focus/` (README there).
  - **LASY** laser files (merged, PR #1): `laser.file` (import, exact propagation, `inject = boundary`) and
    `output.lasy` (export), `src/laser_file.*`, ctest `lasy`, docs/input_deck.md. Motivation: flying focus /
    structured light prepared in LASY. Not yet run on a real GPU; MPI checked in tests/mpi.sh (3 ranks = serial).
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
       Known: dE/E in bunches.txt differs by up to 0.4 % between rank counts at gamma = 1000 (cancellation in
       <g^2> - <g>^2; also with level and on main) - a centred variance would fix it (owner's decision: it changes
       printed serial values in the last digits).
    5. next: the GPU engine on a block (march over [a, b), the pipeline exchanges staged through the host, later
       CUDA-aware MPI); develop with the GPU-emulation build against the serial code.
    6. then: real-GPU tests (`ctest -R lasy`, regression_cuda, physics_cuda) and a scaling series (1, 4, 16, 64
       ranks/GPUs) on JUWELS/JUPITER, run by a group member.
  - Open physics issue: density down-ramp injection vs FBPIC (attributed to the GT2 choice; GT1 future work).
