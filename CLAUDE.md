# GEM-PIC-RZ: notes for Claude sessions

Read this first. README.md is the user documentation; docs/ holds the details.

- Owner: Alexander Pukhov (HHU Duesseldorf). Repository is PUBLIC: everything pushed is visible.
- Never call GEM-PIC(-RZ) "quasi-static" or a "quasi-static approximation": it keeps the slow derivative of the
  forward field (quasi-static = dropping ALL slow-time derivatives).
- Build/test: `cmake -B build && cmake --build build -j && (cd build && ctest)`; for many threads use
  `numerics.backend=openmp` and `OMP_PROC_BIND=false` in containers.
- Unsubmitted papers that are not directly about the code go to the private repository CLAUDE-notes, never here.
- Commit at milestones with descriptive messages; push to origin main.
- Ongoing work:
  - **FEL** (separate paper): draft, notes and handoff are in the PRIVATE repository AlexMPukhov/CLAUDE-notes
    (`gem-pic-rz/`); read them there if that repository is in the session. Benchmarks here:
    `benchmarks/harmonics/` (emission test, FEL deck), `benchmarks/genesis_seeded/`; the FEL result data are private too.
  - Code paper: arXiv:2610.31622 (v1 source in `paper/arxiv_v1/`); journal version for J. Comput. Phys. in
    `paper/gempic_rz.tex` (elsarticle, highlights, Elsevier end statements). Code DOI 10.5281/zenodo.23600178.
  - GPU check: `benchmarks/gpu/` (not yet run on a GPU).
  - **Flying focus** (luminal Bessel pulse fed through the radial boundary, `flying_focus` block): handover in
    `docs/flying_focus_handoff.md`; tests and production decks in `benchmarks/flying_focus/` (README there).
  - **LASY** laser files (branch `lasy`): `laser.file` (import, exact propagation, `inject = boundary`) and
    `output.lasy` (export), `src/laser_file.*`, ctest `lasy`, docs/input_deck.md. Motivation: flying focus /
    structured light prepared in LASY. Not yet run on a real GPU; MPI checked in tests/mpi.sh (3 ranks = serial).
  - **Multi-GPU plan** (MPI pipeline + GPU engine): step 1 done = distributed field storage and per-block initialization/restart in the MPI pipeline
    (`Fields`/`Currents` store [n0, n0 + ns) with global-index accessors; `-DGEMPIC_RANGE_CHECK` checks every
    access). Next: output and checkpoints without assembling the window on rank 0, then the mean-field beam push (`beam_field = mean`) in MPI, then the block-restricted GPU march with the
    pipeline exchanges, then a scaling test on JUPITER.
  - Open physics issue: density down-ramp injection vs FBPIC (attributed to the GT2 choice; GT1 future work).
