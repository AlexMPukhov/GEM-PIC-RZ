# GEM-PIC-RZ: notes for Claude sessions

Read this first. README.md is the user documentation; docs/ holds the details.

- Owner: Alexander Pukhov (HHU Duesseldorf). Repository is PUBLIC: everything pushed is visible.
- Never call GEM-PIC(-RZ) "quasi-static" or a "quasi-static approximation": it keeps the slow derivative of the
  forward field (quasi-static = dropping ALL slow-time derivatives).
- Build/test: `cmake -B build && cmake --build build -j && (cd build && ctest)`; for many threads use
  `numerics.backend=openmp` and `OMP_PROC_BIND=false` in containers.
- Commit at milestones with descriptive messages; push to origin main.
- Ongoing work:
  - **FEL** (separate paper): `docs/fel_handoff.md` (state, lessons, open items), `paper_fel/` (draft + NOTES.md),
    `benchmarks/harmonics/`, `benchmarks/high_rho/`, `benchmarks/genesis_seeded/`.
  - Code paper: `paper/gempic_rz.tex`.
  - GPU check: `benchmarks/gpu/` (not yet run on a GPU).
  - Flying focus (luminal Bessel pulse fed through the radial boundary, `flying_focus` block): first tests in
    `benchmarks/flying_focus/` (README there); production runs planned on Alexander's cluster.
  - Open physics issue: density down-ramp injection vs FBPIC (attributed to the GT2 choice; GT1 future work).
