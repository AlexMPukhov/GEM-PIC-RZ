# FEL work with GEM-PIC-RZ: handoff (state of 2026-10-21)

This file is for continuing the FEL work in a new Claude session (or by hand). It collects what is in the
repository, what was learned, and what is open. Start a session with: *"Read CLAUDE.md and docs/fel_handoff.md,
then continue with the open items."*

## 1. Goal and message
- Separate publication on FELs with GEM-PIC (the code paper `paper/gempic_rz.tex` only mentions FELs).
  Draft: `paper_fel/gempic_fel.tex` (PRAB style). Status and open items are in `paper_fel/NOTES.md`; open
  points in the text are marked `\todo{}`.
- Message: averaged FEL codes (GENESIS) are accurate for the fundamental and the odd harmonics. Even at large
  rho the period averaging is NOT the limitation (1D test: gain length within 3 % up to rho = 0.15).
  What they leave out by construction are the **even harmonics of planar undulators**, which vanish on the
  axis. GENESIS 4 skips their source term (`src/Core/FieldSolver*.cpp`: `harm % 2 == 1`, "do not need to
  calculate for even harmonics"). GEM-PIC-RZ gets them from first principles.
- Terminology rule (Alexander): never call GEM-PIC(-RZ) "quasi-static" or a "quasi-static approximation"; it
  keeps the slow derivative of the forward field.

## 2. Results so far (all in the repository)
| What | Where | Result |
|---|---|---|
| Helical seeded amplifiers vs GENESIS 4 (10 nm and 0.8 nm) | `benchmarks/genesis_seeded/` | gain length 8.33 / 38.7 cm vs GENESIS with the same physics (LSC + beam widening) 7.92 / 38.4 cm; plain GENESIS 7.45 / 37.4 cm |
| SASE (helical) | `examples/fel_sase.deck`, README "Test: high-gain FEL" | 6.3 periods vs Ming Xie 6.4 |
| Coherent emission per azimuthal mode, planar undulator | `benchmarks/harmonics/` (`emission.deck`, `radint.py`, `compare_emission.py`, `spectra/`, `conv.txt`) | m = 0..3 within 0.3-1.8 % of the radiation integral; converged in dz, dr, dzeta |
| Seeded planar FEL, harmonics 1-3, vs GENESIS 4 | `benchmarks/harmonics/` (`fel.deck`, `genesis/`, `fel_tables/`, `compare_fel.py`, `compare_fel.png`) | P1 sat 10.4 GW (GENESIS 9.8 GW), Lg 7.9 vs 7.0 cm; **P2 = 9.6 MW (P2/P1 = 0.9e-3), GENESIS 0**; P3 2.1 % vs 2.8 % of P1; b1/b2/b3 0.63/0.18/0.14 vs 0.66/0.25/0.13 |
| 1D averaged vs unaveraged, rho scan | `benchmarks/high_rho/` (`fel1d.py`, `scan.py`, `scan.txt`) | Lg within 3 % up to rho = 0.15, any K; Psat of the averaged model low by ~1.2 rho |
| Paper figures | `paper_fel/scripts/fig_*.py` -> `paper_fel/figs/` | regenerate from the data above |
| GPU test package | `benchmarks/gpu/` | not run yet (no GPU here); includes 96 steps of the FEL deck |

Details and reproduction commands: `benchmarks/README.md`, sections "Seeded FEL against GENESIS 1.3 v4" and
"Harmonics of a planar undulator".

## 3. FEL-relevant code (where things are)
- `src/external.hpp`: undulator fields (planar with/without natural focusing, helical, taper K(z), linear
  ramps), `sincHalf` kick correction (correction 7), `velocityAtPosition` for the current deposit (correction 8).
  Corrections 7 and 8 are described in README "Corrections"; without them the gain length depends on dz to
  first order (9.5 -> 8.3 cm at 10 nm).
- `src/simulation.cpp` `initBunches`: Fawley beamlets (`shot_noise = fawley`, `beamlet` = macroparticles per
  beamlet, quiet for harmonics h < beamlet; use `beamlet = 8` for harmonics up to 3), finite-gamma
  co-moving bunch fields at the start.
- `src/diagnostics.cpp`: `fel.txt`, `fel_power.txt` (demodulated power at lambda, modes m >= 1, per slice) and
  **`fel_harmonics.txt`** (new: `output.fel_harmonics = N`, power per harmonic h and mode m = 0..m_max per slice,
  plus bunching per harmonic, line m = -1). Power of mode m >= 1 = pi * sum A (|R+-|^2 - |L+-|^2), m = 0 half.
  Merged over MPI ranks (`src/pipeline.cpp` putDiag/getDiag); checked in `tests/undulator.sh`.
- `src/config.cpp`: warning for planar undulators with dz > period/32 and a bunch >= 1 pC.
- `tools/h5_extract.cpp`: `spectrum zlo zhi k1 k2 nk` now prints dW/dk per mode m = 0..M after the m >= 1 total.
- `tools/undulator_radint.py` (test reference) and `benchmarks/harmonics/radint.py` (radiation integral split
  by parity and by azimuthal mode; ~3 min).

## 4. Lessons / pitfalls (learned the hard way)
1. **Planar undulators need dz <= lambda_u/32** with intense beams. At lambda_u/16 a 4 kA beam blew up
   transversely from the tail forward (eps_n 0 -> 2.8 um in 25 periods; lambda_u/24: 0.13, /32: 0.04,
   /64: 0.02 um), and the gain length came out 1.8x too long. More beam sub-steps (`beam_dz_max`) or
   `iterations = 2` did not help: the field step is the cause. Helical at /16 is far less affected
   (eps 0.33 um in the same test). Root cause NOT found. Hypothesis: the s-discretization of the large co-moving
   field of the oscillating transverse current j_x ~ cos(k_u s) breaks the E + v x B cancellation at second
   order in k_u dz. Cheap test setup: `fel.deck` with `grid.z_end=0.25m grid.dzeta=1.328nm
   bunch.macroparticles=128000` (~3 min), watch `eps_nx` in `bunches.txt`.
2. **Window ahead of a short bunch for spontaneous emission:** radiation at angle theta advances in zeta by
   theta^2 (L - z)/2 (11 um at 1.5/gamma for gamma = 100, 10 periods). Too short a gap ahead loses the
   large-angle part (m = 2, 3 came out 12-14 % low). FEL radiation is narrow-angle, so it is not affected.
3. **fel_slice must be a whole number of lambda_r**, otherwise the slice sum of a uniform beam gives a spurious
   bunching ~ 1/(number of wavelengths) (0.008 for 18.8 lambda_r).
4. **Radial absorber inside gamma*L eats the co-moving field** of the bunch: keep `r_absorber >= r_max` (off)
   for FEL runs unless needed. (It was not the cause of item 1, but it drains energy from the beam.)
5. **m = 0 at even harmonics also holds the co-moving field of the bunching.** The 2nd harmonic radiation
   (E_x ~ x) has P(m=0) = P(m=2), so P2 = 2 P(m=2) is the clean measure.
6. **Parity:** odd harmonics -> odd azimuthal modes, even harmonics -> even modes (rotation by pi = half-period
   time shift). m_max = 2 contains the 2nd harmonic; the 3rd needs m_max = 3 for its off-axis part.
7. **Resolution:** dzeta about 40 points per harmonic wavelength for ~2 % (emission test: 21 points -> -4 %).
8. **Memory:** about 800 bytes per (radial node x slice x mode). `fel.deck` (122 x 18 372 x 3) takes 5.6 GB;
   the container had 7 GB, which is why the 3rd-harmonic rerun (dzeta = lambda_r/48, m_max = 3) is still open.
9. **Run times (2 cores, OpenMP):** `fel.deck` 1.8 h (3200 steps); emission deck 11 min; GENESIS steady state
   30 s.
10. The steady-state slices of `fel.deck` are zeta = 2.76 ... 3.40 um (radiation that met only the flat top);
    `compare_fel.py` uses them.

## 5. GENESIS 4 build (reference code)
- Source: github.com/svenreiche/Genesis-1.3-Version4, release 4.6.15 (commit 35f2e14), built against serial HDF5
  with OpenMPI. The serial HDF5 lacks the MPI-IO calls, so a shim header was force-included:
  `cmake -DCMAKE_CXX_FLAGS="-include $PWD/benchmarks/genesis_shim.h" ..` (file in the repo). Run on 1 rank.
- Harmonics: extra `&field` blocks with `harm = 2`, `harm = 3` (power 0), `nbins = 8`, `&track bunchharm = 3`
  (see `benchmarks/harmonics/genesis/fel.in`). Planar undulator: `aw = K/sqrt(2)`, `helical = false`,
  `kx = 0, ky = 0` for no natural focusing.

## 6. Open items (suggested order)
1. Independent check of P2: Puffin (unaveraged, radiates even harmonics) on the seeded planar case, or Huang-Kim
   theory (PRE 62, 7295 (2000)).
2. Planar step instability (section 4, item 1): cause and fix; a fix halves the cost.
3. 3rd harmonic with dzeta = lambda_r/48, m_max = 3 (needs ~8 GB), and the post-saturation difference to GENESIS.
4. Transverse profile of the 2nd harmonic (zero on the axis) as a paper figure: rerun `fel.deck` with an HDF5
   dump near saturation (step ~3000), demodulate at 2 k_r per mode.
5. GENESIS planar case with the same physics (LSC via `&efield`, beam widening) to explain the 13 % in gain length.
6. Cost vs time-dependent GENESIS on the same bunch.
7. Optional: superradiant short bunch; harmonics with OAM (helical); the combined LPA -> transport -> FEL run
   planned by Alexander (two chained runs: LWFA with down-ramp, beam written at the plasma exit; then the
   witness re-imported with quiet start + shot noise, plasma lens or solenoid focusing (quadrupoles need
   m_max >= 2), thin chicane as R56 map, undulator; zeta grid lambda_r/16 around the witness, free-hand z grid).
8. References marked `\todo{check}` in the draft (Freund 2000, Tremaine 2002), the GENESIS 4 citation, arXiv
   number of the code paper.

## 7. Housekeeping
- The repository is public (github.com/AlexMPukhov/GEM-PIC-RZ): everything pushed is visible.
- Big outputs are not in git (`benchmarks/harmonics/em_*/`, `fel_out*/`, HDF5 dumps of 0.7 GB each); the thinned
  tables in `fel_tables/` and the spectra in `spectra/` are enough for all figures.
- Build and tests: `cmake -B build && cmake --build build -j && (cd build && ctest)` (26 tests, ~20 min on 2 cores).
