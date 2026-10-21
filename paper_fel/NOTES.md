# FEL manuscript: status and open items (started 2026-10-18)

See also `docs/fel_handoff.md` (full handoff: results, code locations, lessons, open items).

Working title: *First-principles simulation of free-electron lasers with the Galilean electromagnetic PIC method:
even harmonics of planar undulators*. Target: PRAB (like the code paper). Source: `gempic_fel.tex`, open
points marked `\todo{}` in red.

## Message
Averaged FEL codes are accurate for the fundamental and the odd harmonics. Even at large rho, period averaging
is not what limits them (Sec. VI: gain length within 3 % up to rho = 0.15). What they leave out by construction
is the even harmonics of planar undulators, which are zero on the axis. GENESIS 4 skips their source term.
GEM-PIC-RZ gets them from first principles, at moderate cost.

## What is in the draft (results already exist)
| Section | Content | Data / script |
|---|---|---|
| III | parity: odd harmonics -> odd m, even -> even m; E_x ~ x has P(m=0) = P(m=2) | analytic |
| IV A | coherent emission per mode vs radiation integral, 0.3-1.8 % (Fig. 1, Table I) | `benchmarks/harmonics` (`radint.py`, `emission.deck`, `spectra/`), `scripts/fig_emission.py` |
| IV B | helical seeded 10 nm / 0.8 nm vs GENESIS; SASE vs Ming Xie | `benchmarks/genesis_seeded`, `examples/fel_sase.deck` |
| V | seeded planar FEL 10.6 nm: P2 = 9.6 MW (GENESIS: none), P3, bunching (Fig. 2, Table II); dz <= lambda_u/32 for planar | `benchmarks/harmonics` (`fel.deck`, `genesis/`, `fel_tables/`), `scripts/fig_fel.py` |
| VI | 1D averaged vs unaveraged, rho scan (Fig. 3) | `benchmarks/high_rho` (`fel1d.py`, `scan.txt`), `scripts/fig_averaging.py` |

## Open items (my suggested order)
1. **Independent check of P2**: Puffin (unaveraged, radiates even harmonics) on the same case, or the
   Huang-Kim (PRE 62, 7295) harmonic theory. Without it, P2 rests only on the emission validation.
2. **Planar step instability** (dz = lambda_u/16 blows the beam up): find the cause, ideally fix it.
   Candidate: s-discretization of the co-moving field of j_x ~ cos(k_u s). If it is fixed, the cost halves.
3. **3rd harmonic**: rerun with dzeta = lambda_r/48 and m_max = 3, about 3 h and 8 GB. That needs a machine with
   more memory than this container (7 GB). The present 20 % deficit at saturation is probably resolution.
4. **Transverse profile of the 2nd harmonic** (zero on the axis, E_x ~ x) from the data files: a figure panel.
   The HDF5 dumps of the run are not in git (684 MB each). Rerun `fel.deck` with `output.every` at the
   saturation step.
5. **GENESIS with the same physics** (LSC + widening) for the planar case, as in the helical benchmark, to
   explain the 13 % in gain length.
6. Optional, for a stronger paper: superradiant short bunch (time-dependent); the LPA -> transport -> FEL chain
   (planned combined run); OAM harmonics in a helical undulator.
7. Cost comparison with time-dependent GENESIS on the same bunch.
8. Check the references marked `\todo{check}` (Freund 2000, Tremaine 2002), the GENESIS 4 citation, and
   WarpX/boosted-frame FEL references; fill in the arXiv number of the code paper.

## Note
The repository is public: the draft and the benchmark results are visible once pushed.
