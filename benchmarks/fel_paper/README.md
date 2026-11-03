# FEL paper: input decks

Inputs of the runs in A. Pukhov, *First-principles simulation of free-electron lasers with the Galilean
electromagnetic particle-in-cell method: even harmonics of planar undulators* (arXiv:2611.01640, https://arxiv.org/abs/2611.01640). Run with `gempic_rz <deck> [key=value ...]`; the overrides on the
command line change single deck entries. Many threads need `numerics.backend = openmp` (set in the decks here).

| Paper | Deck | Notes |
|---|---|---|
| Fig. 1, Table I (coherent spontaneous emission, modes m = 0-3) | `../harmonics/emission.deck` | reference: `../harmonics/radint.py`, comparison `../harmonics/compare_emission.py` |
| Helical seeded amplifiers at 10 nm and 0.8 nm | `../genesis_seeded/gem.deck`, `gem_x.deck` | GENESIS 4 inputs in `../genesis_seeded/genesis*/`, comparison `compare.py` |
| Helical SASE to saturation (focused, shot noise) | `sase_focus.deck` | random loading: override line in the deck header |
| Table II, Fig. 2 (seeded planar FEL, reference at lambda_u/32) | `fel_clean.deck` + `grid.dz=0.3125mm numerics.beam_dz_max=0.3125mm` | 9 h on 64 cores, 44 GB |
| Fig. 3 (transverse structure at z = 0.9 m) | `fel_clean.deck` + `output.format=both output.every=1440` | lambda_u/16 |
| Step size lambda_u/16 vs /32 | `fel_clean.deck` | |
| GENESIS 4 steady state with h = 1, 2, 3 (Table II, Fig. 2) | `../harmonics/genesis/fel.in`, `fel.lat` | |
| Runs on the stretched radial grid (cost, quiet start, 3rd harmonic with M = 3) | `../harmonics/fel.deck` | see below |
| GENESIS 4 time-dependent run (cost) | `genesis_td/td.in` (ngrid 401), `td201.in` (201) | lattice `../harmonics/genesis/fel.lat` |
| Fig. 4 (period averaging in 1D) | `fel1d/scan.py` (uses `fel1d.py`) | numpy, scipy; writes `scan.txt` |

Runs on the stretched grid (`../harmonics/fel.deck`, 122 radial nodes, offset orbit), all with
`numerics.backend=openmp`:
- cost (58 min on two cores without, 91 min with the refinement near the axis): the deck as it is, with and without
  `bunch.axis_refine=1`;
- quiet start: `grid.dz=0.3125mm numerics.beam_dz_max=0.3125mm bunch.shot_noise=quiet` (and without the last entry
  for the comparison with shot noise);
- third harmonic with M = 3: `grid.m_max=3 grid.dzeta=0.22135nm` (lambda_r/48), at lambda_u/16 and /32.

Why the reference deck differs from `../harmonics/fel.deck` (centred orbit, uniform radial grid): see the header of
`fel_clean.deck` and the paper (Sec. V). A radial grid must not stretch where radiation with k_perp dr >~ 1 still has
to travel outwards.

The Puffin comparison runs are not included.
