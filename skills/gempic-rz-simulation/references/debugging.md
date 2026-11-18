# Debugging a GEM-PIC-RZ deck or run

Work from the top: most problems are caught by the first three steps.

## 1. The deck does not pass `--check`

The message has the file, the line and the reason. Typical ones:

| message contains | cause and fix |
|---|---|
| `unknown key` / key never read | the key does not exist in this block, or does not apply to the chosen option (e.g. `grid.dr` with `grid.radial_grid = stretched`). Look the block up in `input_deck.md`; do not invent keys |
| `needs key '…' (found '…': misspelled?)` | a required key is missing or misspelled |
| unit missing / wrong kind | every dimensional value needs a unit (`8 um`, `25 fs`, `1e18 cm^-3`, `5 pC`, `500 MeV`); a list has one unit at its end (`z_points = 0, 50, 5000 um`). Units are case-sensitive (`MeV`) |
| expression errors | quote expressions (`"min(z / (100 um), 1)"`); write every dimensional number inside with its unit, the variables (r, z, zeta, x, y) are in code units |
| conflicting keys | two alternatives given together (`grid.nr` and `grid.r_max`, `laser.tau_fwhm` and `laser.length_1e`, `grid.window_length` and `grid.nzeta`) |
| HDF5 not available | the build has no HDF5: install the HDF5 C library and rebuild (skill gempic-rz-setup), or use `output.format = text` |
| needs `grid.m_max` ≥ 1 | a laser, an offset bunch or `output.lasy` needs mode 1 |

## 2. Warnings of the summary

| warning | what to do |
|---|---|
| pulse front closer than 3 length_1e to the front of the window | move `laser.zeta_center` back or lengthen the window; tolerable (small energy cut) only if deliberate |
| pulse resolved by fewer than ~12 points per wavelength | reduce `grid.dzeta` (λ/20–λ/25) |
| largest cell-to-cell change above 25 % | add reference points to the free-hand grid so the step grows more slowly (< 20 %) |
| the window holds N plasma periods: long wakes need `numerics.iterations = 4` | set it |
| bunch: the radial boundary is closer than the reach of its fields / its space-charge field reaches the ends of the window | bunch with γ < 10⁴ in vacuum: enlarge `grid.r_max` and the gaps to the values of the line "box for bunch"; set `grid.r_absorber` beyond `grid.r_max` |
| undulator: `grid.dz` should be ≤ period/16; resonant wavelength is only N ζ steps | refine as stated |
| `numerics.beam_dz_max` exceeds 0.1/k_β; element shorter than 4 sub-steps | reduce `numerics.beam_dz_max` |
| no plasma block: vacuum propagation | fine if intended |
| output volume exceeds the free space / directory on tmpfs | larger `output.every`, fewer formats, `output.particles = false`, another `output.dir` |

## 3. Smoke run (a few steps) before the long run

`gempic_rz my.deck grid.z_end=<3–10 steps> output.dir=smoke output.every=0`, then look at the console lines and
`smoke/axis_00000.txt`:

- `|Ex|max` at step 0 equals a0 (or a0/√(1 + z_focus²/z_R²) for a focus ahead). If it is much smaller, the laser is
  cut by the window or by the absorber, or the wrong `laser.zeta_center` puts it outside the window.
- The bunch lines show the charge, energy and sizes of the deck.
- The time per step × number of steps is the run time.

## 4. The run finishes but the result looks wrong

| symptom | likely cause | fix |
|---|---|---|
| laser amplitude decays or the pulse lags in vacuum | `grid.dzeta` too coarse; the absorber cuts the beam (`grid.r_absorber` < 2.5 w at the widest point); `grid.r_max` too small | λ/25; enlarge the box |
| fringes near the axis with a fixed radial period; radiation comes back to the axis | the radial grid reflects: cells larger than the transverse wavelength, or a grid that coarsens where radiation propagates | keep k_⊥ dr ≲ 1 in that region; extend the fine region and compare (physics_rules.md) |
| noise stripe on the axis | too few macroparticles in the small axis cells | `plasma.ppc_r_axis`, `plasma.ppc_phi_axis`, `plasma.r_axis_refine` ≈ radius of the fine region |
| spike and ripples of E_z and density behind the bubble | singular sheath closure on the axis | `plasma.smooth_length` ≈ 0.1 c/ω_p; judge convergence by the fields ahead of the spike, not by its height. On-axis n_e/n0 of 10–50 at the closure is normal, also with smoothing |
| laser amplitude falls in the plasma faster than expected; `plot_run.py` shows more than ~1 % of the peak field at the absorber | the pulse is not (fully) guided and its wing reaches the absorber and the plasma edge | enlarge `grid.r_max`, `grid.r_absorber`, `plasma.r_max` (vacuum spot at the end of the run × 2.5) and compare |
| wake amplitude grows from step to step in a long window | passes per slice not converged | `numerics.iterations = 4` |
| "trapped" charge appears although nothing should be injected | electrons quivering in a strong laser move ahead of their level and are counted | `numerics.trap_gamma` somewhat above 1 + a0²/2; `numerics.trap_advance` > 1 |
| self-injected or down-ramp charge disagrees with a full PIC code | known limit of the model | do the injection with a full PIC code, import with `bunch.file` |
| witness emittance or size grows without physical reason | bunch not resolved radially (dr > σ_r/5) or the betatron motion not resolved | fine cells near the axis; `numerics.beam_dz_max` ≤ 0.1/k_β; check with half `grid.dz` |
| bunch starts with a wrong size | `bunch.z_waist` ≠ 0 drifts it back to z = 0; mismatched emittance | compare with the summary line of the bunch |
| energy gain differs by ~1 % between runs with different `grid.dz` | steps in the density up-ramp | smaller steps there (`grid.z_grid = table`) |
| results change with `numerics.seed` behind the bubble | chaotic sheath crossing | not an error; compare ahead of the bubble closure or average |
| FEL: no gain or a wrong gain length | `grid.dzeta` > λ_r/12; `grid.dz` > λ_u/16; beam filter on; `numerics.theta` > 0.5; bunch shorter than a few cooperation lengths; window shorter than the slippage; `grid.lambda0` far from λ_r | see "undulator / FEL" in physics_rules.md |
| FEL: beam blows up in a planar undulator | beam pushed in the fields of one level | leave `numerics.beam_field` at its default (`mean` with undulators) |
| numerical Cherenkov noise around a bunch at the grid scale | beam current at the ζ-Nyquist wavenumber | `numerics.beam_filter = 1` (default; 0 only for FEL in vacuum) |
| radiation of an FEL or a flying focus distorted near the outer boundary | absorber acting on the wanted field | flying focus and `laser.inject = boundary`: `grid.absorber_strength = 1` and no foreign plasma in the layer |

## 5. The run stops or is too slow

- **Stops with a message:** read it; wall and dielectric geometry are checked at every step (the wall must stay
  inside the grid, a dielectric inside the absorber radius).
- **Killed / out of memory:** memory is about 0.25 MB per ζ slice (150 radial cells, m ≤ 1) plus 0.1 kB per
  macroparticle, per MPI rank only its block. Output on tmpfs counts as memory.
- **Disk full:** the output, not the run; see the line `output volume`.
- **Too slow:** count slices × radial cells × steps. Use `numerics.backend=openmp` with `OMP_NUM_THREADS`, or MPI for
  long windows; `grid.m_max = 0` without a laser; free-hand ζ and z grids (fine only where needed);
  `numerics.iterations = 1` in vacuum; larger `output.every`.
- **MPI:** at least 4 slices per rank; the efficiency is N_steps/(N_steps + P − 1), so few-step runs gain nothing.
  `numerics.beam_consistent` is not available with several ranks.
- **Restart refuses:** the grid (modes, radial and ζ grid, `grid.lambda0`) and the ion blocks must be those of the
  checkpoint; `grid.z_end`, `grid.dz`, output settings and the backend may change.

## 6. Still wrong

Reduce the case: switch off one block at a time (bunch, ions, undulator), run in vacuum, compare `serial` with
`openmp` (they must agree), halve the steps. If a small deck still shows a behaviour that contradicts the
documentation, report it at https://github.com/AlexMPukhov/GEM-PIC-RZ/issues with the deck, the version
(`gempic_rz --version`) and the console output; do not patch the source to make a deck run.
