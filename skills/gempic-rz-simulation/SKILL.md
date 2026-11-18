---
name: gempic-rz-simulation
description: Write, check, run and debug input decks for the plasma simulation code GEM-PIC-RZ (Galilean electromagnetic particle-in-cell, quasi-cylindrical RZ geometry, window moving at c). Use for laser wakefield acceleration (LWFA), beam-driven wakefields (PWFA), witness bunches, ionization injection, undulators and free-electron lasers (FEL), lasers from LASY files, flying focus, plasma lenses and chicanes, dielectric or metal structures, and for choosing the grid resolution, reading the output (bunches.txt, axis, maps, HDF5, openPMD) or finding out why a run gives wrong or noisy results. Not for full-PIC codes (FBPIC, WarpX, Smilei, EPOCH) or the quasi-static code QUARZ. If gempic_rz is not built yet, use the skill gempic-rz-setup first.
license: BSD-3-Clause
compatibility: Needs the gempic_rz executable (skill gempic-rz-setup). The helper scripts need python3; plotting needs numpy and matplotlib, LASY files need h5py and a build with HDF5.
metadata:
  repository: https://github.com/AlexMPukhov/GEM-PIC-RZ
---

# GEM-PIC-RZ: input decks

GEM-PIC-RZ reads one ASCII deck: blocks `name { key = value ... }`, every dimensional value with a unit. The code
checks the deck strictly (unknown keys, missing or wrong units and conflicts are errors with the line number), so
**let the code check your deck instead of guessing**:

```
gempic_rz my.deck --check        # reads the deck, prints the summary in physical units and all warnings, exits
```

## Is GEM-PIC-RZ the right code?

| the user wants | use |
|---|---|
| laser- or beam-driven wakes in underdense plasma, external or ionization injection, beam transport, undulator radiation and FEL with the radiation on the grid, structured lasers (LASY, flying focus), dielectric or metal structures | **GEM-PIC-RZ** (this skill) |
| long beam-driven stages (metres of plasma, proton or electron drivers) where neither a laser carrier nor radiation has to be resolved | a quasi-static code, e.g. **QUARZ** (github.com/AlexMPukhov/QUARZ): much larger steps. GEM-PIC-RZ can do these cases too, at a higher cost |
| density down-ramp injection or self-injection with a quantitative charge; reflection, backward Raman scattering, colliding pulses; overdense or solid targets | a **full PIC code** (FBPIC, WarpX, Smilei, ...). A beam injected there can be imported here for the acceleration, transport and FEL stages (`bunch.file`, openPMD) |
| geometry far from a few azimuthal modes around one axis | a 3D code |

Say so before writing a deck if the request falls into rows 2 to 4, and propose the split. GEM-PIC-RZ itself is
not quasi-static (rule 8 below).

## Start here

Paths in this skill (`examples/`, `scripts/`, `references/`) are relative to the directory of this `SKILL.md`,
not to the working directory of the user. Copy decks into the working directory and call the scripts with their
full path.

1. `gempic_rz --version` works? If not: skill `gempic-rz-setup`. (The executable may also be `<repo>/build/gempic_rz`
   or wherever the user built it.)
2. **New simulation:** copy the closest deck from `examples/` (table below), change the physics, then follow
   "Verify before delivering". Do not write a deck from memory: keys that do not exist are errors.
   For a plain "LWFA, a0, density, length" request: `examples/lwfa_plain.deck` with the numbers of
   `scripts/plasma_numbers.py`. Add a witness bunch only if the user asks for one.
3. **Which keys exist, their units and defaults:** `references/input_deck.md` (the reference of the code, complete
   and long: search it for the block you need, e.g. `## \`laser\``, instead of reading it through).
4. **How fine must the grid be, how large the box, what the code cannot model:** `references/physics_rules.md`.
   Read it before you choose `grid` values for a new case.
5. **Reading and plotting the output:** `references/diagnostics.md`.
6. **A run fails or looks wrong:** `references/debugging.md`.

Ask the user only for physics that you cannot infer (species, density profile, what they want to measure); choose
the numerical parameters yourself from the rules and say which ones you chose.

## Examples (all run as they are)

| deck | case | run time |
|---|---|---|
| `examples/lwfa_plain.deck` | plain LWFA: a0 = 2, 0.8 um, 30 fs, matched spot, 1e18 cm^-3, 1 mm; every number from `scripts/plasma_numbers.py` | PLAIN_TIME |
| `examples/lwfa_witness.deck` | LWFA, a0 = 2, 3.5e18 cm^-3, external witness bunch, free-hand radial grid (fine near the axis) | a few minutes on 2 cores |
| `examples/pwfa_driver_witness.deck` | beam-driven wake: 1 GeV drive and witness bunches, 1e17 cm^-3, `m_max = 0` | minutes |
| `examples/fel_sase.deck` | high-gain FEL from noise: helical undulator, 4 kA, 10 nm; gain length against Ming Xie | ~8 min on 2 cores |
| `examples/undulator_emission.deck` | coherent emission of a short bunch in a planar undulator against the radiation integral | ~2 min on 2 cores |
| `examples/lasy_laser.deck` | LWFA with the laser from a LASY envelope file, propagated exactly to its focus (`scripts/make_lasy_gaussian.py` writes the file) | LASY_TIME |

## The deck in one page

```deck
constants { n0 = 1e18 cm^-3; w0 = 15 um }     # optional: named values for expressions
grid {                         # once: code units, modes, the three grids
  lambda0 = 0.8 um             # reference wavelength: the laser's (FEL: the resonant wavelength)
  m_max   = 1                  # azimuthal modes 0..m_max: 1 for a linearly polarized laser, 0 for a centred beam driver
  radial_grid = uniform; dr = 0.53 um; r_max = 64 um       # or stretched / table (free-hand)
  dzeta = 32 nm; nzeta = 2503                              # window co-moving at c (80 um); zeta = 0 is its REAR end
  dz = 9.57 um; z_end = 1 mm                               # steps along the propagation distance z
  r_absorber = 53.5 um                                     # radial absorber outside this radius
}
laser  { a0 = 2; waist = w0; tau_fwhm = 30 fs; zeta_center = 56.77 um; z_focus = 100 um }   # zeta_center: from the rear end
plasma { density = n0; z_points = 0, 100 um; z_values = 0, 1; r_max = 51 um }   # n = density * f(z) * g(r)
bunch  { name = witness; charge = 5 pC; energy = 500 MeV; sigma_x = 1 um; sigma_z = 1 um; zeta_center = 30 um }
output { dir = run1; every = 25 }
numerics { iterations = 2 }
```

## Rules that decide whether the result is right

Details, numbers and the reasons are in `references/physics_rules.md`.

1. **Laser carrier resolved:** `grid.dzeta` <= lambda / 20 over the whole pulse (the code warns below 12 points).
   Without a laser: resolve the bunch (sigma_z / 10) and the plasma wavelength.
2. **Window:** the laser centre at least 3 `length_1e` behind the front (`length_1e` = c tau_fwhm / 1.177, the 1/e
   half length of the field; tau_fwhm is the FWHM of the intensity); behind it as many plasma wavelengths as you
   want to see (1.7 for the first bucket). Bunches with gamma < 1e4 in vacuum need the gaps printed in the summary
   line "box for bunch".
3. **`grid.dz` is not a Courant limit** (the step in z is implicit). It must resolve what changes along z: the
   density ramps, the laser focus, betatron oscillations (PWFA: betatron length / 15), an undulator period
   (period / 16). LWFA: start with about 1.8 c/omega_p as in `examples/lwfa_witness.deck`, then halve it and compare.
4. **The radial grid must not reflect radiation.** A cell larger than the transverse wavelength throws the wave
   back, and a grid that coarsens outward acts like a falling refractive index (numerical total internal
   reflection). Keep k_perp dr <~ 1 wherever radiation must propagate freely (laser: uniform dr out to the absorber;
   undulator radiation: its cone, also of the harmonics); coarsen only outside that region, slowly (< 20 % per cell).
   Symptom of a violation: fringes near the axis with a fixed radial wavenumber.
5. **Box:** `plasma.r_max` < `grid.r_absorber` < `grid.r_max`, with `grid.r_absorber` >= 2.5 w, where w is the laser
   spot where it is largest in the run. A matched spot or P/P_c near 1 does not guarantee guiding: over more than
   a Rayleigh length size the box with the vacuum spot w0 sqrt(1 + (z - z_focus)^2 / z_R^2), or use the smaller box
   and check in the result that |E_x| at the absorber stays below ~1 % of the peak (`scripts/plot_run.py` prints it).
6. **Undulators / FEL:** `grid.dzeta` <= lambda_r / 12 (better 16-20), `grid.dz` and `numerics.beam_dz_max` <=
   period / 16, `numerics.beam_filter = 0` in vacuum, `numerics.theta` left at 0.5, the window long enough for
   the slippage N_u lambda_r, the bunch near the front.
7. **Model limits:** forward-moving physics in underdense plasma only. Radiation moving backwards (reflection,
   backward Raman, colliding pulses) and solid targets are outside the model. **Density down-ramp injection and
   self-injected charge are not quantitative**: do that stage with a full PIC code and import the beam
   (`bunch.file`, openPMD). Ionization injection, acceleration, transport and FEL are fine.
8. Never call GEM-PIC-RZ "quasi-static": it keeps the slow derivative of the forward field.

`python3 scripts/plasma_numbers.py --density 1e18 --a0 2 [--waist 15 --tau 30] --length 2` prints lambda_p,
c/omega_p, P/P_c, the matched spot and starting values for the grid, the box, the ramp, the focus,
`numerics.trap_gamma` and `plasma.smooth_length` that follow these rules. Without `--waist` and `--tau` it proposes
the matched spot and c tau_fwhm = lambda_p / 3; tell the user if you take them.

## Verify before delivering

A deck is done only when these hold. Do them in this order.

1. **Check:** `gempic_rz my.deck --check` exits without error. Read the summary: are the numbers (lambda_p,
   Rayleigh length, power, number of slices and steps, bunch energy) what the user meant?
2. **Warnings:** every `WARNING:` line is either fixed or explained to the user. Do not deliver a deck with an
   unexplained warning.
3. **Smoke run:** a few steps in a scratch directory, with the backend of the full run, e.g.
   `OMP_NUM_THREADS=8 gempic_rz my.deck numerics.backend=openmp grid.z_end=50um output.dir=smoke output.every=0`.
   It must end with "peak memory"; at step 0 the console `|Ex|max` equals a0 (a0 / sqrt(1 + z_focus^2 / z_R^2) with
   the focus ahead); `smoke/axis_*.txt` holds finite numbers. Without a bunch, `bunches.txt` has only `trapped` lines
   with zeros: that is right.
4. **Cost:** the line `done: N steps in T s` of the smoke run gives the time per step; times the number of steps
   in the summary it is the run time (cell counts of other decks are a poor guide: the particle number differs).
   Tell the user; ask before starting anything above ~1 hour.
5. **Full run and a look at the result:** `python3 scripts/plot_run.py <output.dir>` prints a table per output
   step (laser amplitude, wake field, field at the absorber, bunches) and, with matplotlib, writes a figure.
   Compare the wake amplitude, the bunch energy or the gain with the estimates in `references/diagnostics.md` and
   report the numbers, not impressions. If the field at the absorber exceeds ~1 % of the peak, enlarge the box.
6. **Convergence**, for any number the user will rely on: one more run with `grid.dz` halved (then `grid.dzeta`,
   `grid.dr`). For a first look it may be skipped; then say that it was not done.

Parallel runs: `OMP_NUM_THREADS=8 gempic_rz my.deck numerics.backend=openmp`, or `mpirun -np P gempic_rz my.deck`
(skill `gempic-rz-setup`). Results do not depend on the backend.
