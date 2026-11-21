# Choosing the parameters of a GEM-PIC-RZ run

Keys are written `block.key`; their units and defaults are in `input_deck.md` (same directory).

## The model in five lines

- Quasi-cylindrical geometry: fields are azimuthal modes m = 0 … `grid.m_max` on an (r, ζ) grid. Mode 0 holds the
  wake and axisymmetric beams, mode 1 a linearly or circularly polarized laser and the dipole motion of offset beams.
- The window moves with c: ζ = z − ct. Inside it the fields are marched from the front to the rear, slice by slice
  (`grid.dzeta`); the whole window then advances along the propagation distance by `grid.dz`.
- Plasma particles stream through the window from the front to the rear once per step. Particles that move forward
  with the window (bunches, trapped electrons) are kept in a beam list and pushed along z.
- The laser carrier is resolved (no envelope model); the field equations keep the slow derivative of the forward
  wave, so diffraction, depletion, red-shift and self-focusing are described. It is **not** a quasi-static code.
- The slow derivative of the backward wave is dropped: everything that propagates backwards in the lab is outside
  the model.

## What it can and cannot model

| works (tested, see the repository's README and `benchmarks/`) | outside the model or not quantitative |
|---|---|
| LWFA in underdense plasma, linear to bubble regime; laser evolution over many Rayleigh lengths (against FBPIC) | reflection, backward Raman scattering, counter-propagating or colliding pulses |
| beam-driven wakes (electrons, positrons, protons), drive + witness, hosing and offsets with `grid.m_max` ≥ 1 | overdense plasma, solid targets, anything not moving roughly with the window |
| external injection of a witness; beam loading; emittance and energy spread | **density down-ramp injection and self-injection: the trapped charge is too small by a factor 2–5 and depends on the step** → use a full PIC code (FBPIC) for that stage and import the beam with `bunch.file` |
| ionization injection (`ion.ionization = adk`), against FBPIC to 10 % in charge | energy taken from the laser by ionization; collisional ionization |
| mobile ions, ion collapse | strongly three-dimensional structures that need many azimuthal modes |
| undulators and FEL with the radiation on the grid (against GENESIS), harmonics, space charge | wavelengths the ζ grid does not resolve |
| plasma lenses, quadrupoles, solenoids, thin chicanes; chained stages through openPMD files | synchrotron radiation and dispersion inside a chicane |
| lasers from LASY files, flying focus fed through the radial boundary | |
| metal walls of slowly varying radius, dielectric linings, a resonant medium (`bound_electrons`) | reflections at steps and irises; dispersive dielectrics |

If the user asks for something in the right column, say so before writing a deck, and propose the split (e.g.
injection with FBPIC, then acceleration here; worked example in the repository: `examples/chain_lwfa_fel/`).

## Step by step: a laser-driven wake

1. `python3 scripts/plasma_numbers.py --density <cm^-3> --a0 <a0> [--waist <um> --tau <fs>] --length <mm>` gives
   λ_p, c/ω_p, P/P_c, the matched spot and starting values for everything below. If the user gave no spot or
   duration, it proposes the matched spot w0 = 2√a0 c/ω_p and c τ_FWHM = λ_p/3 (τ_FWHM: FWHM of the intensity, the
   `laser.tau_fwhm` of the deck; the pulse then fits into the first half of the plasma wave). Say that you chose them.
2. `grid.lambda0` = the laser wavelength; `grid.m_max = 1`.
3. **ζ grid.** `grid.dzeta` = λ/25 (20 is the minimum for quantitative work; the code warns below 12).
   Window = distance from the rear end to the laser + 3 `length_1e` ahead of the laser centre
   (`length_1e` = c τ_FWHM / 1.177, the 1/e half length of the field). Behind the laser about 1.7 λ_p shows the
   first bucket; add λ_p per further bucket. Give the window as `grid.nzeta` (a whole number of slices); a window
   that is not a multiple of `grid.dzeta`, or a cut pulse front, gives a warning.
   With `grid.zeta_grid = table` the slices may be coarser where there is no laser (behind the bubble), which saves
   time; see `input_deck.md`, "Free-hand grids".
4. **Radial grid.** In the region the laser occupies use a uniform step, `grid.dr` ≤ min(w0/20, 0.1 c/ω_p). A thin
   witness or trapped bunch needs σ_r/10 near the axis: then use `grid.radial_grid = table` with the fine step out
   to a few σ_r and a slow growth (< 20 % per cell, printed in the summary) to the outer step, as
   `examples/lwfa_witness.deck` does. Box: `plasma.r_max` < `grid.r_absorber` < `grid.r_max`, absorber at ≥ 2.5 w,
   where w is the largest spot during the run. In vacuum w = w0 √(1 + (z − z_focus)²/z_R²). In plasma the pulse is
   guided only in part, even with the matched spot and P ≈ P_c (the head of the pulse always diffracts): over
   more than one Rayleigh length either size the box with the vacuum formula (safe, more cells) or keep 2.5 w0
   and check afterwards that |E_x| at the absorber stays below ~1 % of the peak (`scripts/plot_run.py` prints
   this for every output step). A laser wing that reaches the absorber and the plasma edge is cut, which
   changes the evolution of the pulse from then on.
5. **z steps.** `grid.dz` has no stability limit tied to `grid.dzeta`. It is an accuracy parameter: the examples use
   about 1.8 c/ω_p (5 µm at 3.5·10¹⁸ cm⁻³), which gives the energy gain to about 1 %; the steps in a density
   up-ramp matter most (`grid.z_grid = table` allows small steps there only). Halve it once and compare.
6. **Plasma profile.** Always start with a ramp, not a density step at z = 0: `plasma.z_points`, `plasma.z_values`
   (or `plasma.density_function`), about 3 λ_p long if the user gave none; the profile is constant beyond the
   last point. Put the focus (`laser.z_focus`) at the end of the ramp unless the user says otherwise.
   A measured or simulated profile (gas jet from interferometry, nozzle or hydro simulation) is read from a file
   with `plasma.density_file`: a text table whose header states the geometry and the columns (`# geometry = z`,
   `# columns = z[mm] n[cm^-3]`; or `# geometry = rz`, `# columns = z[mm] r[um] n[cm^-3]`), or an openPMD mesh. A
   file with n(z) may also carry the transverse shape (`# transverse = "exp(-(r / (50 um))^2)"`). A density with a unit in the file is the absolute density (`plasma.density` then only
   normalizes the output); without unit it is a factor on `plasma.density`. Do not retype such data into
   `plasma.z_points`. `plasma.file_z` shifts the file along z; outside the file the density is zero unless
   `plasma.file_outside = hold`. Read the summary of `--check`: it prints the range of the file against the run
   and the peak density, and warns if the file is cut at a finite density or has gradients shorter than 4 steps
   (then reduce `grid.dz` there, or smooth the data). `ion.density_file` does the same for a species with
   `ion.profile = own`.
7. **Macroparticles.** The defaults (`plasma.ppc_r = 2`, `plasma.ppc_phi = 8`, refined near the axis) are those of
   the benchmarks. With a fine axis region set `plasma.r_axis_refine` to about its radius.
8. **Bubble regime (a0 ≳ 2):** the density spike where the sheath closes on the axis is singular. For clean
   fields behind the bubble use `plasma.smooth_length` ≈ 0.1 c/ω_p (not smaller than the axis cell). The on-axis
   n_e/n0 still shows spikes of tens there; judge by the fields ahead of the spike. To count only accelerated
   electrons as trapped set `numerics.trap_gamma` to about twice 1 + a0²/2 (6 for a0 = 2; use the largest a0 of
   the run if the pulse self-focuses).
9. **Long windows:** more than ~30 plasma periods in the window need `numerics.iterations = 4` (the code warns).

## Step by step: a beam-driven wake

- `grid.m_max = 0` for centred, round beams (4 particles per ring are enough: `plasma.ppc_phi = 4`). Offsets and
  hosing need `grid.m_max = 1`; then load centred bunches with `bunch.azimuthal_copies = 4` to keep the noise out of
  mode 1.
- `grid.lambda0` only defines the code units; 0.8 um is fine.
- `grid.dzeta` ≤ σ_z/10 of the shortest bunch and ≤ λ_p/200; `grid.dr` ≤ σ_r/10 at the bunch (free-hand grid: fine
  near the axis, ~0.1–0.4 c/ω_p in the blowout region).
- `grid.dz` ≤ 1/15 of the shortest betatron length λ_β = 2π √(2γ) c/ω_p of the bunches; the beam itself is pushed
  with sub-steps `numerics.beam_dz_max`.
- Window: driver near the front with ≥ 4 σ_z ahead of it; the witness ~λ_p/2 … λ_p behind it.
- Every bunch starts with its vacuum self-field; start the plasma with a short ramp. For γ < 10⁴ bunches in vacuum
  sections obey the summary line "box for bunch".

## Step by step: undulator / FEL

- `grid.lambda0` = the resonant wavelength λ_r = λ_u (1 + K²/2)/(2γ²) (helical: 1 + K²); the summary prints it.
- `grid.m_max = 1` (the radiation of a centred beam is in mode 1; even harmonics of a planar undulator need 2).
- `grid.dzeta` ≤ λ_r/12, better λ_r/16–20. `grid.dz` ≤ λ_u/16 and `numerics.beam_dz_max` ≤ λ_u/16.
- `numerics.beam_filter = 0` in vacuum (the filter damps the current at λ_r). Do not change `numerics.theta`.
- `numerics.iterations = 1` without plasma.
- Window: bunch near the front; behind and ahead of it the gaps of "box for bunch"; length ≥ bunch + slippage
  N_u λ_r + gaps. The bunch must be several cooperation lengths long for SASE.
- Radial box: the space-charge field needs `grid.r_max` ≥ γ L (bunch length L) and no absorber
  (`grid.r_absorber` > `grid.r_max`); the radiation needs k_⊥ dr ≲ 1 inside its cone (next section).
- Shot noise: `bunch.shot_noise = fawley` with `bunch.noise_wavelength` = λ_r for the noise of the real electrons;
  `quiet` for a seeded amplifier. Diagnostics: `fel.txt`, `fel_power.txt` (written automatically).

## The radial grid must not reflect radiation

This is the rule that is easiest to violate with a free-hand radial grid (`grid.radial_grid = table`).

1. **Cutoff.** A cell larger than the transverse wavelength cannot carry the wave (discrete cutoff
   k_⊥ = 2/dr): the wave is thrown back.
2. **Numerical total internal reflection.** Below the cutoff the numerical dispersion is
   k_⊥,eff = (2/dr) sin(k_⊥ dr/2). It falls as dr grows outward and acts like a refractive index that decreases
   with r: waves at grazing incidence are reflected back to the axis.

How to set the grid:
- Estimate the largest k_⊥ = k θ that must leave the beam region: the diffraction of the smallest feature
  (k_⊥ ≈ 2/w for a spot w), the radiation cones of all relevant harmonics, and off-axis harmonics that land at the
  diagnosed frequency. For a planar undulator the harmonic h′ is seen at the frequency of h under the angle
  γθ = √((1 + K²/2)(h′/h − 1)). A flying focus with cone angle θ needs dr ≤ λ/(10 sin θ) everywhere.
- Keep k_⊥ dr ≲ 1 wherever that radiation must propagate freely.
- Let dr grow only outside that region, or behind an absorber that removes the radiation first, and stretch
  slowly (the summary prints the largest cell-to-cell change; keep it below 20 %).
- **Symptom of a violation:** fringes near the axis with a fixed radial wavenumber that disappear when the fine
  region is extended. Test: extend the fine region and compare.

## Output volume and cost

- Cost ∝ slices × radial cells × steps × (particles per cell). Use the run times in the headers of the example
  decks to scale.
- Memory does not grow with the steps; the output files do. `output.every` field files of a large window can fill
  a disk: read the line `output volume: …` printed at the first output step. `bunches.txt` holds the beam moments
  of every step anyway, so `output.every` can be large. For big runs: `output.format = hdf5`,
  `output.hdf5_precision = single`, `output.hdf5_compression = 6`, `output.particles = false`.
- Long runs: `output.checkpoint_every` and `gempic_rz deck --restart latest`.

## Convergence check (before trusting a number)

Change one parameter at a time and compare the quantity of interest (bunch energy and spread in `bunches.txt`, the
on-axis E_z, the FEL gain length): halve `grid.dz`; reduce `grid.dzeta` by 1.5–2; halve `grid.dr` near the axis;
extend the fine radial region; another `numerics.seed` shows the statistical scatter (the region behind the bubble
closure is chaotic: differences there between two seeds are not a resolution error).
