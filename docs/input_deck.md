# GEM-PIC-RZ input deck reference

The input deck is a plain ASCII file made of blocks:

```
# comment ("!" works too), anywhere on a line
grid {
  lambda0  = 0.8 um
  dr_min   = 40 nm              # a unit is required for every dimensional value
  z_points = 0, 50, 5000 um     # lists: comma- or blank-separated, one unit at the end
}
laser { a0 = 2; waist = 8 um; tau_fwhm = 25 fs; zeta_center = 30 um }   # ";" separates statements
```

| block      | count   | content                                          |
|------------|---------|--------------------------------------------------|
| `constants`| 0 or 1  | named values for the expressions of all blocks   |
| `grid`     | 1       | reference wavelength, modes, radial grid, window, steps, absorber |
| `laser`    | 0 … N   | one block per laser pulse                        |
| `plasma`   | 0 or 1  | background plasma: density and profiles, macroparticles |
| `ion`      | 0 … N   | one block per mobile ion species of the plasma   |
| `bunch`    | 0 … N   | one block per particle bunch                     |
| `undulator`| 0 … N   | static undulator fields acting on the particles  |
| `quadrupole`, `plasma_lens`, `solenoid` | 0 … N | transport elements: static fields acting on the particles |
| `chicane`  | 0 … N   | thin magnetic chicane (R56, delay) applied to the beam |
| `output`   | 0 or 1  | diagnostics                                      |
| `numerics` | 0 or 1  | algorithm parameters                             |

## General rules

- **Names are case-insensitive:** block names and keys can be written in any case. Units are case-sensitive, so `MeV` is not `meV`.
- **Units are required:** every dimensional value needs a unit. A number without a unit, or with a unit of the wrong kind (e.g. `waist = 8 fs`), is an error.
- **Unknown keys are errors:** every key that the code does not read is reported with its line number, so typos cannot pass silently. For a missing required key, the message names a similar key found in the block (`needs key 'zeta_center' (found 'zeta_centre' in line 15: misspelled?)`).
- **Repeated blocks are numbered:** blocks of the same type are numbered in order of appearance, `laser[0]`, `laser[1]`, …
- **Strings:** quote a string that contains blanks, e.g. `dir = "run 1"`.
- **Unit suffix:** a unit may follow the last number directly, e.g. `5nm`.
- **Expressions:** a single value may be an expression, e.g. `waist = "2 * w0"` (see [Functions and expressions](#functions-and-expressions)).

### Units

| kind | accepted units |
|------|----------------|
| length | `nm`, `um` (`µm`, `micron`), `mm`, `cm`, `m`; normalized `k0^-1` (= `1/k0` = `c/w0`), `lambda0` |
| emittance (normalized) | any length unit; also `mm*mrad` (= um) |
| time | `as`, `fs`, `ps`, `ns`; normalized `1/w0`, `T0` (laser periods) |
| density | `cm^-3` (`cm-3`, `1/cm3`, `/cm3`), `m^-3`; normalized `n_c` |
| charge | `fC`, `pC`, `nC`, `uC`, `C`, `e` (number of elementary charges) |
| energy | `eV`, `keV`, `MeV`, `GeV`, `TeV` (kinetic energy) |
| angle | `rad`, `mrad`, `urad`, `deg` |
| mass | `u` (= `amu`), `m_p`, `m_e`, `kg` |
| magnetic field | `T`, `mT`, `kG`, `G` |
| field gradient | `T/m`, `T/mm`, `kT/m` |
| dimensionless | no unit, or `%` |

**Code units.** The code uses normalized units derived from `grid.lambda0`:

| quantity | code unit |
|----------|-----------|
| length | 1/k0 = λ0/2π |
| time | 1/ω0 |
| field | m c ω0/e (so the laser field amplitude equals a0) |
| density | n_c = 1.1148·10²¹ cm⁻³ / λ0[µm]² |
| charge | e n_c k0⁻³ |

## Functions and expressions

Profiles of plasmas, bunches, lasers and undulators can be given as mathematical functions of the position, as in
HiPACE++ and QUARZ. A function is a quoted expression:

```
constants {                        # named values, usable in every expression below
  n0   = 1e18 cm^-3
  w0   = 8 um
  Lup  = 200 um
}
plasma {
  density          = n0
  density_function = "min(z / Lup, 1) * (1 + 0.2 * (r / (20 um))^2)"
}
laser {
  a0 = 2
  profile_function  = "exp(-(r / w0)^4)"                                # super-Gaussian, at z = 0
  envelope_function = "sin(pi * (zeta - 20 um) / (30 um))^2 * rect(zeta - 35 um, 30 um)"
  zeta_center       = 35 um                                           # reference of the carrier phase
}
bunch {
  charge = 20 pC; sigma_x = 1 um; macroparticles = 100000
  distribution_z = function
  profile_z      = "zeta < 12 um ? (zeta - 10 um) / (2 um) : 1"     # a ramped current
  zeta_min = 10 um;  zeta_max = 14 um
  gamma_function = "400 * (1 + 0.01 * (zeta - 12 um) / (1 um))"      # energy chirp
}
undulator { period = 2 cm; periods = 100; K_function = "1.2 * (1 - 0.05 * max(z - 1 m, 0) / (1 m))" }
```

**Syntax.**

| element | forms |
|---------|-------|
| numbers | `2`, `1.5e-3`, and numbers with a unit: `10 um`, `10um`, `3e18 cm^-3`, `5 fs`, `0.2 k0^-1`, `50%`. Every unit of the deck is accepted and converted to code units |
| variables | depend on the key (table below), in code units (lengths in 1/k0) |
| constants | `pi`, and the entries of the block `constants` |
| operators | `+ - * /`, `^` (also `**`, right-associative, `-x^2 = -(x^2)`), comparisons `< <= > >= == !=` (give 1 or 0), `&& || !`, `cond ? a : b` |
| functions | `sin cos tan asin acos atan sinh cosh tanh exp log` (= `ln`) `log10 sqrt abs floor ceil sign`, `step(x)` (1 for x ≥ 0), `atan2(y, x) min(a, b) max(a, b) pow(a, b)`, `if(cond, a, b)`, `gauss(x, s)` = exp(−x²/2s²), `rect(x, w)` = 1 for \|x\| ≤ w/2 |

- **Units in expressions:** write every dimensional number with its unit, e.g. `r / (20 um)`. The variables are in
  code units, so `r / 20` would mean r/(20 k0⁻¹). The result of an expression is taken in code units as it is.
- **Quotes:** expressions must be quoted when they contain blanks, commas, `=`, `!` or `#` (the deck syntax uses
  these characters). Names are case-insensitive.
- **Errors:** an unknown name, an unknown function or a syntax error is reported with the line of the deck and
  the position in the expression.
- **Command line:** expressions can be overridden like any value, e.g. `'plasma.density_function="exp(-z/(1 mm))"'`.

**Keys that take functions.**

| block | key | variables | meaning |
|-------|-----|-----------|---------|
| `constants` | any name | — | a value (number with unit, or expression of the constants above it), code units |
| any | any single value | — | an expression of the constants, e.g. `waist = w0` or `length_1e = "L / 4"` |
| `plasma` | `density_function` | x, y, z, r | factor of the density: n = density · f(z) · g(r) · density_function. It may break the axial symmetry (x, y): the weights of the macroparticles then vary along φ. Values < 0 count as 0 |
| `laser` | `profile_function` | r | transverse field profile at z = 0, replaces `waist` and the Gaussian beam (no `z_focus`: the field is given where the simulation starts) |
| `laser` | `envelope_function` | zeta | longitudinal envelope of the field, replaces `tau_fwhm` / `length_1e`. `zeta_center` is then optional: the reference of the carrier phase (default: the centre of the intensity) |
| `bunch` | `profile_z` (with `distribution_z = function`) | zeta | current profile on [`zeta_min`, `zeta_max`], any normalization, ≥ 0. Sampled by the inverse of its cumulative distribution |
| `bunch` | `profile_r` (with `distribution_r = function`), `ux_function`, `uy_function` | x, y, r, zeta | transverse density, sampled by rejection in r < `r_extent`; mean transverse momenta [m c] (the energy is kept) |
| `bunch` | `gamma_function` | zeta | γ(ζ) (energy chirp), replaces `energy` / `gamma`; `energy_spread` is relative to the local value. The bunch's γ (self-field, messages) is the mean over the profile |
| `undulator` | `K_function` | z | K(z) (taper), replaces `K`; z is the lab position. The ramps multiply it. Tabulated with 32 points per period |

The laser field is a0 ω · profile_function(r) · envelope_function(ζ) · cos(ω(ζ − ζ_c) + cep) (linear) and the
peak of the product should be 1 for the amplitude a0. The summary prints every function and the estimates it
uses for the checks (w0 from the 1/e radius of the profile, τ from the rms length of the envelope).

## `grid` (required, once)

| key | kind | default | meaning |
|-----|------|---------|---------|
| `lambda0` | length (physical) | wavelength of `laser[0]` | reference wavelength, defines the code units |
| `m_max` | integer | 1 | azimuthal modes 0…m_max. A laser needs ≥ 1; axisymmetric beam drivers need only 0 |
| `radial_grid` | `table` \| `stretched` \| `uniform` | `stretched` | type of radial grid |
| `r_points` | list of lengths | — | table: reference radii, ≥ 0, strictly increasing |
| `dr_points` | list of lengths | — | table: radial step at each reference radius (as many values as `r_points`) |
| `dr_min`, `dr_max` | length | — | stretched: radial step near the axis and far out |
| `r_transition` *or* `transition_cell` | length / number | — | stretched: radius (or cell index ξ0) where the step changes |
| `transition_cells` | number | 8 | stretched: width of the transition, in cells |
| `dr` | length | — | uniform: radial step |
| `nr` *or* `r_max` | integer / length | — | number of radial cells, or the outer radius (rounded up to a whole cell) |
| `zeta_grid` | `uniform` \| `table` | `uniform` | type of the longitudinal grid of the window |
| `dzeta` | length | — | uniform: longitudinal step in the window. Use at least ~20 points per laser wavelength |
| `zeta_points` | list of lengths | — | table: reference positions in the window (ζ = 0 is the rear end), ≥ 0, strictly increasing |
| `dzeta_points` | list of lengths | — | table: slice spacing at each reference position |
| `window_length` *or* `nzeta` | length / integer | — | window length, or number of slices. ζ = 0 is the rear end of the window. Table: the length is rounded up to a whole cell |
| `z_grid` | `uniform` \| `table` | `uniform` | type of the grid along the propagation distance |
| `dz` | length | — | uniform: step along the propagation distance (Δ in the paper) |
| `z_points` | list of lengths | — | table: reference positions along z (z = 0: start of the simulation) |
| `dz_points` | list of lengths | — | table: step at each reference position |
| `z_end` | length | — | total propagation distance |
| `r_absorber` | length | 0.85 r_max | the radial absorber acts on R for r > r_absorber |
| `absorber_strength` | number | 0.2 | damping per step at r_max |

The stretched grid is `r = R(ξ)` with `dR/dξ = dr_min + (dr_max − dr_min)(1 + tanh((ξ − ξ0)/w))/2`, where ξ is the cell index and w = `transition_cells`.

### Free-hand grids (`radial_grid`, `zeta_grid`, `z_grid` = `table`)

The user gives the radial step at reference radii; between them the step is
interpolated linearly in r, before the first point and beyond the last one it
is constant:

```
radial_grid = table
r_points    = 0,  1.5, 4,   24 um      # one unit applies to the whole list
dr_points   = 40, 40,  400, 400 nm
r_max       = 24 um                    # or: nr
```

- **Mapping:** the node positions follow from `dr/dξ = h(r)` with ξ the cell index, so the local cell size is h(r). On each segment `h = h_k + b (r − r_k)` and the mapping is analytic: `ξ = ξ_k + ln(h/h_k)/b`, `r = r_k + h_k (exp(b(ξ − ξ_k)) − 1)/b`. The grid therefore grows geometrically by the factor e^b per cell inside a segment.
- **Reference points:** if the first point is not 0, the step `dr_points[0]` is used from the axis on. `r_points` need not coincide with cell faces; `r_max` (or `nr`) sets the outer edge independently.
- **Smoothness:** h(r) is continuous, its slope jumps at the reference points. The scheme stays second-order accurate (checked with the Laplacian eigenvalues of modes 0–2 under refinement); only the cell-to-cell change |b| = |Δh/Δr| matters. The summary prints the largest change, and a warning is issued above 25 %. Keep it below ~20 %.
- **Resolution:** a uniform grid is `r_points = 0` with a single `dr_points` value.
- **Summary:** `--check` prints Nr, the step range, the actual `r_max`, the largest cell-to-cell change and the reference points in µm.

The window and the propagation distance use the same mapping:

```
zeta_grid     = table              # slices: fine at the witness and at the laser
zeta_points   = 0,   8,   12,    20,    22.5, 25,    38.2 um
dzeta_points  = 127, 127, 31.83, 31.83, 63.7, 31.83, 31.83 nm
window_length = 38.2 um            # or: nzeta
z_grid        = table              # steps: small in the density up-ramp
z_points      = 0,   51,  150 um
dz_points     = 2.5, 2.5, 10.2 um
z_end         = 0.5 mm
```

- **Window:** slice n lies at ζ_n = X(n), and the half-cell n has the length h_n = ζ_{n+1} − ζ_n. Laser pulses need about 20 points per wavelength over their whole length. Coarser steps are fine where the fields are smooth, e.g. behind the bubble or in a long linear wake. A warning is issued if a pulse is resolved by fewer than ~12 points per wavelength.
- **Steps along z:** the level s^1 follows from s^0 by Δ = X(ξ(s^0) + 1) − s^0, so the levels are the nodes of the mapping. A run ends at the first level ≥ `z_end`. Small steps help where the plasma or the laser changes quickly (ramps, focus).
- **Restart:** `dz` or the z table may be changed at a restart; the ζ grid must be the one of the checkpoint.
- **Summary:** `--check` prints Nz, the number of steps, the step ranges, the largest change from one step to the next, and the reference points.


## `laser` (0 … N blocks)

| key | kind | default | meaning |
|-----|------|---------|---------|
| `a0` | number | — | normalized peak amplitude (per component for circular polarization) |
| `wavelength` | length | `lambda0` | laser wavelength (another harmonic or colour is allowed) |
| `waist` | length | — | w0: 1/e field radius at focus (or `profile_function`, see [Functions](#functions-and-expressions)) |
| `tau_fwhm` *or* `length_1e` *or* `envelope_function` | time / length / function | — | FWHM of the intensity, or 1/e half-length of the field envelope (L = τ_FWHM·c/√(2 ln 2)) |
| `zeta_center` | length | — | envelope centre, measured from the rear end of the window |
| `z_focus` | length | 0 | position of the focal plane; z = 0 is the start of the simulation. Before focus the pulse is initialized with wavefront curvature and Gouy phase |
| `cep` | angle | 0 | carrier-envelope phase |
| `polarization` | `linear` \| `circular` | `linear` | |
| `pol_angle` | angle | 0 | linear polarization: angle of E to the x axis |
| `helicity` | +1 \| −1 | +1 | circular polarization: sense of rotation |

The field at z = 0 is the paraxial Gaussian beam
`E = a0 ω/q · exp(−r²/(w0² q)) · exp(−(ζ−ζc)²/L²) · exp(i(ω(ζ−ζc) + cep))`, with `q = 1 + i(0 − z_focus)/zR` and `zR = ω w0²/2`. The pulses of all laser blocks are added.

### Lasers from LASY / openPMD envelope files (`file`)

`laser { file = laser_00000.h5; zeta_center = 30 um }` takes the pulse from an envelope file written by
[LASY](https://github.com/LASY-org/lasy) (or by any code that follows its layout, including GEM-PIC-RZ's own
`output.lasy`): the openPMD mesh `laserEnvelope`, the complex envelope of the transverse electric field at a plane
z = const versus time, E_x,y = Re[p_x,y env e^{−iω_L t}], in `thetaMode` (axes t, r; any number of azimuthal
modes) or `cartesian` (axes t, y, x) geometry, with the attributes `angularFrequency`, `polarization` and
`envelopeField` (`electric_field` or `normalized_vector_potential`). Any spatio-temporal structure is kept: chirps,
pulse-front tilt in mode 1, aberrations, axiparabola or flying-focus fields, Laguerre–Gauss modes, … .

Why it fits: a level s of the code is the plane z = s and its slice ζ the time t = s − ζ at which the field passes
that plane (GT2). An envelope file at a plane is therefore a level of the code; no conversion from time to space and
no assumption about the propagation are needed. The time t_center of the file is put at `zeta_center`, i.e. the slice ζ
holds the file's time τ = t_center + zeta_center − ζ. The carrier keeps its absolute phase e^{−iω_L τ} (the CEP of the file).

- **At the file's plane** (`file_z` = 0, the default): the envelope is interpolated linearly in r and t to the nodes,
  and R = 2E⊥, L = 0 (forward wave) is set in all modes m = 0 … `m_max`. A cartesian file is first resampled on
  rings around the axis (`n_theta` angles) and decomposed into azimuthal modes.
- **At another plane** (`file_z` ≠ 0): the field at z = 0 is computed by the exact propagation of every plane wave,
  per azimuthal mode a Hankel transform and in time an FFT, k_z = (ω² − ω_p² − k⊥²)^{1/2}, without paraxial
  approximation. E.g. a file at the focus with `file_z = 1 mm` starts the pulse 1 mm before its focus. The forward
  transform integrates the interpolant exactly (no quadrature error at the axis).
- **Through the boundary** (`inject = boundary`): as for the `flying_focus` block, the field is also fed through
  the radial absorber layer at every step (total field / scattered field): in r_absorber < r < r_max, R relaxes to the
  incident field of the level, computed with the same propagator (in vacuum or in `background_density`). This brings
  converging (axicon, axiparabola, flying-focus) waves in from the side over any distance. Use
  `grid.absorber_strength = 1`, and a file that reaches beyond r_max (the incident field is zero outside the file).

| key | kind | default | meaning |
|-----|------|---------|---------|
| `file` | string | — | the envelope file (HDF5; relative to the working directory or to the deck) |
| `file_mesh` | string | `laserEnvelope` | name of the mesh record |
| `file_iteration` | integer | the last | iteration of the file |
| `zeta_center` | length | — | window position of the file's time `t_center` |
| `t_center` | time | intensity centroid | time of the file put at `zeta_center` (in the file's time base, which includes LASY's position/c) |
| `file_z` | length | 0 | plane of the file in the coordinates of the simulation (z = 0: its start) |
| `inject` | `box` \| `boundary` | `box` | `box`: initial field only; `boundary`: also fed through the absorber layer at every step |
| `background_density` | density | 0 | uniform plasma of the propagation (`file_z` ≠ 0 and `inject = boundary`); the layer must hold this plasma |
| `a0` | number | from the file | rescales the field so that the peak of \|env\| is a0 ω_L (for circular polarization each component then has a0/√2) |
| `n_theta` | integer | max(32, 4 m_max + 8) | cartesian files: angles of the resampling on rings |
| `time_padding` | integer | 2 | length of the FFT window / length of the file's time axis (the propagated pulse must stay inside it) |

The keys of the built-in pulse (`waist`, `tau_fwhm`, `polarization`, …) are given by the file and not allowed.
Messages: the summary lists the file, its grid, ω_L, polarization and the placement; the run reports whether the
field was interpolated or propagated, the band and grid of the propagator, the share of the energy in azimuthal modes
beyond `m_max` + 1 (an envelope mode k contributes to the code's modes |k ± 1|), and a warning if the envelope is not
small at the edge of the file. The time span of the file must fit into the window at `zeta_center` (warning).

Continuing a GEM-PIC-RZ run from its own `output.lasy` file: the exported times are t = s − ζ + T, with T = t_center +
zeta_center − file_z of the first file laser of that run (else 0). With the same `zeta_center` and `file_z = 0`, set
`t_center` = (T + s − zeta_center)/c to put every slice back where it was (for a run without a file laser:
t_center = (s − zeta_center)/c).

Checks (`tests/lasy.sh`, ctest `lasy`): a Gaussian at focus from a thetaMode file equals the built-in Gaussian to
1.2·10⁻⁴ of the peak (cartesian file 1·10⁻³, circular polarization 1.2·10⁻⁴); propagated from its focus at 0.3 mm
back to z = 0, its peak agrees with the paraxial Gaussian started before focus to 0.1 %; the exported level read back
gives the same field to 7·10⁻⁴; the luminal Bessel pulse of the `flying_focus` block written as a file and fed with
`inject = boundary` reproduces the built-in flying focus after 0.3 mm (peak 0.25 %, rms of the on-axis E_x
difference 1.3·10⁻³, 0.5 mm: 0.3 % / 1.4·10⁻³). Serial, OpenMP, the emulated GPU and the MPI pipeline give the same fields.

## `flying_focus` (0 or 1 block)

A luminal Bessel pulse (focus-wave mode): a converging conical wave whose intensity peak moves with the window at c,
so that the driver neither dephases nor depletes. It is fed in from the side.
- **Field:** each frequency ω has the cone angle that gives k_z = ω − κ, with κ = ω0(1 − cos θ). Then
  k_⊥² = 2ωκ − κ² − ω_p², and the pattern moves at exactly c:
  E_x = Re[e^{−iκs} Σ_ω a(ω) J0(k_⊥(ω) r) e^{iω(ζ−ζc)+i·cep}].
- **Spectrum and polarization:** a(ω) is Gaussian, so on the axis E_x = a0 exp(−(ζ−ζc)²/L²) cos(…). The
  polarization is linear along x (mode m = 1) by default, or circular (`polarization = circular`, `helicity` as
  in the `laser` block; a0 is then the amplitude of each of E_x, E_y).
- **Several components (beats):** `beat_cone_angle`, `beat_order`, `beat_amplitude`, `beat_phase` add components
  U_c = A_c e^{iψ_c} e^{−iκ_c s} Σ_ω a(ω) J_|l_c|(k_⊥,c(ω) r) e^{iω(ζ−ζc)+i·cep}, κ_c = ω0(1 − cos θ_c), of Bessel order
  l_c (azimuthal phase e^{i l_c φ}), with the same spectrum, each moving at c. The physical field is
  E_x = Re Σ_c U_c e^{i l_c φ} (linear) or E_x + i h E_y = Σ_c U_c e^{i l_c φ} (circular, helicity h); the main
  component has `order` (default 0), `cone_angle` and amplitude 1. Two components beat with the period
  2π/|κ_1 − κ_0| ≈ 4π/(k(θ_1² − θ_0²)) in z. A J0 core with a weak J1 e^{iφ} component (amplitude A ≪ 1) has its
  intensity maximum displaced by x_b ≈ A θ_1/(k θ_0²) from the axis, turning around the axis with that period: the
  driver of a helical bubble (a bubble whose axis winds helically in the lab).
- **Azimuthal modes:** with circular polarization and helicity h a component of order l lives in the single mode
  m = |l − h|; linear polarization puts it into m = |l − 1| and |l + 1|. So a J0 core plus J1 e^{ihφ} needs only
  m_max = 1 (core in m = 1, J1 part in m = 0); the summary lists the modes of every component, and modes above
  `m_max` are left out (warning).
- **Boundary condition (total field / scattered field):** in the radial absorber layer
  r_absorber < r < r_max, R relaxes to the incident field, R ← R_inc + d(r)(R − R_inc). The layer therefore feeds
  the converging wave and absorbs only what the plasma adds.
- **Choice of absorber strength:** use `grid.absorber_strength = 1`, so that R = R_inc at r_max. With the default
  0.2 the field near the axis fell by 9 % in 2 mm.

| key | kind | default | meaning |
|-----|------|---------|---------|
| `a0` | number | — | peak normalized amplitude on the axis |
| `cone_angle` | angle | — | θ at ω0 (core radius 2.405/(k sin θ), transverse wavelength λ/sin θ) |
| `tau_fwhm` *or* `length_1e` | time / length | — | on-axis intensity FWHM, or 1/e half length of the field |
| `zeta_center` | length | — | centre of the on-axis envelope |
| `cep` | angle | 0 | carrier-envelope phase |
| `background_density` | density | 0 | uniform plasma in which the incident wave propagates (the plasma in the layer must be this one); needs ω_p < ω sin θ for every component |
| `order` | integer | 0 | Bessel order l of the main component (J_\|l\| e^{ilφ}) |
| `polarization` | `linear` \| `circular` | `linear` | linear along x, or circular |
| `helicity` | +1 \| −1 | +1 | circular polarization: sense of rotation (as in the `laser` block) |
| `beat_cone_angle` | angles | — | cone angles of up to 3 further components |
| `beat_order` | integers | 0 | their Bessel orders l_c |
| `beat_amplitude` | numbers | 1 | their amplitudes relative to a0 |
| `beat_phase` | angles | 0 | their phases ψ_c (at s = 0) |

Requirements and checks:
- **Radial grid:** dr ≤ λ/(10 sin θ) everywhere (θ the largest cone angle). The conical wave has the same transverse wavelength at every
  radius, so the grid cannot coarsen outward.
- **Window:** it must reach ζc + 3L + r_max θ/2. The rings at radius r run ahead of the axial focus by about rθ/2.
- **Layer:** it must not contain plasma that differs from `background_density`, such as a neutral gas that gets
  ionized there.
- **Tests:** `benchmarks/flying_focus/` (vacuum: on-axis amplitude within 2 % over 2 mm, constant phase offset
  0.08 rad, group velocity c to 3·10⁻⁵). Several components: ctest `flying_focus_beat`
  (`benchmarks/flying_focus/beat_vacuum.deck`, results there): E_x within 0.15 % of a0, the J1 part within 8 %, the
  turning intensity maximum within 0.04 µm after 120 µm.
- **Half-step lag:** the field fed at level s is, after the step, the incident field of s − Δz/2: each component
  carries the extra phase κ_c Δz/2, so a beat pattern is turned by the constant angle (κ_1 − κ_0)Δz/2.
- **Backends:** serial, OpenMP and MPI give identical results; the GPU kernel is the OpenMP one and compiles with nvcc, but has not been run on a GPU yet.

## `plasma` (0 or 1 block)

The density is `n(r, z) = density · f(z) · g(r)` for `r < r_max`, and 0 outside.

| key | kind | default | meaning |
|-----|------|---------|---------|
| `density` | density | — | reference density n0 (also the normalization of the `n_e/n0` output) |
| `z_points`, `z_values` | length list, number list | constant 1 | piecewise-linear f(z), constant beyond the end points |
| `density_function` | function of x, y, z, r | 1 | a further factor, see [Functions](#functions-and-expressions) |
| `radial_profile` | `uniform` \| `parabolic` \| `table` | `uniform` | form of g(r) |
| `channel_depth`, `r_channel` | number, length | — | parabolic: g = 1 + depth · r²/r_channel² |
| `r_points`, `r_values` | length list, number list | — | table: piecewise-linear g(r) |
| `r_max` | length | open | radius of the plasma column |
| `smooth_length_m0` | length | = `smooth_length` | smoothing length of the mode m = 0 (the wake) if different from the modes m ≥ 1 (laser-driven currents). In a dense plasma the wake needs about 0.3 c/ω_p, the laser coupling only about 0.1 c/ω_p (`benchmarks/laser_evolution`) |
| `smooth_length` | length | 0 (off) | a > 0: radial smoothing of all sources of the streaming plasma (electrons, ions, freed electrons), j ← (1 − a²∇²)⁻¹ j per mode, which caps the singular density and E_z spike where the electrons close the bubble (README, "Bubble closure"). The beam list and the fields acting on the particles are not smoothed. The filter commutes with the discrete divergence, so the smoothed currents conserve the smoothed charge. Choose a at least the axis cell and the slice spacing; 2 k0⁻¹ (0.25 µm, 0.09 c/ω_p) in `benchmarks/lwfa_fbpic` |
| `free_electrons` | true/false | true | `false`: no free electrons at the start (neutral gas, ionized by the fields) |
| `ppc_r`, `ppc_phi` | integer | 2, 8 | macroparticles per radial cell (quiet start in r²) and per ring |
| `ppc_r_axis`, `ppc_phi_axis` | integer | 4, 16 | the same, for cells inside `r_axis_refine` |
| `r_axis_refine` | length | 10 k0⁻¹ | radius of the refined loading near the axis. It avoids sampling noise in the tiny cells on the axis |

Plasma electrons are always mobile. Ions are an immobile background unless `ion` blocks make them mobile (next section). Without a `plasma` block the propagation is in vacuum.

## `ion` (0 … N blocks, requires `plasma`)

Each block represents one ion species of the plasma by macroparticles. It serves two purposes:
- **Ion motion:** it makes the species mobile. This matters for high-amplitude lasers and for dense or very low-emittance bunches (ion collapse).
- **Field ionization:** it lets the species be field-ionized, level by level, with the ionization energies of the actual element from NIST.

| key | kind | default | meaning |
|-----|------|---------|---------|
| `element` *or* `mass` | chemical symbol / mass | — | `H`, `D`, `He`, `Li`, `Be`, `B`, `C`, `N`, `O`, `F`, `Ne`, `Na`, `Mg`, `Al`, `Si`, `Ar`, `Cu`, `Kr`, `Ag`, `Xe`, `Au` (standard atomic mass and NIST ionization energies of all levels), or an explicit mass such as `1 m_p` or `4.0026 u` |
| `charge_state` | integer | see below | initial charge state Z0. Without ionization: 1 for `H`, `D` and `mass`, required for heavier elements. With ionization: 0 (neutral atoms) |
| `fraction` *or* `density` | number / density | `fraction = 1` | amount of the species. `fraction` is the share of the positive charge that neutralizes the free electrons, Z0·n_atom = fraction·n_e (needs Z0 ≥ 1). `density` is the atom density at the reference density of the plasma profile; it is required for neutral atoms |
| `mobile` | true/false | true | `false` keeps the ions fixed. This saves time for heavy dopants; they are still field-ionized |
| `ionization` | `none` \| `adk` | `none` | ADK tunnel ionization, see below |
| `max_charge_state` | integer | atomic number | highest charge state reachable by ionization |
| `ionization_potentials` | energy list | NIST | ionization energies of the charge states 0, 1, 2, … (e.g. `552.07, 667.05 eV`). For an `element` they replace the NIST values (all Znuc values); for an ion given by `mass` they define the levels |
| `profile` | `plasma` \| `own` | `plasma` | `own`: the species has its own density profile (see below) |
| `name` | string | e.g. `N5+`, `He` | label |

- **Own profile (`profile = own`):** a dopant that is only in part of the plasma (e.g. a doped section for
  localized ionization injection, a gas jet in front of a cell) gets its own profile. Then `density` is the atom
  density at the reference point, and the profile keys of the `plasma` block are available in the `ion` block:
  `z_points`/`z_values`, `radial_profile` (with its keys), `r_max` and `density_function(x, y, z, r)`. The ions are
  loaded wherever this density is positive, also where the plasma has none. With `charge_state` Z0 > 0 the
  code also loads Z0 free electrons per atom at the same place (the species is neutral); `fraction` is not
  used. Example:
  ```
  ion { element = N;  ionization = adk;  charge_state = 5;  profile = own;  density = 1e17 cm^-3
        z_points = 0.5, 0.55, 0.95, 1.0 mm;  z_values = 0, 1, 1, 0 }   # N only between 0.5 and 1 mm
  ```
  Checks (`tests/dopant.sh`): an own profile equal to the plasma's reproduces the run bit by bit, also for a
  pre-ionized species with its electrons in place of the plasma electrons; a dopant at 3 < z < 6 k0⁻¹ is
  ionized only there, as strongly as a uniform one.
- **Profile and loading:** the atom density follows the plasma profile, n_atom(r, z) = density·f(z)·g(r) (or fraction/Z0 · n_e). Ions are loaded at the positions of the plasma electrons, so the loaded plasma is exactly neutral.
- **Neutrality:** the charge of all ion blocks, Σ Z0·n_atom, must not exceed the free-electron density. The rest of the positive charge stays an immobile background. Example: in a He/N mixture pre-ionized to He²⁺ and N⁵⁺, only the nitrogen needs a block, and He²⁺ is the background.
- **Push:** mobile ions are pushed through the window like the plasma electrons, with the relativistic Boris pusher and q/m = Z/M. Each block doubles the number of streaming macroparticles.

### Field ionization (`ionization = adk`)

- **Rate:** for every charge state k the code uses the ADK tunnel-ionization rate (Ammosov, Delone, Krainov 1986; s-state form with l* = n* − 1, as in Chen et al., JCP 236, 220 (2013), and FBPIC):

  w = ω_a C²_{n*} (U/2U_H) [2(U/U_H)^{3/2} E_a/E]^{2n*−1} exp[−(2/3)(U/U_H)^{3/2} E_a/E],

  with n* = (k+1)·(U_H/U)^{1/2} and C²_{n*} = 2^{2n*}/(n* Γ(n*+1) Γ(n*)). Here U is the NIST ionization energy of the level, and E is the instantaneous |E| at the ion. The laser carrier is resolved, so the static-field rate applies; there is no cycle averaging. For hydrogen the formula reduces to the Landau rate. Beyond the maximum of the ADK expression (barrier suppression) the rate stays at its maximum.
- **Monte Carlo:** each ion macroparticle has its own charge state. In every step of the ζ-march the rate is averaged over the two slices, (w(E_{n+1}) + w(E_n))/2, the trapezoidal rule in time. Sequential ionization k → k+1 → … is sampled with exponential waiting times, so several levels can ionize in one step.
- **Freed electrons:** each event creates an electron macroparticle at the ion position with the ion velocity, carrying the charge of one electron per represented atom. Freed electrons are tracked separately: trapped ones appear as `inj_<name>` in `bunches.txt`, and as species −30−s in the particle files.
- **Charge-state output:** `ionization.txt` gives, per step, the charge-state distribution of the ions leaving the window and the charge freed in that step.
- **Not included:** the energy that ionization takes from the laser (ionization current) is neglected, which is fine for dilute dopants. Collisional ionization is not modelled.

The plasma block key `free_electrons = false` starts from a neutral gas. All electrons then come from field ionization, and all ion blocks must have `charge_state = 0`.

## `bunch` (0 … N blocks)

| key | kind | default | meaning |
|-----|------|---------|---------|
| `name` | string | `bunchI` | label in the output |
| `species` | `electron` \| `positron` | `electron` | |
| `charge` | charge | — | total charge. Only the magnitude is used; the sign comes from `species` |
| `energy` *or* `gamma` | energy / number | — | kinetic energy, or the Lorentz factor |
| `energy_spread` | number | 0 | relative rms spread of p_z (e.g. `0.1 %`) |
| `sigma_x`, `sigma_y` | length | —, `sigma_x` | rms transverse size at the waist (Gaussian) |
| `distribution_r` | `gaussian` \| `function` | `gaussian` | transverse shape. `function`: `profile_r` and `r_extent` replace `sigma_x`, `sigma_y` |
| `profile_r` | function of x, y, r, zeta | — | transverse density (any normalization, ≥ 0; x, y relative to the bunch axis before `x0`, `y0`), sampled by rejection in the disk r < `r_extent` at the ζ of each particle. Rings, hollow or flat-top beams, profiles that change along the bunch. The summary prints the resulting rms sizes |
| `r_extent` | length | — | radius of the sampling disk of `profile_r` |
| `ux_function`, `uy_function` | function of x, y, r, zeta | 0 | mean transverse momenta [m c], added to the random part of the emittance; the energy of each particle is kept (u_z is reduced). E.g. a helical orbit u_φ = K: `ux_function = "-K*y/r"`, `uy_function = "K*x/r"` |
| `emittance_x`, `emittance_y` | length | 0, `emittance_x` | normalized rms emittance (`0.3 um` = 0.3 mm mrad) |
| `distribution_z` | `gaussian` \| `flattop` \| `function` | `gaussian` | longitudinal shape. `function`: `profile_z`, `zeta_min`, `zeta_max` (and optionally `gamma_function`), see [Functions](#functions-and-expressions); `zeta_center` is not used |
| `sigma_z` | length | — | gaussian: rms length |
| `length` | length | — | flat-top: full length |
| `zeta_center` | length | — | centre, measured from the rear end of the window |
| `x0`, `y0` | length | 0 | transverse offset. It excites modes m ≥ 1, which need m_max ≥ 1 |
| `z_waist` | length | 0 | position of the beta waist. The bunch is drifted ballistically back to z = 0 |
| `z_start` | length | 0 | the bunch is loaded only when the run reaches z_start (at the end of that step), at its deck position in the window and with its co-moving self-field, as at the start; then the beam current of the level is deposited anew. E.g. a witness loaded after the vacuum–plasma transition of the driver's wake has settled. A restart from a checkpoint that already holds the bunch does not load it again. Not yet available in the MPI pipeline |
| `macroparticles` | integer | 10000 | number of macroparticles |

**Helical ring** (`tests/bunch_functions.sh`): every electron on a circular orbit of radius R0 in a linear focusing
channel (ion column, plasma lens), as for an ion-channel laser with a single betatron amplitude:
```
constants { R0 = 1.06 um; K = 0.71 }
bunch {
  gamma = 50; charge = 3 pC; distribution_z = flattop; length = 1 um; zeta_center = 30 um
  distribution_r = function; profile_r = "gauss(r - R0, 20 nm)"; r_extent = 1.5 um
  ux_function = "-K * y / R0"; uy_function = "K * x / R0"     # rigid rotation: each radius on its own orbit
}
```
u_φ = K r/R0 puts every radius on its matched circle (k_β r = u_φ/γ); u_φ = K at all radii (`-K*y/r`) gives the
off-centre electrons a radial oscillation. In the test (R0 = 100 µm, γ = 50, matched plasma lens) the ring keeps its
radius to 0.03 % and its width to 0.1 % over 1/6 of a betatron period.

**Shot noise.** By default the macroparticles are placed at random, so the bunching noise at any wavelength is that
of N_macro independent particles, much larger than that of the real electrons. For FEL start-up from noise use

| key | kind | default | meaning |
|-----|------|---------|---------|
| `shot_noise` | `random` \| `fawley` \| `quiet` | `random` | `fawley`: beamlets of `beamlet` macroparticles spread evenly over one `noise_wavelength` (no bunching), displaced by random amounts such that the bunching at the wavelength and its harmonics h < beamlet/2 equals that of the real number of electrons, ⟨\|b_h\|²⟩ = 1/N_e (W. Fawley, PRST-AB 5, 070701 (2002)). `quiet`: the same beamlets without the displacements (quiet start: no bunching at h < beamlet; the amplifier starts from the seed only, as in codes without shot noise). Checked in ctest `beam_import` |
| `noise_wavelength` | length | — | the wavelength of the noise model, usually the resonant wavelength λ_r |
| `beamlet` | integer ≥ 3 | 4 | macroparticles per beamlet; `macroparticles` must be a multiple of it |
| `axis_refine` | integer 1 … 1000 | 1 | more macroparticles near the axis: the centres (particles, or beamlets with `fawley`) at r < `axis_radius` are sampled `axis_refine` times as often, each with 1/`axis_refine` of the charge. The charge density is unchanged, only its macroparticle noise there drops. Useful for the modes m ≥ 1, which vanish on the axis while the noise of the small cells there is largest (FEL harmonics). Bunches from files need resampling (`macroparticles` ≠ the number in the file) |
| `axis_radius` | length | 4 radial cells | radius of the refined region (with `axis_refine` > 1) |

With `axis_refine = 16` inside 4 µm (the 4 kA, σ_r = 10 µm beam of `benchmarks/harmonics/fel.deck`, 2.2 times the
macroparticles) the noise of E_z and E_φ in m = 2 in the first radial cell halves, as with 4 times the
macroparticles everywhere; the outer cells, the charge and the beam size are unchanged (`tests/axis_refine.sh`).

The beamlet members share x, y and the momentum. Checked in `tests/beam_import.sh`: N_e⟨\|b\|²⟩ = 1.08 (Fawley) and
N_macro⟨\|b\|²⟩ = 0.93 (random), against 1 ± 0.1.

### Bunches from files (`file`)

A bunch can be read from the particle output of an earlier stage, of this code or of any code that writes openPMD
(FBPIC, HiPACE++, WarpX, …). Then the shape keys (`charge`, `energy`, `sigma_*`, `emittance_*`, `distribution_z`)
are not used: everything comes from the file.

```
bunch {
  name        = witness
  file        = stage1/diag/openpmd/openpmd_2000.h5   # relative to the working directory or to the deck
  file_species = witness                  # needed if the file holds several species
  select      = "uz > 300 && abs(z) < 5 um"   # keep these (z: relative to the charge centroid)
  zeta_center = 25 um                      # centroid of the selected particles in the new window
  macroparticles = 200000                  # resample (default: the particles of the file as they are)
  shot_noise  = fawley;  noise_wavelength = 0.4 um
}
```

| key | kind | default | meaning |
|-----|------|---------|---------|
| `file` | path | — | openPMD file (one iteration per file or all in one) or a native `data_NNNNN.h5` |
| `file_species` | name | the only species | openPMD species; native files: a bunch name, `trapped` or `inj_<ion>` |
| `file_iteration` | integer | the last | openPMD iteration |
| `select` | function of x, y, z, ux, uy, uz, gamma | all | particles with a non-zero value are kept. z is the lab position relative to the charge centroid of the species, momenta in m_e c |
| `zeta_center` | length | — | the charge centroid of the selected particles is placed here |
| `x0`, `y0` | length | 0 | transverse shift |
| `charge` | charge | from the file | rescales the weights to this total charge |
| `macroparticles` | integer | as in the file | resample to this number: systematic resampling with probabilities ∝ charge, equal weights |
| `resample_smoothing` | number | 0 | Gaussian jitter of resampled particles, in units of the rms of each coordinate (avoids identical copies when the number grows; keep it small, it adds to the slice energy spread) |

- **Coordinates:** positions and momenta are converted to this run's units (the files carry SI factors). The
  longitudinal position is taken relative to the centroid, so the time of the file does not matter.
- **Species:** electrons and positrons only (charge and mass are checked).
- **Moments:** the summary prints the charge, energy, spread, rms sizes and emittances of the imported bunch, which
  also define its initial self-field (γ: the mean).
- **Restart:** a restarted run reads the beam from its checkpoint; the file is opened again only for the checks.
- **Check:** a bunch written by a run (openPMD or native) and read into a new run has the same moments to all
  printed digits (`tests/beam_import.sh`).

Each bunch starts with its relativistic vacuum self-field: E⊥ = −∇⊥ψ with −Δ⊥ψ = j_z, and B = ẑ × E. This assumes γ ≫ 1; a warning is printed for γ < 10. The field is initialized even if the plasma begins at z = 0: the bunch then enters the plasma with its vacuum field.

**Self-field of a bunch.** Every bunch starts with the field of a charge distribution moving with its own β:
−(Lap⊥ + γ⁻²∂ζ²)ψ = j_z/β, E⊥ = −∇⊥ψ, E_z = −γ⁻²∂ζψ, B = β ẑ × E. (For γ > 10⁵ this is the v = c field,
R = 2E⊥.) For γ ≲ 10⁴ in vacuum the box must hold this field, but not much more. The summary prints the smallest box for
every such bunch ("box for bunch …") and warns if the box is smaller:
- **Radial boundary** (a conducting wall): r_max ≥ γL, the reach of the field of a bunch of length L. A closer wall
  shields the longitudinal space charge like a beam pipe, g = 1 + 2 ln(min(r_max, γL)/a). For a γ = 100 bunch the
  energy spread changes by 1.4 % at r_max = γL and by 0.3 % at 2γL. With undulators also r_max ≥ 5w, where w is the
  size of the radiation at the end of the undulator. Switch the absorber off (`r_absorber > r_max`): it would damp
  the self-field.
- **Gaps at both ends of the window:** the scheme sets E_z = 0 at the front, so a field cut off there acts as a
  uniform E_z on the whole bunch. The gap needed for E_z(d)/E_z(edge) = 1 % is
  d = min( 10 √(a·max(a, γL))/γ , ln(100) r_max/(2.405 γ) ). The first term is for free space, the second for the
  field inside the wall. The gap behind the bunch is checked at the end of the run, after the slippage
  ∫(1 + K_eff²)/(2γ²) dz.
- **No check in a plasma:** there the field is screened over c/ω_p.

## `undulator` (0 … N blocks)

A static undulator field that acts on all particles (bunches, trapped and plasma particles). The field solver does not see it. Several blocks add up, e.g. two sections with a gap.

| key | kind | default | meaning |
|-----|------|---------|---------|
| `type` | `planar` \| `helical` | `planar` | planar: B_y = B0 cosh(k_u y) cos(ψ), B_z = −B0 sinh(k_u y) sin(ψ). Helical: B_x = B0 cos ψ, B_y = B0 sin ψ (on-axis field). ψ = k_u (z − z_start), B0 = K k_u m_e c/e |
| `period` | length | — | undulator period λ_u |
| `K` *or* `K_function` | number / function of z | — | peak deflection parameter, K = 0.934 B0[T] λ_u[cm]; `K_function`: a taper K(z) |
| `periods` *or* `length` | integer / length | — | number of periods, or magnetic length (including the ramps) |
| `z_start` | length | 0 | entrance |
| `ramp_periods` | integer | 1 | the field amplitude rises and falls linearly over this many periods at both ends. With whole periods the orbit has no net angle; a planar orbit is offset by K/(γ k_u). Helical: ≥ 1 |
| `natural_focusing` | true/false | true | planar only: the off-axis field (focusing in y with k_β = K k_u/(√2 γ)). `false`: B_y = B0 cos ψ everywhere |

**Resolution.** The summary prints the resonant wavelength λ_r = λ_u (1 + K²/2)/(2γ²) for every bunch (helical: 1 + K²). Warnings are issued when the resolution is too coarse:
- `dz` and `numerics.beam_dz_max` should be at most λ_u/16; the transverse current oscillates with the period.
- With undulators the beam is pushed in the mean of the fields of s⁰ and s¹ (`numerics.beam_field = mean`, the default). The self-field of an intense beam turns with the wiggle; pushed in the fields of s⁰ only, a 4 kA beam in a planar undulator (γ = 1000, K = 1.5, σ_r = 10 µm) blew up at λ_u/16 (ε_n 0 → 2.8 µm in 25 periods) and needed λ_u/32. With `mean` it stays at 0.04 µm at λ_u/16 (λ_u/64 converged value: 0.02 µm), and the push stays stable at λ_u/8.
- The particle push uses the undulator field averaged over its step along z (the field times sinc(k_u Δz/2), exact for a sinusoid). The orbit amplitude K/γ is then exact for any step; without this factor it is too large by (k_u Δz)²/24, 0.6 % at λ_u/16, which shifts the resonance by 0.6 %·K²/(1 + K²). The beam current is deposited with the momentum at the time of the position (the leapfrog push has it half a sub-step earlier, a phase lag k_u Δz/2 of the transverse current). With both, the FEL gain changes by less than 1 % in gain length when Δz or Δζ is halved from λ_u/16 and λ_r/16 (`benchmarks/genesis_seeded`).
- `dzeta` should be at most λ_r/12, better λ_r/20 or finer.

The (1,2,1) beam filter (`numerics.beam_filter`) damps the current at λ_r by cos²(k_r h/2), i.e. the radiated energy by cos⁴(k_r h/2) (−2 % at λ_r/32, −13 % at λ_r/12). In vacuum it can be switched off (`beam_filter = 0`).

**Geometry.** The radiation stays at fixed ζ, and the bunch slips back by one λ_r per period. Start the bunch near the front of the window; the window must hold the slippage N_u λ_r. Linear polarization lives in mode 1 (`m_max ≥ 1`), as for a laser.

## Transport elements: `quadrupole`, `plasma_lens`, `solenoid`, `chicane`

Between two stages (plasma → undulator) the beam is matched by magnets. These blocks add static fields that act on
all particles, like the undulators; the field solver does not see them. Several blocks of each kind are allowed.

| block | keys | field |
|-------|------|-------|
| `quadrupole` | `z_start`, `length`, `gradient` (T/m) | B_x = −G y, B_y = −G x. G > 0 focuses electrons in x and defocuses them in y (positrons: the reverse); use G < 0 for the other plane |
| `plasma_lens` | `z_start`, `length`, `gradient` (T/m), `radius` (optional) | active plasma lens: azimuthal B_φ = −G r inside `radius`, −G R²/r outside. G > 0 focuses electrons in both planes. Only the magnetic field of the discharge is modelled, not its plasma |
| `solenoid` | `z_start`, `length`, `field` (T), `ramp` (default 0) | B_z = B f(z) and B_r = −(r/2) B f′(z), with linear fringes of length `ramp` inside the length. Use `ramp > 0`: the radial fringe field gives the focusing (strength (B/2Bρ)² in the Larmor frame) |
| `chicane` | `z`, `r56`, `delay` (default \|r56\|/2), `energy` or `gamma` (reference) | thin chicane at z: every beam particle moves in ζ by r56·δ − delay, δ = γ/γ_ref − 1. r56 > 0: particles with more energy move forward (towards the head), as in a four-dipole chicane; a negative chirp (more energy at the tail) is then compressed. `delay` is the path difference of the reference particle (moves the beam back with respect to the radiation; ≈ R56/2 for a symmetric chicane) |

- **Thin chicane:** the shift is applied when the beam passes `z` (at the end of that step), without transverse
  dispersion or synchrotron radiation. The co-moving self-field of the beam moves with it (subtracted at the old,
  added at the new positions; vacuum is assumed at the chicane); the radiation in the window stays. About 0.6 % of the
  field energy of a γ = 1000 bunch remains at its old position. Not available in MPI runs with several ranks yet.
- **Resolution:** the beam push must resolve the elements: `numerics.beam_dz_max` ≤ 0.1/k_β (k_β = √(G/Bρ) or
  B/(2Bρ)), and an element should span ≥ 4 sub-steps. Warnings are printed otherwise. The step `grid.dz` may be larger.
- **Check (`tests/transport.sh`):** a cold bunch (γ = 1000) through a 100 T/m quadrupole, a plasma lens and a 5 T
  solenoid and a drift: the rms sizes at the end agree with the linear optics to 0.1 % (quadrupole x 0.38409 vs
  0.38374, y 1.64493 vs 1.64531; solenoid r 0.96466 vs 0.96466). A chicane with R56 = 100 µm halves a chirped 4 µm
  bunch (rms 0.5761 vs 0.5766 µm).

## `output` (0 or 1 block)

| key | kind | default | meaning |
|-----|------|---------|---------|
| `dir` | string | `diag` | output directory |
| `every` | integer | 25 | field/particle files every N steps (0: only at the start and end) |
| `map_stride` | integer | 4 | every n-th slice in the 2D maps |
| `axis`, `maps`, `particles` | true/false | true | which files to write (`particles` also controls the beam list in the HDF5 files) |
| `format` | one or more of `text`, `hdf5`, `openpmd` (`both` = `text, hdf5`) | `text` | `text`: `axis_`, `map_`, `beam_NNNNN.txt`; `hdf5`: one file `data_NNNNN.h5` per output step with all modes on the native grid (see below); `openpmd`: openPMD 1.1.0 files `openpmd/openpmd_<step>.h5` (see below). E.g. `format = text, openpmd` |
| `openpmd_dr`, `openpmd_dzeta` | length | smallest dr / dzeta of the grid | spacing of the uniform r–ζ grid of the openPMD fields |
| `openpmd_grid` | `uniform` \| `native` \| `both` | `both` | openPMD fields interpolated to a uniform grid (for every openPMD reader), on the code's own nodes (records `native_*`, exact), or both |
| `hdf5_precision` | `double` \| `single` | `double` | floating-point precision of the HDF5 datasets (single: half the size) |
| `hdf5_compression` | integer 0 … 9 | 0 | gzip level of the HDF5 datasets (0: none; the files of the examples shrink by about 3.5× with `single` + 6) |
| `fel` | `auto` \| `true` \| `false` | `auto` | FEL table `fel.txt` every step (auto: on when there are undulator and bunch blocks) |
| `fel_wavelength` | length | resonance of `undulator[0]` for `bunch[0]` | wavelength of the bunching factor in `fel.txt` |
| `fel_slice` | length | 20 `fel_wavelength` | slice length for the largest bunching `b_max` in `fel.txt` and for the power profile `fel_power.txt` (best a whole number of `fel_wavelength`) |
| `fel_harmonics` | integer 1…9 | 1 | > 1: also write `fel_harmonics.txt`, the power at the harmonics h = 1 … `fel_harmonics` of `fel_wavelength` per azimuthal mode m = 0 … `m_max` and slice, and the bunching at every harmonic per slice |
| `checkpoint_every` | integer | 0 | write a checkpoint `checkpoint_NNNNN.h5` every N steps and at the end of the run (0: none; needs HDF5) |
| `checkpoint_keep` | integer | 2 | keep only the newest N checkpoints (0: all) |
| `walltime_hours` | number | 0 | stop with a checkpoint (exit code 3) before this wall-clock time is exceeded (0: no limit) |
| `lasy` | true/false | false | write LASY envelope files `lasy/laser_NNNNN.h5` of the laser (see [LASY files](#lasy-files)) |
| `lasy_every` | integer | 0 | every N steps (0: only at the end of the run) |
| `lasy_dr`, `lasy_dt` | length, time | smallest dr, dzeta | uniform r and t grid of the file |
| `lasy_r_max` | length | `r_absorber` (r_max without absorber) | outer radius of the file |
| `lasy_wavelength` | length | wavelength of `laser[0]` (else `lambda0`) | carrier of the envelope |
| `lasy_bandwidth` | number 0 … 1 | 0.5 | spectrum kept: \|ω − ω_L\| < `lasy_bandwidth` ω_L (removes the wake, harmonics) |

Every run writes:

| file | content |
|------|---------|
| `deck_used.txt` | the deck after command-line overrides |
| `summary.txt` | the parameters in physical units |
| `bunches.txt` | one line per step and bunch: charge, energy, energy spread, sizes, emittances, bunch length; also the self-trapped plasma electrons (`trapped`) and the trapped electrons freed by ionization (`inj_<ion name>`). The energy spread (rms of γ over its mean, charge-weighted) is accumulated centred, so it is exact also for dE/E ≪ 1 at large γ (`tests/energy_spread.sh`: the printed value is the rounded value from the particles down to dE/E = 2e-6 at γ = 1000) and does not depend on the number of MPI ranks |
| `fel.txt` | (with `output.fel`) per step: field energy of the modes m ≥ 1 (the radiation of a centred beam) [J], the largest power through a slice, ∫S_z dA over the modes m ≥ 1 [W], and its ζ, the energy of mode 0 [J] (self-fields, wake), the resonant power `P_res_max` [W] and its ζ (see below), the bunching factor \|⟨e^{ikζ}⟩\| of every bunch at `fel_wavelength`, and its largest value over slices of `fel_slice` (slices with at least half the charge of the fullest one). Checked against the data files (energy) and the form factor of a short bunch (bunching) in `tests/undulator.sh`; MPI runs give the same table |
| `fel_power.txt` | (with `output.fel`) per step, the power of the radiation at `fel_wavelength` in every slice of length `fel_slice` [W]: the modes m ≥ 1 of R and L are demodulated at ±k over the slice, so the fields that move with the bunch (its Coulomb field, displaced by the helical or planar orbit, carries a Poynting flux of MW at kA currents) drop out. A seed in vacuum gives its power to 1e-4. This is the gain curve P(z) at fixed ζ, e.g. in the steady-state part of a long bunch |
| `fel_harmonics.txt` | (with `output.fel_harmonics` > 1) per step, harmonic h and mode m, the power [W] at wavelength `fel_wavelength`/h in every slice of `fel_slice`, demodulated at ±h k as in `fel_power.txt` (the total is the sum over m; h = 1 summed over m ≥ 1 is `fel_power.txt`); lines with m = −1 hold the bunching \|⟨e^{ihkζ}⟩\| of all bunches in the slice. Odd harmonics of a planar undulator radiate on the axis in m = 1; even harmonics are zero on the axis and radiate in m = 0, 2, … (an E_x ∝ x field has equal power in m = 0 and m = 2); m = 0 also holds the field of the bunching that moves with the beam, which is not radiation. Use a `fel_slice` that is a whole number of `fel_wavelength`: otherwise the slice sum of a uniform beam gives a spurious bunching of ~1/(number of wavelengths) |
| `ionization.txt` | per step and ionizing species: mean charge state and charge-state fractions of the ions leaving the window, freed charge |
| `axis_NNNNN.txt` | on-axis E_x, E_z, n_e/n0, Z n_i/n0 (mobile ions) versus ζ [µm] |
| `map_NNNNN.txt` | E_x, E_z, n_e/n0, Z n_i/n0 in the x–ζ plane (gnuplot `pm3d` format) |
| `beam_NNNNN.txt` | x, y, ζ [µm], u_x, u_y, u_z, q [pC], species (bunch index; −1 for self-trapped plasma electrons) |

### openPMD files

`format = openpmd` writes files that follow the openPMD standard 1.1.0 (file-based, `<dir>/openpmd/openpmd_<step>.h5`),
for openPMD-viewer, openPMD-api, VisualPIC, ParaView and other tools:

```python
from openpmd_viewer import OpenPMDTimeSeries
ts = OpenPMDTimeSeries('diag/openpmd')
Ex, info = ts.get_field('E', 'x', iteration=100, theta=0)        # plane theta = 0, all modes
z, uz, w = ts.get_particle(['z', 'uz', 'w'], species='witness', iteration=100)
```

- **Meshes:** `E`, `B` (components `r`, `t`, `z`) and `n_e` (density of all electrons; `Zn_i`: charge density of the mobile ions / e, with ion blocks), in SI via `unitSI`. Geometry `thetaMode` with dataset index 0 for m = 0 and 2m−1, 2m for the cos(mθ), sin(mθ) parts (`imag=+`).
- **Interpolation:** the standard needs a uniform grid. The fields are interpolated from the staggered (and possibly non-uniform) r–ζ grid to the nodes r_j = j·`openpmd_dr`, ζ_k = k·`openpmd_dzeta`, linearly as in the particle gather, with the regular behaviour on the axis. The defaults (the smallest steps) keep the resolution but make large files on stretched grids; choose coarser values for overviews. The native records (next item) and the `hdf5` format store the exact data.
- **Native records** (`openpmd_grid = native` or `both`): every mode on the code's own staggered, possibly non-uniform nodes, without interpolation, as scalar meshes `native_E_r`, `native_E_t`, `native_E_z`, `native_B_r`, `native_B_t`, `native_B_z` (faces or centres in r, slices in ζ) and `native_n_e`, `native_Zn_i` (densities / n0, centres × half-cells). They follow the proposal for rectilinear meshes in openPMD 2.0 ([`openpmd_rectilinear.md`](openpmd_rectilinear.md)): the node positions are data sets `<basePath>gridCoordinates/r_faces`, `r_centers`, `zeta_slices`, `zeta_halfcells`, referenced by the mesh attribute `gridCoordinates` (one path per axis). While the files are openPMD 1.1 the records carry `geometry = "other"` with `geometryParameters = "thetaMode;m=…;imag=+"`, so readers that know only regular grids skip them. E_r and B_φ on the axis hold the regular limit. Readers: `h5_extract` (every command except `beam` works on these files exactly as on `data_NNNNN.h5`) and `tools/openpmd_rectilinear.py` (Python; also reads QUARZ's native records, which keep the radial nodes in the attribute `gridCoordinates_r`).
- **Coordinates:** the lab position is z = ζ + s, and the time t = s/c (`gridGlobalOffset` and particle `positionOffset/z` hold s). Each step is an iteration.
- **Particles:** one species per bunch (its `name`), `trapped` (self-injected plasma electrons) and `inj_<ion>` (electrons from field ionization), if present: `position`, `positionOffset`, `momentum` (per real particle), `weighting` (real particles per macroparticle), constant `charge` and `mass`, and a `particlePatches` bounding box.
- **Precision and compression:** `hdf5_precision` and `hdf5_compression` apply to these files too.
- **Check:** the files pass the openPMD validator (`openPMD_check_h5`) without errors or warnings; the `openpmd` test also compares the fields rebuilt from the thetaMode data with the code's own reconstruction.

### LASY files

`output { lasy = true }` writes the laser at the plane z = s of a level as a LASY envelope file
`<dir>/lasy/laser_NNNNN.h5` (iteration NNNNN = step), at the end of the run (z_end, so the plane where the next stage
starts) or every `lasy_every` steps. The file can be read by LASY (`lasy.profiles.FromOpenPMDProfile(path,
"laserEnvelope")`), propagated, analysed or handed on to FBPIC, WarpX, HiPACE++, Wake-T, QUARZ, or back to
GEM-PIC-RZ (`laser.file`):

```python
from lasy.laser import Laser
from lasy.profiles.from_openpmd_profile import FromOpenPMDProfile
p = FromOpenPMDProfile("diag/lasy/laser_00200.h5", "laserEnvelope")
```

- **Field:** the forward field R/2 (the laser; the backward field L of the plasma is left out), interpolated as in
  the particle gather to r_j = j `lasy_dr` and to the times t_l = s − ζ + T (T: see the file lasers above), at 2 m_max + 4
  angles, demodulated at ω_L (analytic signal, spectrum |ω − ω_L| < `lasy_bandwidth` ω_L: the wake and the low
  frequencies drop out), and projected on the dominant polarization (eigenvector of the coherence matrix of E_x, E_y;
  the console line gives the share of the energy it captures, also stored in the root attribute `comment`).
- **Layout:** as LASY writes it: mesh `laserEnvelope`, complex, thetaMode components [2(m_max + 1) + 1][t][r] (the
  envelope modes |k| ≤ m_max + 1 that the code's modes 0 … m_max make), axes (t, r), SI units, attributes
  `angularFrequency`, `polarization` (real, or complex for elliptic polarization), `envelopeField = electric_field`.
- **MPI:** the fields are collected on rank 0 at these steps as for the other output files.

### HDF5 files

`data_NNNNN.h5` (NNNNN = step) contains the complete state needed for analysis:

| path | shape | content |
|------|-------|---------|
| `/` (attributes) | | `step`, `z`, `z_m`, `lambda0_m`, `k0inv_m`, `omega0_per_s`, `n0_m3`, `E0_V_per_m`, `B0_T`, `charge_unit_C`, `modes`, `Nr`, `Nzeta`, `dzeta`, `dz`, `backend` |
| `/grid/r_faces`, `/grid/r_centres` | [Nr+1], [Nr] | radii of the faces and centres of the (non-uniform) radial grid |
| `/grid/zeta_slices`, `/grid/zeta_halfcells` | [Nz+1], [Nz] | ζ of the field slices and of the half-cells |
| `/fields/Er`, `Bphi`, `Bz` | [M+1][Nz+1][Nr+1] | complex modes on the faces |
| `/fields/Ephi`, `Br`, `Ez` | [M+1][Nz+1][Nr] | complex modes on the centres |
| `/densities/ne_over_n0`, `Zni_over_n0` | [Nz][Nr] | electron density and ion charge density / n0 (m = 0), half-cells × centres |
| `/particles/beam/x`, `y`, `z`, `zeta`, `ux`, `uy`, `uz`, `charge`, `species` | [Np] | beam list: bunches and trapped particles (species as in `beam_NNNNN.txt`; attribute `bunch_names`) |
| `/input_deck` | string | the deck after command-line overrides |

- **Mode convention:** F(r, φ) = F₀(r) + 2 Re Σₘ Fₘ(r) e^{imφ}.
- **Complex numbers:** stored as the compound {r, i}. h5py reads them as `complex128`; C++ reads them into `std::complex<double>` (see `tools/h5_extract.cpp`).
- **Units:** all data are in code units. Every dataset has the attributes `unitSI` (multiply to get SI), `unit` and `description`.
- **Ghost slice:** slice Nz of the fields is the empty ghost ahead of the window.
- **Axis:** face 0 of E_r and B_φ holds 0. On the axis the fields follow from regularity, as in the particle gather.
- **Reconstruction:** `h5_extract <file> lineout <r_um> <phi_deg>` reconstructs E and B at any (r, φ).

`bunches.txt` and `ionization.txt` (one line per step) are always written as text.

### Checkpoints and restart

A checkpoint `checkpoint_NNNNN.h5` holds the complete state after step NNNNN,
in double precision:

- the raw field characteristics of level s⁰;
- all currents;
- the beam list, and the particles that continue streaming in the next step;
- the state of the random generator of the plasma loading;
- the step counters;
- the deck.

A restarted run continues exactly as the uninterrupted one. The tests check this value by value, including re-streaming, trapping and field ionization.

```
gempic_rz run.deck --restart diag/checkpoint_00100.h5    # from this file
gempic_rz run.deck --restart latest                      # newest checkpoint in output.dir
gempic_rz run.deck grid.z_end=2mm --restart latest       # continue further than planned
```

- **What must match:** the grid of the deck must be the one of the checkpoint (m_max, radial grid, ζ grid, lambda0), and so must the ion blocks. The restart stops with a message otherwise.
- **What may change:** `z_end`, `dz` or the z table, the output settings, the plasma profile ahead of the window, `iterations`, and the backend. A GPU checkpoint can be continued on the CPU and vice versa.
- **Ignored:** laser and bunch blocks; the state comes from the file.
- **Seed:** a different `seed` changes the ionization random numbers (warning).
- **Tables:** `bunches.txt` and `ionization.txt` are cut back to the restart step, so steps are never duplicated.
- **Files written:** the restart writes `deck_used_restart.txt` and `summary_restart.txt`.
- **Wall-time limit:** `walltime_hours` stops the run 1.5 mean step times before the limit, writes a checkpoint and exits with code 3. A job script can then resubmit:

```
# job.sh (e.g. SLURM), deck with  output { checkpoint_every = 100; walltime_hours = 23.5 }
R=""
if ls diag/checkpoint_*.h5 > /dev/null 2>&1; then R="--restart latest"; fi
gempic_rz run.deck $R
if [ $? -eq 3 ]; then sbatch job.sh; fi     # time limit reached: continue in a new job
```

- **Safe writing:** a checkpoint is written to a temporary file and then renamed, so a job killed while writing leaves the previous checkpoints intact.

## `numerics` (0 or 1 block)

| key | kind | default | meaning |
|-----|------|---------|---------|
| `backend` | `serial` \| `openmp` \| `cuda` | `cuda` in a GPU build, else `serial` | `serial`: reference implementation, one thread. `openmp`: parallel engine on CPU threads (CPU build; threads from `OMP_NUM_THREADS`). `cuda`: parallel engine on an NVIDIA GPU (build with `-DGEMPIC_CUDA=ON`). With MPI (several ranks): `serial` or `openmp` (threads on every rank); `cuda` is not yet available with MPI. See "Backends" below |
| `iterations` | integer | 2 | passes per slice (1 = predictor only) |
| `seed` | integer | 12345 | random seed (plasma loading, bunches, ionization) |
| `tridiag` | `auto` \| `thomas` \| `partition` | `auto` | tridiagonal solves of the field update in the parallel engine: `thomas`, one thread per azimuthal mode; `partition`, the Thomas/PCR partition method (Laszlo, Giles & Appleyard 2016) on all threads of the block. `auto` uses `partition` on the GPU and `thomas` on CPU threads. Both agree to rounding (1e-12); the serial backend always uses Thomas |
| `beam_dz_max` | length | 10 k0⁻¹ | largest sub-step of the beam push along s |
| `beam_pusher` | `vay` \| `boris` | `vay` | pusher of the beam list. Vay keeps the cancellation E + v × B = 0 of a relativistic beam in its own field exactly; Boris adds a spurious force ~ (qEΔt)² (`boris`: the behaviour of earlier versions) |
| `beam_field` | `auto` \| `level` \| `mean` | `auto` | fields of the beam push. `level`: the fields of s⁰ over the whole step (first order in the step). `mean`: after the march of s¹ the beam is pushed again from s⁰ in (F⁰ + F¹)/2 and its current of s¹ is deposited anew (the fields of s¹ are not solved again). `auto`: `mean` if the deck has an undulator, else `level`. Needed for planar undulators with intense beams (see the undulator section). In MPI runs the second push of a step is done at the start of the next step, when F¹ below the block of a rank is known; results as on one rank (tests/mpi.sh) |
| `vz_beam` | number | 0.9 | a streaming particle that reaches the next s-level joins the beam list (pushed along s) if v_z ≥ vz_beam; slower particles continue streaming in the next step, entering the ζ-march at their ζ |
| `trap_advance` | number ≥ 1 | 1 | a streaming particle leaves its s-level (and becomes trapped) when its lab position has advanced this many steps beyond the level. Values > 1 keep electrons that move forward by more than Δ without being trapped (laser, sheath) in their level and so avoid the switch of such particles between levels, a source of noise behind the bubble at small Δ (README, "Bubble closure") |
| `trap_gamma` | number ≥ 0 | 0 | a streaming particle that advances beyond its level becomes trapped only if its γ is at least this; slower ones stay in the streaming list and finish their path through the window in the fields of the level (as in a quasi-static code). Electrons quivering in a laser of amplitude a0 reach γ ≈ 1 + a0²/2 (v_z ≈ 0.9 c at a0 = 3) and advance far beyond Δ while crossing the pulse, but are not trapped; with 0 (default, as before) they are moved to the beam list and counted as trapped. A value somewhat above 1 + a0²/2 counts only electrons that are really accelerated by the wake |
| `theta` | 0.5 … 1 | 0.5 | implicitness of the R update in s, R¹ = R⁰ + Δ[(1 − θ)F⁰ + θF¹]: 0.5 is the trapezoid (second order); θ > 0.5 damps the spurious fast modes of the reduced equations by (1 − θ)/θ per step, at the cost of a damping (2θ − 1)(ωΔ)²/2 per step of the physical ones. Do not use with undulators (the phase of the FEL coupling becomes first order in Δz) |
| `theta_stiff` | number ≥ 0 | 0 | implicitness per radial eigenmode for the fields E_z, B_z in the R update: θ_k = (θ + μ_k)/(1 + μ_k), μ_k = c h Δ \|λ_k\| with c = `theta_stiff`, h the ζ step, λ_k the eigenvalues of the radial Laplacian of the mode. Modes the step resolves (c h Δ\|λ\| ≪ 1, all radiation) keep θ (no FEL damping); the stiff near-field modes of thin beam structures get θ_k → 1 (L-stable), which removes the undamped odd-even mode in s of the trapezoid there. 1 is a good value (results hardly change between 0.3 and 3). Cost: one pentadiagonal and four tridiagonal solves per slice and mode instead of two tridiagonal ones (about +50 % run time on CPU threads). 0 = off. See README, "Thin, intense beams" |
| `beam_consistent` | `true` \| `false` | `false` | beam sources consistent with the slice equations: the half-cell equations (E_z, L_r, L_φ) get the mean of the slice values at ζ_n, ζ_{n+1} of the beam current that the R update uses ((¼, ½, ¼) along ζ), and the R update gets the charge-conserving radial flux of the beam over the step (all modes, with the j_φ of the R_φ update) instead of (j_r⁰ + j_r¹)/2. Bunches also start with the co-moving field of their azimuthal current (B_z, E_φ of a rotating beam). Use together with `theta_stiff` for thin, intense, low-γ beams. Not yet available in the MPI pipeline |
| `beam_smooth_length` | length | 0 | radial smoothing of the beam current, as `plasma.smooth_length` (stopgap from before `beam_consistent`, kept for studies; one MPI rank only). 0 = off |
| `vz_min` | number | 0.8 | beam particles slower than this return to streaming in the same way. No particle is removed because it is slow |
| `beam_filter` | integer 0–8 | 1 | passes of a (1,2,1)/4 filter along ζ applied to the beam current (beam-list particles) before the field solve. It removes the source of numerical Cherenkov radiation at the ζ-Nyquist wavenumber (see README, "Numerical Cherenkov radiation"). Response cos^(2p)(kh/2): one pass attenuates kh = 0.2 by 1%. 0 = off |

### Backends

All backends run the same algorithm and read the same deck; `serial` is the
reference against which the others are tested.

- **Plasma loading and bunches:** these use the same random numbers in every backend.
- **`openmp`:** reproduces `serial` bit for bit, for any number of threads, when both are compiled without FMA contraction (`-DGEMPIC_NATIVE=OFF`). With `-march=native` there are last-bit differences in the complex arithmetic.
- **`cuda`:** the order of the atomic current sums varies from run to run, so results agree with `serial` to rounding (1e-9 relative on the regression deck) and are not bitwise reproducible.
- **Field ionization:** the random numbers come from a counter-based generator keyed by (seed, particle, step, slice, pass), in every backend. They do not depend on the thread schedule or the backend.
- **Timing:** the environment variable `GEMPIC_PROFILE=1` prints the time spent in each phase of the parallel engine at the end of the run.

## Command line

```
gempic_rz deck.txt                                  # run
gempic_rz deck.txt --check                          # read the deck, print the summary, exit
gempic_rz deck.txt numerics.backend=cuda           # run on the GPU (GPU build)
gempic_rz deck.txt --restart latest                # continue from the newest checkpoint
mpirun -np 8 gempic_rz deck.txt                    # MPI: 8 ranks along zeta (any number of ranks can restart a checkpoint)
gempic_rz deck.txt grid.z_end=2mm laser[1].a0=1.5   # override values (parameter scans)
gempic_rz deck.txt "bunch[0].name=witness 2" output.dir=scan_a
gempic_rz --template > my.deck                      # annotated example deck
```

An override replaces the value of an existing key, or adds the key to the block. Block indices start at 0.
