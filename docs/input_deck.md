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
| `name` | string | e.g. `N5+`, `He` | label |

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
| `sigma_x`, `sigma_y` | length | —, `sigma_x` | rms transverse size at the waist |
| `emittance_x`, `emittance_y` | length | 0, `emittance_x` | normalized rms emittance (`0.3 um` = 0.3 mm mrad) |
| `distribution_z` | `gaussian` \| `flattop` \| `function` | `gaussian` | longitudinal shape. `function`: `profile_z`, `zeta_min`, `zeta_max` (and optionally `gamma_function`), see [Functions](#functions-and-expressions); `zeta_center` is not used |
| `sigma_z` | length | — | gaussian: rms length |
| `length` | length | — | flat-top: full length |
| `zeta_center` | length | — | centre, measured from the rear end of the window |
| `x0`, `y0` | length | 0 | transverse offset. It excites modes m ≥ 1, which need m_max ≥ 1 |
| `z_waist` | length | 0 | position of the beta waist. The bunch is drifted ballistically back to z = 0 |
| `macroparticles` | integer | 10000 | number of macroparticles |

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
- `dzeta` should be at most λ_r/12, better λ_r/20 or finer.

The (1,2,1) beam filter (`numerics.beam_filter`) damps the current at λ_r by cos²(k_r h/2), i.e. the radiated energy by cos⁴(k_r h/2) (−2 % at λ_r/32, −13 % at λ_r/12). In vacuum it can be switched off (`beam_filter = 0`).

**Geometry.** The radiation stays at fixed ζ, and the bunch slips back by one λ_r per period. Start the bunch near the front of the window; the window must hold the slippage N_u λ_r. Linear polarization lives in mode 1 (`m_max ≥ 1`), as for a laser.

## `output` (0 or 1 block)

| key | kind | default | meaning |
|-----|------|---------|---------|
| `dir` | string | `diag` | output directory |
| `every` | integer | 25 | field/particle files every N steps (0: only at the start and end) |
| `map_stride` | integer | 4 | every n-th slice in the 2D maps |
| `axis`, `maps`, `particles` | true/false | true | which files to write (`particles` also controls the beam list in the HDF5 files) |
| `format` | one or more of `text`, `hdf5`, `openpmd` (`both` = `text, hdf5`) | `text` | `text`: `axis_`, `map_`, `beam_NNNNN.txt`; `hdf5`: one file `data_NNNNN.h5` per output step with all modes on the native grid (see below); `openpmd`: openPMD 1.1.0 files `openpmd/openpmd_<step>.h5` (see below). E.g. `format = text, openpmd` |
| `openpmd_dr`, `openpmd_dzeta` | length | smallest dr / dzeta of the grid | spacing of the uniform r–ζ grid of the openPMD fields |
| `hdf5_precision` | `double` \| `single` | `double` | floating-point precision of the HDF5 datasets (single: half the size) |
| `hdf5_compression` | integer 0 … 9 | 0 | gzip level of the HDF5 datasets (0: none; the files of the examples shrink by about 3.5× with `single` + 6) |
| `checkpoint_every` | integer | 0 | write a checkpoint `checkpoint_NNNNN.h5` every N steps and at the end of the run (0: none; needs HDF5) |
| `checkpoint_keep` | integer | 2 | keep only the newest N checkpoints (0: all) |
| `walltime_hours` | number | 0 | stop with a checkpoint (exit code 3) before this wall-clock time is exceeded (0: no limit) |

Every run writes:

| file | content |
|------|---------|
| `deck_used.txt` | the deck after command-line overrides |
| `summary.txt` | the parameters in physical units |
| `bunches.txt` | one line per step and bunch: charge, energy, energy spread, sizes, emittances, bunch length; also the self-trapped plasma electrons (`trapped`) and the trapped electrons freed by ionization (`inj_<ion name>`) |
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
- **Interpolation:** the standard needs a uniform grid. The fields are interpolated from the staggered (and possibly non-uniform) r–ζ grid to the nodes r_j = j·`openpmd_dr`, ζ_k = k·`openpmd_dzeta`, linearly as in the particle gather, with the regular behaviour on the axis. The defaults (the smallest steps) keep the resolution but make large files on stretched grids; choose coarser values for overviews. The native `hdf5` format stores the exact data.
- **Coordinates:** the lab position is z = ζ + s, and the time t = s/c (`gridGlobalOffset` and particle `positionOffset/z` hold s). Each step is an iteration.
- **Particles:** one species per bunch (its `name`), `trapped` (self-injected plasma electrons) and `inj_<ion>` (electrons from field ionization), if present: `position`, `positionOffset`, `momentum` (per real particle), `weighting` (real particles per macroparticle), constant `charge` and `mass`, and a `particlePatches` bounding box.
- **Precision and compression:** `hdf5_precision` and `hdf5_compression` apply to these files too.
- **Check:** the files pass the openPMD validator (`openPMD_check_h5`) without errors or warnings; the `openpmd` test also compares the fields rebuilt from the thetaMode data with the code's own reconstruction.

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
| `backend` | `serial` \| `openmp` \| `cuda` | `cuda` in a GPU build, else `serial` | `serial`: reference implementation, one thread. `openmp`: parallel engine on CPU threads (CPU build; threads from `OMP_NUM_THREADS`). `cuda`: parallel engine on an NVIDIA GPU (build with `-DGEMPIC_CUDA=ON`). See "Backends" below |
| `iterations` | integer | 2 | passes per slice (1 = predictor only) |
| `seed` | integer | 12345 | random seed (plasma loading, bunches, ionization) |
| `beam_dz_max` | length | 10 k0⁻¹ | largest sub-step of the beam push along s |
| `beam_pusher` | `vay` \| `boris` | `vay` | pusher of the beam list. Vay keeps the cancellation E + v × B = 0 of a relativistic beam in its own field exactly; Boris adds a spurious force ~ (qEΔt)² (`boris`: the behaviour of earlier versions) |
| `vz_beam` | number | 0.9 | a streaming particle that reaches the next s-level joins the beam list (pushed along s) if v_z ≥ vz_beam; slower particles continue streaming in the next step, entering the ζ-march at their ζ |
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
