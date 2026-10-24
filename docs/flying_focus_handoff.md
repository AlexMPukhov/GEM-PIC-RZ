# Flying focus with GEM-PIC-RZ: handover (state of 2026-10-21)

For the student continuing this project, and for his or her Claude instance. Start a session with:
*"Read CLAUDE.md and docs/flying_focus_handoff.md, then continue with the open tasks."*
The project leader is Alexander Pukhov (HHU Duesseldorf, TP I); physics decisions and publication questions go to him.

## 1. Goal
Simulate laser-wakefield acceleration driven by a **flying focus**. In such a driver the intensity peak moves along
the axis at a velocity v_f that can be set by the optics, e.g. v_f = c, so that the wake does not dephase from the
electrons. The energy comes in from the side through a converging (axicon-like) wave, which continuously refills
the region near the axis, so the driver does not deplete. The questions are:
- Does the wake stay stationary over long distances?
- What energy gain is reachable?
- How do the off-axis rings, the ionization of the gas and the plasma response limit the scheme?

Related literature (check before citing):
- Palastro et al., PRL 124, 134802 (2020), "Dephasingless laser wakefield acceleration";
- Caizergues et al., Nat. Photon. 14, 475 (2020), "Phase-locked laser-wakefield electron acceleration"
  (axiparabola);
- Froula et al., Nat. Photon. 12, 262 (2018) (flying focus);
- Smartsev et al., Opt. Lett. 44, 3414 (2019) (axiparabola).

## 2. The code in two paragraphs
**GEM-PIC-RZ** (README.md; code paper in `paper/gempic_rz.tex`) is a Galilean electromagnetic PIC code in
quasi-cylindrical geometry.
- **Variables:** ζ = z − ct is the window coordinate and s = z is the slow coordinate (GT2): a level s is the plane
  z = s of the laboratory.
- **Fields:** they are expanded in azimuthal modes m. The forward field R = E_⊥ + B_⊥ (rotated) obeys ∂_s R = …,
  so the laser carrier is resolved in ζ. The backward and longitudinal fields follow from their ζ-equations.
- **What it is not:** the code is **not** quasi-static. It keeps the slow derivative of the forward field.
  Alexander insists on this wording: never call it "quasi-static" or a "quasi-static approximation".
- **Radial grid:** a non-uniform finite-element grid.
- **Particles:** plasma electrons stream through the window; beams and trapped particles move with it.
- **Ionization:** ADK ionization of neutral gases and ions is built in.

**The flying-focus extension** (added 2026-10-21, `flying_focus` block, docs/input_deck.md):
- **Incident field:** an ideal **luminal Bessel pulse** ("focus-wave mode"). Each frequency ω gets the transverse
  wavenumber that makes k_z = ω − κ, with κ = ω0 (1 − cos θ). Then k_⊥² = 2ωκ − κ² − ω_p², and in the window
  variables
  `R_r,1 = E_x = Re[ e^{-iκs} Σ_ω a(ω) J0(k_⊥(ω) r) e^{iω(ζ-ζc)+i·cep} ]`, `R_φ,1 = i E_x`
  (linear polarization along x, mode m = 1; a(ω) is Gaussian). The whole pattern moves with the window at c.
  |E| does not depend on s; only the carrier phase slips by κs (phase velocity c / cos θ).
- **Boundary condition, total field / scattered field:** in the radial absorber layer r_absorber < r < r_max,
  every step sets `R ← R_inc + d(r) (R − R_inc)`, with d = 1 − strength·x² and x = (r − r_a)/(r_max − r_a). The
  layer thus feeds the converging wave, and only the field that the plasma adds (scattered or radiated) is absorbed.
  With `absorber_strength = 1`, R = R_inc at r_max exactly. Use 1: with 0.2 the on-axis field fell by 9 % in 2 mm.
- **Background plasma:** `background_density` computes the incident wave in a uniform plasma; it must equal the
  plasma density in the layer.

## 3. Where it is implemented
| file | what |
|---|---|
| `src/config.hpp` | `FlyingFocusConfig` (a0, theta, length_1e, zeta_center, cep, background; `kappa()`) |
| `src/config.cpp` | `readFlyingFocus` (keys, checks: m_max ≥ 1, absorber layer present, dr ≤ λ/(10 sin θ) warning, window warning, ω_p < ω sin θ), summary line |
| `src/simulation.cpp` | `setupFlyingFocus()` builds G(r, ζ) on faces and centres (241 frequencies over ±6/L, Bessel J0 via `std::cyl_bessel_j`); initial field in the constructor; `absorber(double s)` with TF/SF (serial reference) |
| `src/simulation.hpp` | `ffGf`, `ffGc`, `ffInc(face, n, i, s)` |
| `src/pipeline.cpp` | `absorbSlice(F, n, s)`: the same for the MPI pipeline |
| `src/engine.cpp` | the same in the OpenMP/GPU kernel (end of the step, "absorber"); G uploaded as `E.ffGf/ffGc` |
| `tests/flying_focus.sh` | ctest `flying_focus`: vacuum, 0.5 mm, peak on the axis 2 ± 2 %, envelope stays at ζc ± 0.05 µm |

Serial, OpenMP and 2 MPI ranks give **identical** results. The GPU kernel compiles with nvcc but has **not yet run
on a GPU**: run `ctest` in a GPU build first (`benchmarks/gpu/`, `docs/gpu.md`).

## 4. Input decks
Tests from 2026-10-21 (`benchmarks/flying_focus/`, README there): λ = 0.8 µm, a0 = 2 on the axis, 50 fs,
θ = 50 mrad (core radius 6.1 µm, transverse wavelength 16 µm), dr = 0.5 µm, Δζ = λ/24, 2 CPU cores.
| deck | content | result |
|---|---|---|
| `vacuum.deck` | vacuum, 2 mm | amplitude within 2 %, constant carrier phase offset 0.08 rad (numerical dispersion), envelope moves at c to 0.07 µm / 2 mm |
| `gas.deck` | neutral H with an **artificial 200 eV** potential, gas r < 55 µm, plasma-free layer 60–80 µm, 1 mm | ionized channel r ≈ 10–20 µm; laser on the axis within 3 % (no depletion); first wake period stationary; behind it irregular (hot electrons from a 200 eV potential: artefact) |
| `plasma.deck` | preformed plasma 1e18 cm⁻³, incident wave in plasma (`background_density`), 1 mm | first wake period stationary from 0.4 mm on (zero crossing 60.3–61.4 µm); on-axis peak 2.0 → 1.8 in 0.5 mm, then constant; behind the first period the ring wakes cross |

Production decks (`benchmarks/flying_focus/production/`, every line commented):
- `ff_plasma.deck`: **start here**. Hydrogen as a preformed plasma, r_max 150 µm, Δζ = λ/32, 5 mm,
  ~3 GB, ~12 s per step on 2 cores (2000 steps): cluster (MPI) or one GPU.
- `ff_gas_highIp.deck`: neutral gas ionized near the axis only (200 eV test gas), layer in neutral gas.
- `convergence.sh`: one-parameter variations (dr, Δζ, Δz, r_max, particles per cell) over 2 mm.

Analysis (`benchmarks/flying_focus/`):
- `check_axis.py <dir> a0 theta tau zc`: laser on the axis against the ideal pattern;
- `check_wake.py <dir> [zc]`: laser amplitude and phase, peak E_z, position of the first zero crossing of E_z
  behind the pulse, density on the axis;
- `check_map.py`: field map against the incident field by radius band;
- `plot_gas.py <dir>`: 4-panel figure.

## 5. Lessons (read before changing parameters)
1. **Neutral boundary vs. real hydrogen.** The incident intensity falls off only as I(r)/I_axis ≈ 2/(π k_⊥ r).
   - Hydrogen (13.6 eV) is ionized by 50 fs at a ≈ 0.005, so at a0 = 2 the gas would be ionized out to ~0.25 m.
     A boundary in neutral gas is impossible.
   - Because that ionization is early and cold, real H = preformed plasma, with `background_density`.
   - A neutral-gas boundary works only for high ionization potentials. Estimate where the gas ionizes with ADK
     before choosing r_max; the ADK estimate script is reproduced in `benchmarks/flying_focus/README.md`.
2. **Radial resolution.** The conical wave has the same transverse wavelength λ/sin θ at every radius, so the
   non-uniform grid cannot coarsen outward. Use dr ≤ λ/(10 sin θ) everywhere. The number of radial cells grows
   with r_max.
3. **Window.** At radius r the converging wave runs ahead of the axial focus by about rθ/2. The window needs
   ζc + 3L + r_max θ/2 in front (the code warns).
4. **Layer.**
   - Use `absorber_strength = 1`.
   - The layer must not contain plasma that differs from `background_density`.
   - Its response should stay linear (a ≲ 0.3).
   - Make it several transverse wavelengths wide.
5. **Numerical dispersion.** In vacuum the carrier phase settles at an offset of ~0.08 rad (λ/24, dr = 0.5 µm),
   about 1 % in κ. It should fall with Δζ² and dr²: check with `convergence.sh`.
6. **Ring wakes.** Every Bessel ring drives its own plasma response; behind the first period these wakes cross and
   the wake on the axis becomes irregular. This is probably physics, not numerics; check it with
   `convergence.sh`, and compare with the literature.
7. **Shell.** Never run `pkill -f <pattern>` in the session shell: it can kill the shell itself. Use
   `ps aux | awk '/pattern/ && !/awk/ {print $2}' | xargs -r kill`.

## 6. Open tasks (suggested order)
1. **GPU / cluster check:** `ctest` in the GPU build (includes `flying_focus`); one short `ff_plasma.deck` run on
   GPU vs. CPU.
2. **Convergence** of `ff_plasma.deck` (`convergence.sh`): laser on the axis, wake amplitude and phase, ring wakes.
3. **Physics run:** 5 mm or more in `ff_plasma.deck`. Measure the stationarity of the wake (E_z(ζ) at all z) and
   the laser energy budget. Compare with a Gaussian driver of the same on-axis a0 and energy in the core, which
   dephases and depletes: that is the key figure.
4. **Witness / injection:** a test bunch (`bunch` block, `passive`?) or ionization injection (dopant with
   `profile = own`) to measure the energy gain without dephasing.
5. **Focus velocity v_f ≠ c:** generalize the incident field to k_z = ω/v_f − κ. Then
   k_⊥² = ω² − ω_p² − (ω/v_f − κ)², and the pattern moves at v_f relative to the window. Code: `setupFlyingFocus`
   and `ffInc`, where the s-phase becomes frequency dependent, so G(r, ζ, s) must be summed per step in the layer.
   Then study v_f slightly above c (to compensate for the wake phase velocity).
6. **Finite flying focus from a real optic** (axiparabola + echelon or chirp): now possible with a LASY file and
   `laser { file = ...; inject = boundary }` (branch `lasy`, docs/input_deck.md "Lasers from LASY files"): LASY models
   the optic, GEM-PIC-RZ propagates the file exactly (Hankel + FFT) into the box and the layer at every step. Checked
   against the built-in luminal Bessel pulse (ctest `lasy`).
7. **Radial polarization** (m = 0, strong on-axis E_z) as an option.
8. Write a validation section and a figure for the README, as for the other benchmarks.

## 7. Working rules (from CLAUDE.md)
- Build and test: `cmake -B build && cmake --build build -j && (cd build && ctest)` (27 tests, ~20 min on
  2 cores). In containers use `OMP_PROC_BIND=false` and `numerics.backend=openmp` for threads.
- **The repository is public** (github.com/AlexMPukhov/GEM-PIC-RZ): everything pushed is visible. Do not push
  unpublished results without Alexander's consent; work on a branch (e.g. `flying-focus`) and open pull requests
  only when Alexander agrees.
- Commit at milestones with descriptive messages. Every new feature gets a test (ctest) against an independent
  reference, and an entry in docs/.
- Big outputs (HDF5 dumps, `*_out/` directories) stay out of git.
