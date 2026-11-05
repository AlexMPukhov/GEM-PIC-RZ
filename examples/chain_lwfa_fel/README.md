# Chained example: plasma injector → accelerator → chicane → FEL

A laser-wakefield injector, a long plasma accelerator, beam transport and an undulator, chained with openPMD
particle files:

| stage | code | file | what |
|---|---|---|---|
| 1 | FBPIC | `stage1_fbpic.py` | density down-ramp injection, 0.2 mm |
| hand-off | python | `handoff.py` | position of the injected beam behind the laser |
| 2 | GEM-PIC-RZ | `stage2_accel.deck` | acceleration in a guiding channel to 1.4 mm, exit ramp |
| 3 | GEM-PIC-RZ | `stage3_fel.deck` | energy cut, shot noise, plasma lens, chicane, helical undulator |

`run_chain.sh` runs all stages in a work directory and prints the beam after each import and the FEL power along the
undulator. The parameters are a generic demonstration sized for a workstation; they do not reproduce a published
case.

## Why two codes

GEM-PIC-RZ traps electrons in a density down-ramp, but the trapped charge is not reliable there
(`benchmarks/downramp_fbpic`). Injection in a down-ramp is therefore computed with a full PIC code, here FBPIC. (Field
ionization injection is computed by GEM-PIC-RZ itself, `examples/ionization_injection.deck`.) GEM-PIC-RZ takes over
for the long acceleration (large steps along the propagation), the transport elements and the FEL without period
averaging, with the full space charge. Any code that writes openPMD particles can be stage 1.

## Stage 1: down-ramp injection (FBPIC)

The case of `benchmarks/downramp_fbpic` with a guiding channel: laser λ0 = 0.8 µm, a0 = 1.8, w0 = 8.02 µm,
τ = 10.6 fs, focus at the entrance; density 0 → 2 n0 over 20 µm, plateau to 60 µm, down-ramp to
n0 = 3.48e18 cm⁻³ (0.002 n_c) at 70 µm; parabolic channel n = n0 (1 + 0.5 r²/w0²), matched to the waist
(Δn(w0) = 1/(π r_e w0²) = 0.5 n0). FBPIC 0.27, 2 modes, Δz = λ0/24, Δr = 0.2 µm, 6000 steps to z = 0.2 mm.
Output: the electrons with u_z > 10 (`stage1/hdf5/data00006000.h5`) and the lab position of the laser centre
(`stage1/laser_centre.txt`).

- Run time: 1 min 49 s on an RTX 5090 (`python3 stage1_fbpic.py 200`; `--cpu` without a GPU).
- Result (γ > 20): Q = 20.7 pC, ⟨E⟩ = 14.0 MeV, rms ΔE/E = 14 % (the head, injected first, has the most energy),
  σ_z = 0.71 µm, σ_x = 0.22 µm, ε_n,x = 0.084 µm, centroid 13.5 µm behind the laser centre.

## Hand-off

```
bunch {
  name        = witness
  file        = stage1/hdf5/data00006000.h5   # openPMD (species 'electrons'); positions in SI, any time
  select      = "gamma > 20"                  # function of x, y, z, ux, uy, uz, gamma (z: from the centroid)
  zeta_center = 17.02 um                      # where the charge centroid of the selected particles goes
}
```

- Positions are taken relative to the charge centroid, so the time of the file does not matter; `zeta_center`
  places the beam in the new window. `handoff.py` computes it from the distance between the laser centre and the
  beam centroid in the file: ζ_beam = ζ_laser − (z_laser − z_beam), with ζ_laser = `laser.zeta_center` of stage 2.
- The laser of stage 2 starts again from the analytic Gaussian with the same parameters, at its focus. (An evolved
  pulse can be handed over as a LASY file, `laser { file = … }`, see `docs/input_deck.md`.)
- The beam starts with its co-moving self-field (γ of the beam).
- `macroparticles` resamples the beam (equal weights), `shot_noise = fawley` loads the shot noise of the real number
  of electrons at `noise_wavelength` (stage 3). `charge` would rescale the weights.

## Stage 2: acceleration (GEM-PIC-RZ)

Plasma n0 with the channel to 1.2 mm, exit ramp to zero at 1.4 mm; the grid of `tests/regression.deck` (2 modes,
Δζ = λ0/25, 100 stretched radial cells), Δz = 2.5 µm, 560 steps. Output: openPMD particles of the beam at the exit
(`stage2/openpmd/openpmd_560.h5`, species `witness`).

- Run time: 245 s on 32 threads (AMD EPYC 7543, `numerics.backend = openmp`), 0.33 GB.
- Result: Q = 20.1 pC, ⟨E⟩ = 141 MeV, rms ΔE/E = 11 %, mostly a chirp (γ from 173 at the tail to 309 at the head),
  ε_n,x = 0.59 µm, σ_z = 0.83 µm. The head is nearly monoenergetic: γ > 290 carries 9.2 pC at γ = 302 with rms
  ΔE/E = 1.9 %, slice spread 0.5 % (50 nm slices), σ_z = 0.30 µm, peak current 3.8 kA.
- The projected emittance of the whole beam rises to 2 µm in the first 0.4 mm, while the beam adapts to the wake of
  the restarted laser, and settles at 0.6 µm.

## Stage 3: transport, chicane and FEL (GEM-PIC-RZ)

- **Import:** `select = "gamma > 290"`, resampled to 300 000 macroparticles with Fawley beamlets of 4 at
  λ_r = 109.6 nm (`shot_noise = fawley`).
- **Plasma lens** (2000 T/m, z = 2–3 cm) catches the diverging beam (3 mrad).
- **Chicane** at z = 4 cm, R56 = +100 µm, `delay = 0` (no radiation yet). R56 > 0 moves the higher energies forward;
  the chirp of this beam (more energy at the head) is therefore stretched: σ_z from 0.48 to 2.37 µm. The slice
  energy spread and the chirp per cooperation length drop by the stretch factor, the current by the same factor.
- **Undulator:** helical, period 1 cm, K = 1, 80 periods from z = 5 cm, λ_r = λ_u (1 + K²)/(2γ²) = 109.6 nm.
- **Grid:** 16 slices per λ_r (6569 slices, 45 µm window: gap, beam, slippage 80 λ_r, gap), 254 radial cells
  (0.5 µm at the axis for the 0.7 µm beam at the import, 3 µm to 0.3 mm for the collimated beam and the radiation,
  60 µm at the wall), r_max = 2.5 mm ≥ γL for the space-charge field, absorber off; 16 steps per period
  (1360 steps).
- Run time: 725 s on 32 threads, 3.1 GB.
- **Result** (`stage3/fel.txt`, P_res: power at λ_r): 0.9 kW at the undulator entrance, 0.43 MW after 80 periods,
  power gain length about 9 cm between z = 0.1 and 0.3 m. Without the chicane (`chicane.r56=0um`): 0.05 MW.

**Limitations.** The transport is crude on purpose (one lens): the chromatic focusing of the 2 % spread raises
ε_n,x from 0.6 to 1.55 µm and the beam keeps diverging in the undulator (σ_x from 84 to 163 µm, natural focusing of
the helical undulator β ≈ 0.7 m), so the gain slows down after 0.3 m. A matched lattice (quadrupole triplet,
focusing along the undulator) is the next step for a quantitative case. Trapped charge and beam quality of stage 1
depend on the injector's resolution.

## Running it

```
mkdir work && cd work
GEMPIC=/path/to/build/gempic_rz OMP_NUM_THREADS=32 sh /path/to/examples/chain_lwfa_fel/run_chain.sh
```

Stage 1 needs FBPIC (a GPU is recommended); if `stage1/hdf5/data00006000.h5` and `stage1/laser_centre.txt` exist
(e.g. computed on a GPU machine), the script starts with the hand-off. `handoff.py` needs numpy and h5py.

## Adapting it to your case

- **Another injector:** any openPMD particle file. Set `file_species` if the file holds several species,
  `file_iteration` if it holds several iterations; choose the beam with `select`.
- **Placement:** keep the beam's distance behind the laser (`handoff.py`) when a laser stage follows; for a vacuum
  stage put it near the front of the window, leaving room for the slippage N_u λ_r behind the radiation front.
- **Vacuum stages:** the summary line "box for bunch" prints the radial size (r_max ≥ γL) and the gaps at both ends
  that the self-field needs; switch the absorber off. Resolve the elements: `numerics.beam_dz_max` ≤ 0.1/k_β.
- **Radial grid in the undulator:** a stretched grid must not reflect the radiation (README "Free-hand
  grids"): keep k_⊥ dr ≲ 1 wherever the radiation must propagate, grow dr only outside.
- **Chicane sign:** r56 > 0 moves higher energies to the head. A beam with more energy at the tail is compressed,
  one with more energy at the head (as here) is stretched.

The reduced chain (2 steps of `tests/regression.deck` → import → chicane → 10 undulator periods) is the ctest
`chain` (`tests/chain.sh`).
