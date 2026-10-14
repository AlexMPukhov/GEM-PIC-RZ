# Benchmarks against other codes

| directory | reference code | case | result |
|---|---|---|---|
| `genesis_seeded/` | GENESIS 1.3 v4 (steady state) | seeded FEL amplifiers at 10 nm and 0.8 nm | power gain length within 4 % (10 nm) and 0.6 % (0.8 nm) of GENESIS run with the same physics |
| `lwfa_fbpic/` | FBPIC 0.27 (full PIC, spectral, quasi-cylindrical) | LWFA, a0 = 2, n0 = 3.5·10¹⁸ cm⁻³ | on-axis wake E_z of the first period to 1.6 % rms; laser field to 1 % |

The runs were made on 2 CPU cores. The decks, input files and scripts are here; the GEM-PIC-RZ output
tables (`fel_power.txt`, `fel.txt`, `bunches.txt`, axis files), the GENESIS outputs (`*.out.h5`) and the FBPIC lineouts are included.

## Seeded FEL against GENESIS 1.3 v4

Both codes amplify a 10 MW Gaussian seed (w0 = 30 µm, circular polarization) in a helical undulator, with a
round beam of 4 kA, σ_r = 10 µm, no emittance and no energy spread. Neither has natural focusing.

| case | γ | λ_u | K | periods | λ_r | ρ |
|---|---|---|---|---|---|---|
| `gem.deck` (10 nm) | 1000 | 1 cm | 1 | 60 | 10 nm | 9.1·10⁻³ |
| `gem_x.deck` (0.8 nm) | 5000 | 2 cm | 1 | 140 | 0.8 nm | 2.9·10⁻³ |

**Steady state in GEM-PIC-RZ.** GENESIS solves a single slice. GEM-PIC-RZ is time dependent, so the bunch is a
long flat top (2 µm = 200 λ_r at 10 nm, 160 nm = 200 λ_r at 0.8 nm) with smooth edges, and the seed is flat over
it. The radiation slips over the electrons by N_u λ_r. The radiation in a ζ range that met only the flat part of
the bunch behaves as in a steady-state run: ζ = 2.6 … 3.6 µm and 0.32 … 0.36 µm. Its power is read from
`fel_power.txt`, which demodulates the fields at λ_r, so the Coulomb field moving with the bunch is excluded.
Shot noise is the Fawley beamlet noise of the real electron number (`shot_noise = fawley`), negligible against
the seed. Resolution: Δζ = λ_r/16, Δz = λ_u/16, `m_max = 1`. The windows follow the box criteria of the code.

**GENESIS** (`genesis/fel.in`, `genesis_x/fel.in`): steady state, 8192 particles, field grid ±0.3 mm with 301
points (±0.2 mm at 0.8 nm), step λ_u/16 (10 nm: within 0.6 % of 601 points and λ_u/32).

**Physics that GENESIS leaves out.** GEM-PIC-RZ computes the full fields of the beam, so two effects that
standard GENESIS runs omit are present:
- **Transverse space charge.** At γ = 1000 and 4 kA it widens the beam from 10 to 12 µm over 0.6 m, as the
  envelope equation gives (σ'' = K/(4σ), K = 2I/(I_A γ³)). The nonlinear part raises the normalized emittance
  to 0.11 µm. At γ = 5000 the beam widens by 5 %.
- **Longitudinal space charge at λ_r** (GENESIS `&efield`, off by default).

For a like-for-like comparison GENESIS was also run with `&efield` and with a beam that widens the same way.
That beam diverges from a waist with β = 0.95 m and ε_n = 0.105 µm at 10 nm, and β = 8.7 m and
ε_n = 0.057 µm at 0.8 nm (`fel_sigma_sc.in`). These are the grey curves.

| power gain length (fit over 0.45 … 0.6 m and 2.0 … 2.8 m) | 10 nm | 0.8 nm |
|---|---|---|
| GENESIS | 7.45 cm | 37.4 cm |
| GENESIS + longitudinal space charge | 7.52 cm | 37.7 cm |
| GENESIS + longitudinal space charge + beam widening | **7.92 cm** | **38.4 cm** |
| GEM-PIC-RZ before corrections 7 and 8 (see the main README) | 9.53 cm | — |
| GEM-PIC-RZ with correction 7 (undulator kick averaged over the step) | 8.72 cm | 39.9 cm |
| GEM-PIC-RZ with corrections 7 and 8 (current deposited at the time of the position) | **8.33 cm** | **38.7 cm** |
| … with Δz = λ_u/32 | 8.26 cm | |
| … with Δζ = λ_r/32 | 8.28 cm | |

![10 nm](genesis_seeded/compare_10nm.png) ![0.8 nm](genesis_seeded/compare_0.8nm.png)

- **0.8 nm (transverse space charge small):** GEM-PIC-RZ agrees with GENESIS to 0.6 % in gain length. The
  power at the end is 0.83 of GENESIS with the same physics, and 0.75 of plain GENESIS.
- **10 nm (transverse space charge strong):** the gain length is 4 % longer than in GENESIS with the same
  physics. The remaining difference is in the nonlinear transverse dynamics, which the widening model only
  approximates (the emittance grows during the run in GEM-PIC-RZ, but is fixed in GENESIS).
- **Lethargy:** the seeded start-up (the first ~2 gain lengths) agrees closely in both cases.
- **How the corrections were found:** before the corrections the power first dropped, the sign of a detuning,
  and the gain length depended on Δz to first order. The two corrections removed both.
- **Seed helicity:** the seed must rotate with the electrons. For the undulator field of the code (B_x = B0 cos ψ,
  B_y = B0 sin ψ) that is `helicity = -1`; the other helicity is not amplified.

Reproduce:
```
cd genesis_seeded
genesis4 genesis/fel.in            # and fel_sc.in, fel_sigma_sc.in; the same in genesis_x/
gempic_rz gem.deck numerics.backend=openmp     # 12 min on 2 cores
gempic_rz gem_x.deck numerics.backend=openmp   # 34 min
python3 compare.py 10nm gem_out_nocorr gem_out gem_out_v res_dzeta_v
python3 compare.py 0.8nm gemx_out_nocorr gemx_out   # gemx_out_nocorr: with correction 7 only
```
GENESIS 4 was built from github.com/svenreiche/Genesis-1.3-Version4 (v1.10.10) against serial HDF5. That needed a
small shim header defining `H5Pset_fapl_mpio` and `H5Pset_dxpl_mpio` as no-ops; it was run on one MPI rank.

## LWFA against FBPIC

`lwfa_fbpic/gem.deck` and `lwfa_fbpic/fbpic_lwfa.py` set up the same case:
- **Laser:** a0 = 2, w0 = 8.02 µm, field 1/e half length 3.18 µm (τ = 10.6 fs), λ = 0.8 µm, linear polarization,
  starting at focus.
- **Plasma:** n0 = 3.48·10¹⁸ cm⁻³ (0.002 n_c) behind a 20 µm linear ramp, column radius 19.1 µm, immobile ions.
- **FBPIC grid:** Δz = λ/32, Δr = 0.2 µm, 2 modes, 16 particles per cell, open boundaries.
- **Comparison:** at ct = 50 and 100 µm.

| | GEM-PIC-RZ | FBPIC |
|---|---|---|
| first wake period, min E_z on axis at 100 µm [m c ω0/e] | −0.0339 | −0.0342 |
| rms difference of E_z over the first period | 1.6 % of the peak | |
| laser peak E_x at 100 µm | 1.971 | 1.954 |
| run time on one core | 1.7 min | ~1.5 h |

![LWFA](lwfa_fbpic/compare.png)

- **At 100 µm,** where the plasma behind the laser is uniform, the first period agrees to 1.6 % rms. Behind it,
  where the wave has steepened and is close to breaking, both codes show the same spiky structure.
- **At 50 µm** the rear of the FBPIC window is still in the up-ramp, so its wake there is weaker. In
  GEM-PIC-RZ the plasma of a step is loaded at the front of the window with the density n(s) of that step and
  streams through the whole window within the step. So the profile is a function of the propagation distance
  s and constant across the window. This is how the code samples the plasma profile, not an error; it matters
  only where the profile changes over a window length.
- **Coordinates:** FBPIC's plasma profile is the lab profile of GEM-PIC-RZ, which loads the plasma at the front of
  the window, so it is shifted by the 7.64 µm between the laser centre and the window front.

Reproduce:
```
cd lwfa_fbpic
python3 fbpic_lwfa.py 50,100     # writes fbpic_D50.txt, fbpic_D100.txt (pip install fbpic pyfftw)
gempic_rz gem.deck numerics.backend=openmp
python3 compare.py
```
