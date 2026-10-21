# Flying focus (luminal Bessel pulse fed through the radial boundary): first tests

Block `flying_focus` (docs/input_deck.md). The driver is a converging conical wave whose on-axis peak moves at c.
It is fed in through the radial absorber layer: R relaxes there to the incident field, and only what the plasma
adds is absorbed. This is the total-field / scattered-field method.

All three cases: λ = 0.8 µm, a0 = 2 on the axis, 50 fs, cone angle 50 mrad (core radius 6.1 µm, transverse
wavelength 16 µm), r_max = 60–80 µm, dr = 0.5 µm, Δζ = λ/24, 2 CPU cores.

| deck | content | result |
|---|---|---|
| `vacuum.deck` | vacuum, 2 mm (`absorber_strength = 1`) | on-axis amplitude within 2 % of the incident field; carrier phase offset settles at 0.08 rad (numerical dispersion of the conical wave, ~1 % in κ); envelope moves at c to 0.07 µm over 2 mm. Serial, OpenMP and 2 MPI ranks are identical. With the default `absorber_strength = 0.2` the amplitude fell by 9 % |
| `gas.deck` | neutral H with an **artificial ionization potential of 200 eV** (ADK: full ionization at a ≳ 0.8), gas column r < 55 µm, layer 60–80 µm plasma-free; 1 mm, 15 min | channel ionized to r ≈ 10–20 µm (⟨Z⟩ 0.03 over the column); laser on the axis within 3 % of vacuum over 1 mm (no depletion); first wake period stationary in ζ. Behind it the wake is irregular: electrons born near the peak of a 200 eV potential are hot (ionization heating), an artefact of the artificial potential |
| `plasma.deck` | preformed plasma 1e18 cm⁻³ in the whole box, incident field computed in that plasma (`background_density`); 1 mm, 40 min | first wake period (ζ = 60–110 µm) stationary from 0.4 mm on (zero crossing 60.3–61.4 µm); on-axis peak 2.0 → 1.8 a0 in the first 0.5 mm, then constant; behind the first period the wakes of the Bessel rings cross and the wake is irregular |

Figures: `gas_out.png`, `plasma_out.png` (`python3 plot_gas.py <dir>`); numbers: `python3 check_wake.py <dir>`,
`python3 check_axis.py <dir> 2 0.05 50 90`.

**Why real hydrogen needs a different boundary.** The incident intensity falls off as I(r)/I_axis ≈ 2/(π k_⊥ r).
- Real hydrogen (13.6 eV) is ionized by a 50 fs pulse already at a ≈ 0.005, about 5·10¹³ W/cm².
- At a0 = 2 a neutral boundary would therefore have to sit at r ≈ 0.25 m, which is impossible.
- At that low intensity the ionization happens early in the pulse front and hardly heats the electrons, so
  real hydrogen behaves as a preformed plasma. The boundary then lies in the plasma, and the incident field
  must be that of the plasma (`plasma.deck`).
- A neutral-gas boundary is only possible for high ionization potentials (inner shells, or the 200 eV test gas).

**Production suggestions** (cluster):
- **Physical hydrogen:** use the setup of `plasma.deck`, with a larger r_max (several transverse wavelengths
  between the ionized core and the layer), a finer Δζ (λ/32) and a longer propagation.
- **Ionization in the core:** keep a high-Ip species for the core physics; the gas test shows that the boundary
  works.
- **Convergence checks:** dr, Δζ and Δz (the carrier phase offset in vacuum is a numerical-dispersion effect).
- **Other drivers:** the incident field is the ideal luminal focus-wave mode (infinite energy). An axiparabola /
  echelon flying focus of finite length needs a different R_inc (e.g. a Hankel propagation of the field after
  the optic); the boundary method stays the same.
