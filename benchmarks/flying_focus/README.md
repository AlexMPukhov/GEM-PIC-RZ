# Flying focus (luminal Bessel pulse fed through the radial boundary): first tests

Handover, open tasks and lessons: `docs/flying_focus_handoff.md`. Production decks: `production/` (start with
`ff_plasma.deck`; every line is commented).

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

**Several components (beats), `beat_vacuum.deck`** (ctest `flying_focus_beat`): J0 at 40 mrad + 0.2 J1 e^{iφ} at
100 mrad, circularly polarized (core in mode 1, J1 part in mode 0), a0 = 1, 30 fs, dr = 0.4 µm, Δz = 5 µm, 0.12 mm.
The intensity maximum near the axis turns around the axis on a circle of 1.15 µm with the beat period
2π/(κ₁ − κ₀) = 190.7 µm: a helically winding driver. Against an independent quadrature of the analytic field
(`check_physics ffbeat`, E_x in the plane y = 0, |ζ − ζc| < 2L, r < 40 µm):

| z | rms(E_x − ref)/a0 | J1 part: rms diff / its rms | intensity maximum x (ref) |
|---|---|---|---|
| 40 µm | 6.0·10⁻⁴ | 3.0 % | 0.575 µm (0.551) |
| 80 µm | 1.0·10⁻³ | 5.4 % | −1.017 µm (−1.032) |
| 120 µm | 1.5·10⁻³ | 7.8 % | −1.024 µm (−0.985) |

- The reference is taken at s − Δz/2: the field fed through the layer at level s equals, after the step, the incident
  field of s − Δz/2 (a half-step lag of the feeding; it also exists for a single component and is part of the
  0.08 rad carrier phase offset of `vacuum.deck`). Against the level s itself the J1 part is off by its phase κ₁Δz/2
  (12 % at 40 µm with Δz = 5 µm, 6 % with 2.5 µm); against s − Δz/2 the difference no longer depends on Δz.
- The growth of the J1 difference with z (~6·10⁻⁴ per µm, independent of Δz) is numerical dispersion of the
  100 mrad cone (~1.5 % in κ₁: radial (k⊥dr)²/12 and ζ discretization). The beat period is therefore right to ~2 %.
- Helicity −1 with Bessel order −1 (again mode 0) gives the same E_x; linear polarization (J1 part in modes 0 and 2,
  m_max = 2) agrees to the same accuracy. Real GPU (RTX 5090): the same numbers.
- Steep cones: at 150 mrad and Δz = 5 µm a single J0 or J1 component differs from the analytic pulse at the level s
  by ~25 % of its rms after 20 µm (half of it at Δz = 2.5 µm); the half-step lag (phase κΔz/2 = 0.22 rad) accounts
  for ~22 % of it. For a beat the lag only turns the pattern by the constant angle (κ₁ − κ₀)Δz/2.

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

**ADK estimate used for the 200 eV test gas** (probability of ionization by a 50 fs, 0.8 um pulse of peak a; hydrogen-like,
static ADK averaged over the carrier):
```
import numpy as np; from scipy.special import gamma as G
def w(F, Ip):  # atomic units, Z = 1, l = m = 0
    k = np.sqrt(2 * Ip); n = 1 / k; C2 = 2**(2 * n) / (n * G(n + 1) * G(n))
    return C2 * Ip * (2 * k**3 / F)**(2 * n - 1) * np.exp(-2 * k**3 / (3 * F))
t = np.linspace(-150e-15, 150e-15, 200001); om = 2 * np.pi * 3e8 / 0.8e-6; Lt = 50e-15 / np.sqrt(2 * np.log(2))
P = lambda a, IpeV: 1 - np.exp(-np.sum(w(np.maximum(a * 7.8 * np.exp(-(t / Lt)**2) * abs(np.cos(om * t)), 1e-9),
                                          IpeV / 27.211)) * (t[1] - t[0]) / 2.4189e-17)
```
(a = 1 at 0.8 um is 7.8 atomic units of field.) Then the radius where the incident wave reaches a: r ~ 2 a0^2 / (pi k_perp a^2).
