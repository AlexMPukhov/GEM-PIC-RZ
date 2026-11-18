# Reading and plotting the output of GEM-PIC-RZ

Everything is in `output.dir` (default `diag`, relative to the directory where the run was started). Every text
file starts with `#` header lines that name its columns: read them instead of assuming.

## The console

```
step    25  z =   0.1273 mm  |Ex|max = 2.031  Ez = [-0.0312, 0.0195]  trapped 0.000 pC (N = 0)
      witness    N = 20000  Q =   5.000 pC  E =    517.93 MeV  dE/E = 1.02e-03  sig_x = 0.371 um  eps_nx = 0.301 um  sig_z = 0.51 um
```

One line per step (largest |E_x| in the window, in m c ω0 / e: the laser amplitude, equal to a0 at the focus; range of
the on-axis E_z; the
charge of trapped plasma electrons) and one per bunch. The run ends with `done: N steps in T s` and
`peak memory: …`. Before the steps the summary of the deck is printed in physical units (also in `summary.txt`);
`WARNING:` lines follow it. The first output step prints `output volume: … GB per output step …`.

## Files written by every run

| file | content |
|---|---|
| `summary.txt`, `deck_used.txt` | the parameters in physical units; the deck after command-line overrides |
| `bunches.txt` | one line per step and bunch: `step  z[mm]  bunch  N  Q[pC]  E[MeV]  dE/E(rms)  sig_x[um]  sig_y[um]  eps_nx[um]  eps_ny[um]  <zeta>[um]  sig_z[um]`. The bunch `trapped` holds the self-injected plasma electrons, `inj_<ion>` the trapped electrons from field ionization |
| `axis_NNNNN.txt` | on the axis, against ζ: `zeta[um]  Ex  Ez  ne/n0  Zni/n0` (NNNNN = step; the header gives z) |
| `map_NNNNN.txt` | the plane y = 0: `zeta[um]  x[um]  Ex  Ez  ne/n0  Zni/n0`, blocks of constant ζ separated by blank lines (gnuplot pm3d), every `output.map_stride`-th slice |
| `beam_NNNNN.txt` | beam particles: `x[um] y[um] zeta[um] ux uy uz q[pC] species` (momenta in m c; species = bunch index, −1 self-trapped) |
| `fel.txt`, `fel_power.txt` | with undulators: radiation energy, peak power, bunching per step; power at λ_r per ζ slice and step (the gain curve) |
| `ionization.txt` | with `ion.ionization = adk`: charge states and freed charge per step |
| `data_NNNNN.h5` | `output.format = hdf5`: all modes on the native grid, densities, beam list, code units with `unitSI` attributes |
| `openpmd/openpmd_<step>.h5` | `output.format = openpmd`: openPMD 1.1 (thetaMode fields E, B, n_e; one particle species per bunch) |
| `checkpoint_NNNNN.h5` | with `output.checkpoint_every` |

`axis`, `map` and `beam` files appear every `output.every` steps and at the start and the end; `bunches.txt` and the
FEL tables have every step. ζ runs from 0 (rear end of the window) to the window length (front).

## Units

- Fields are in m c ω0 / e = 3.2107·10¹² V/m / λ0[µm] (4.013·10¹² V/m at 0.8 µm; B: 1.0709·10⁴ T / λ0[µm]). A laser has
  max |E_x| = a0. The cold wave-breaking field is √(n/n_c) in this unit.
- Densities are relative to `plasma.density`.
- The lab position of a point of the window is z + ζ, where z (in the headers and in `bunches.txt`) is the
  propagation distance of the window.

## Quick look

```
python3 scripts/plot_run.py <output.dir> [--step N]     # table per output step (numpy only); with matplotlib also quicklook_NNNNN.png
gnuplot -e "step=50; dir='<output.dir>'" <repo>/tools/plot_maps.gp   # maps_00050.pdf (rlim, ezlim set the ranges)
python3 <repo>/tools/gempic_gui.py my.deck              # browser GUI: live maps, cuts, bunch table
```

Plain numpy:

```python
import numpy as np
zeta, Ex, Ez, ne = np.loadtxt("diag/axis_00050.txt", usecols=(0, 1, 2, 3), unpack=True)
rows = [l.split() for l in open("diag/bunches.txt") if not l.startswith("#") and l.split()[2] == "witness"]
z_mm, E_MeV, spread = (np.array([float(r[k]) for r in rows]) for k in (1, 5, 6))
```

HDF5 (`<repo>/build/h5_extract` needs no Python; h5py reads the complex modes directly):

```
h5_extract diag/data_00050.h5 info                  # what is in the file
h5_extract diag/data_00050.h5 lineout 2.0 90        # E, B against zeta at r = 2 um, phi = 90 deg
h5_extract diag/data_00050.h5 energy [zeta_lo zeta_hi]      # field energy per mode (m >= 1: radiation)
h5_extract diag/data_00050.h5 spectrum <zeta_lo> <zeta_hi> [k_min k_max nk]
```

openPMD (openPMD-viewer, VisualPIC, ParaView):

```python
from openpmd_viewer import OpenPMDTimeSeries
ts = OpenPMDTimeSeries("diag/openpmd")
Ez, info = ts.get_field("E", "z", iteration=100, theta=0)
z, uz, w = ts.get_particle(["z", "uz", "w"], species="witness", iteration=100)
```

## What to look at

| case | quantity | where | sanity estimate |
|---|---|---|---|
| any laser | amplitude against z | console `\|Ex\|max`, `axis_*.txt` | starts at a0 (less if the focus is ahead: a0/√(1 + z_focus²/z_R²)); in vacuum it follows the Gaussian beam |
| LWFA | wake amplitude | most negative `Ez` in the first bucket of `axis_*.txt` (the positive spike at the bubble closure is not the wake amplitude) | in units of √(n/n_c): upper estimate (a0²/2)/√(1 + a0²/2) for a resonant, wide pulse (1D); 0.84 in `examples/lwfa_plain.deck` (a0 = 2, 30 fs, matched spot); of order √a0 in a fully formed bubble |
| LWFA | laser amplitude and spot against z | console `|Ex|max`, table of `plot_run.py` | in `examples/lwfa_plain.deck` (a0 = 2, matched spot, P ≈ P_c) the amplitude rises to 2.08 at 0.5 mm (self-focusing) and falls to 1.81 at 1 mm (diffraction of the unguided part, depletion). A fall that changes with the box size is the box (field at the absorber) |
| witness | energy, spread, emittance against z | `bunches.txt` | gain = e E_z × distance: E_z[code] × 3.2107 TV/m / λ0[µm] |
| trapped charge | `trapped`, `inj_*` lines | `bunches.txt` | self-injection is **not quantitative** (see physics_rules.md); with a laser set `numerics.trap_gamma` so that quivering electrons are not counted |
| PWFA | accelerating and decelerating field | `axis_*.txt` | of order √(n/n_c) in the blowout regime; transformer ratio ≤ 2 for a symmetric driver |
| FEL | gain curve | `fel_power.txt` (column of a ζ slice inside the bunch's radiation), `fel.txt` | exponential stage: fit ln P against z; 1D gain length λ_u/(4π√3 ρ) |
| undulator emission | radiated energy, spectrum | `h5_extract … energy`, `… spectrum` | N_u-cycle pulse at λ_r ahead of the bunch |

Report numbers with their source file and step. If a quantity is noisy or resolution dependent, say so and show
the comparison (see "Convergence check" in physics_rules.md).
