#!/usr/bin/env python3
"""Quick look at a GEM-PIC-RZ run from its text output (output.format = text, the default).

  python3 plot_run.py diag                # newest output step of the directory `diag`
  python3 plot_run.py diag --step 50

Prints a table for every output step (laser amplitude, wake field on the axis, field at the radial absorber, density
spike) and the last line of every bunch; this part needs numpy only. With matplotlib it also writes
<dir>/quicklook_NNNNN.png: on-axis E_x (laser) and E_z (wake), the electron density on the axis, the map of E_z in
the plane y = 0 with the beam particles, and the bunch table against z.
Fields are in the code unit m c omega0 / e (E_x of a laser = a0); multiply by 3.2107e12 V/m / lambda0[um] for SI.
"""
import argparse
import glob
import os
import re
import numpy as np

p = argparse.ArgumentParser()
p.add_argument("dir")
p.add_argument("--step", type=int)
a = p.parse_args()

steps = sorted(int(os.path.basename(f)[5:10]) for f in glob.glob(os.path.join(a.dir, "axis_*.txt")))
if not steps:
    raise SystemExit(f"no axis_NNNNN.txt in {a.dir}: is output.format = text and output.axis = true?")
step = a.step if a.step is not None else steps[-1]
name = lambda kind, s=None: os.path.join(a.dir, f"{kind}_{step if s is None else s:05d}.txt")


def load_map(s):
    """map file: zeta[um]  x[um]  Ex  Ez  ne/n0  Zni/n0, blocks of constant zeta -> arrays [zeta][x]"""
    m = np.loadtxt(name("map", s))
    nz = len(np.unique(m[:, 0]))
    return [m[:, k].reshape(nz, -1) for k in range(4)]


# ---- absorber radius and wavelength from summary.txt, if present
r_abs = lam0 = None
try:
    txt = open(os.path.join(a.dir, "summary.txt")).read()
    m = re.search(r"absorber r > ([0-9.eE+-]+) um", txt)
    r_abs = float(m.group(1)) if m else None
    m = re.search(r"lambda0 = ([0-9.eE+-]+) um", txt)
    lam0 = float(m.group(1)) if m else None
except OSError:
    pass

# ---- text table: one line per output step
unit = f"1 code unit = {3.2107e3 / lam0:.4g} GV/m" if lam0 else "code units m c omega0 / e"
print(f"# {a.dir}: fields in m c omega0 / e ({unit})")
print("# step   z[mm]   max|Ex| axis   min Ez axis   max Ez axis   max ne/n0 axis" +
      ("   max|Ex| at r_absorber / peak" if r_abs else ""))
for s in steps:
    z, ex, ez, ne = np.loadtxt(name("axis", s), usecols=(0, 1, 2, 3), unpack=True)
    with open(name("axis", s)) as f:
        m = re.search(r"z = ([0-9.eE+-]+) mm", f.readline())
    line = f"{s:6d}  {float(m.group(1)) if m else float('nan'):7.4f}  {np.abs(ex).max():12.4f}  {ez.min():12.5f}  " \
           f"{ez.max():12.5f}  {ne.max():15.2f}"
    if r_abs and os.path.exists(name("map", s)):
        Z, X, EX, EZ = load_map(s)
        j = np.argmin(np.abs(X[0] - r_abs))
        peak = np.abs(EX).max()
        line += f"  {np.abs(EX[:, j]).max() / peak if peak > 0 else 0:20.4f}"
    print(line)
if r_abs:
    print("# the last column should stay below ~0.01: otherwise the absorber (and the plasma edge) cut the laser; "
          "enlarge the box")

# ---- bunches.txt: step z[mm] bunch N Q[pC] E[MeV] dE/E sig_x sig_y eps_nx eps_ny <zeta> sig_z
rows = {}
with open(os.path.join(a.dir, "bunches.txt")) as f:
    for line in f:
        if line.startswith("#") or not line.strip():
            continue
        w = line.split()
        rows.setdefault(w[2], []).append([float(x) for x in w[:2] + w[3:]])
print("# bunch: last step  z[mm]  Q[pC]  E[MeV]  dE/E  sig_x[um]  eps_nx[um]  sig_z[um]   (first step: E[MeV], Q[pC])")
for bunch, r in rows.items():
    r = np.array(r)
    print(f"{bunch:>10s} {int(r[-1, 0]):6d}  {r[-1, 1]:.4f}  {r[-1, 3]:.4g}  {r[-1, 4]:.6g}  {r[-1, 5]:.3g}  "
          f"{r[-1, 6]:.4g}  {r[-1, 8]:.4g}  {r[-1, 11]:.4g}   ({r[0, 4]:.6g}, {r[0, 3]:.4g})")

# ---- figure (optional)
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError:
    raise SystemExit("no figure: matplotlib is not installed (the table above is complete; pip install matplotlib)")

fig, ax = plt.subplots(2, 2, figsize=(12, 8))
z, ex, ez, ne = np.loadtxt(name("axis"), usecols=(0, 1, 2, 3), unpack=True)
with open(name("axis")) as f:
    head = f.readline().strip("# \n")
ax[0, 0].plot(z, ex, lw=0.5, color="C0")
ax[0, 0].set(xlabel="zeta [um] (window; front at the right)", ylabel="E_x on the axis (laser) [m c omega0 / e]",
             title=head[head.find("("):] if "(" in head else head)
axz = ax[0, 0].twinx()                      # the wake field is much smaller than the laser field: own scale
axz.plot(z, ez, color="C1")
axz.set_ylabel("E_z on the axis (wake)", color="C1")
ax[0, 1].plot(z, ne)
ax[0, 1].set(xlabel="zeta [um]", ylabel="n_e / n0 on the axis")

if os.path.exists(name("map")):
    Z, X, EX, EZ = load_map(step)
    lim = np.percentile(np.abs(EZ), 99.5) or 1.0
    im = ax[1, 0].pcolormesh(Z, X, EZ, cmap="RdBu_r", vmin=-lim, vmax=lim, shading="auto")
    fig.colorbar(im, ax=ax[1, 0], label="E_z")
    if os.path.exists(name("beam")) and os.path.getsize(name("beam")) > 0:
        b = np.atleast_2d(np.loadtxt(name("beam")))
        if b.size:
            ax[1, 0].plot(b[:, 2], b[:, 0], "k,", alpha=0.5)
    ax[1, 0].set(xlabel="zeta [um]", ylabel="x [um]", title="E_z in the plane y = 0, beam particles black")

twin = ax[1, 1].twinx()
shown = 0
for bunch, r in rows.items():
    r = np.array(r)
    if r[:, 3].max() <= 0:
        continue
    shown += 1
    ax[1, 1].plot(r[:, 1], r[:, 4], label=f"{bunch}: E [MeV]")
    twin.plot(r[:, 1], r[:, 3], "--", label=f"{bunch}: Q [pC]")
ax[1, 1].set(xlabel="z [mm]", ylabel="mean energy [MeV]", title="bunches.txt (solid: energy, dashed: charge)")
twin.set_ylabel("charge [pC]")
if shown:
    ax[1, 1].legend(loc="upper left", fontsize=8)
else:
    ax[1, 1].text(0.5, 0.5, "no bunch and no trapped charge", ha="center", transform=ax[1, 1].transAxes)
fig.tight_layout()
out = os.path.join(a.dir, f"quicklook_{step:05d}.png")
fig.savefig(out, dpi=130)
print(out)
