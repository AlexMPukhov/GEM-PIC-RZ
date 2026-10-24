"""Energy spread dE/E of bunches.txt against the particle data of the same run (tests/energy_spread.sh).

The reference is computed from /particles/beam of every data_NNNNN.h5: charge-weighted mean and variance of gamma
in two passes and in extended precision. The table prints 4 digits: the printed value must be the reference rounded (within
half a unit of the last digit), and dE/E = 0 for a monoenergetic bunch.
usage: python3 energy_spread.py <output dir> <bunch name>"""
import glob
import os
import sys

import h5py
import numpy as np

out, name = sys.argv[1], sys.argv[2]
table = {}
for line in open(os.path.join(out, "bunches.txt")):
    if line.startswith("#"):
        continue
    f = line.split()
    if f[2] == name:
        table[int(f[0])] = float(f[6])
worst, n = 0.0, 0
for fn in sorted(glob.glob(os.path.join(out, "data_*.h5"))):
    step = int(fn[-8:-3])
    with h5py.File(fn) as h:
        g = h["/particles/beam"]
        names = list(g["species"].attrs["bunch_names"]) if "bunch_names" in g["species"].attrs else []
        idx = [i for i, s in enumerate(names) if (s.decode() if isinstance(s, bytes) else str(s)) == name]
        sel = g["species"][()] == (idx[0] if idx else 0)
        u = [np.asarray(g[c][()][sel], dtype=np.longdouble) for c in ("ux", "uy", "uz")]
        w = np.abs(np.asarray(g["charge"][()][sel], dtype=np.longdouble))
    gam = np.sqrt(1 + u[0] ** 2 + u[1] ** 2 + u[2] ** 2)
    mean = (w * gam).sum() / w.sum()
    ref = float(np.sqrt((w * (gam - mean) ** 2).sum() / w.sum()) / mean)
    got = table[step]
    if ref < 1e-12:  # monoenergetic (gamma from the stored momenta: rounding ~1e-19)
        ok = got == 0.0
        print(f"  step {step}: dE/E = {got:.4e}, reference {ref:.1e} (monoenergetic)")
    else:  # the table prints 4 digits: the printed value must be the reference rounded
        half = 0.5 * 10.0 ** (np.floor(np.log10(ref)) - 3)
        rel = abs(got - ref) / half
        worst = max(worst, rel)
        ok = rel <= 1.0 + 1e-6
        print(f"  step {step}: dE/E = {got:.3e}, reference {ref:.6e}, |difference| = {rel:.2f} half-units of the last digit")
    n += 1
    if not ok:
        print("  FAILED")
        sys.exit(1)
if n < 3:
    print("  too few output steps")
    sys.exit(1)
print(f"  ok: {n} steps, worst |difference| {worst:.2f} half-units of the last digit")
