"""Coherent spontaneous emission of a short bunch in a planar undulator (emission.deck): GEM-PIC-RZ against the
radiation integral (radint.py), energy per azimuthal mode m = 0 ... 3 (odd m: fundamental and odd harmonics, even m:
even harmonics, zero on the axis) and per harmonic band of the odd and even modes.
The m = 0 spectrum of the last step contains the Coulomb field of the bunch; that of step 0 (same bunch, no
radiation yet) is subtracted.
usage: python3 compare_emission.py <out_dir> [<out_dir> ...]   (needs h5_extract in ../../build)"""
import os
import subprocess
import sys
import numpy as np

H5X = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "build", "h5_extract")
R = np.loadtxt("radint_em_modes.txt")       # k/k0, dW/dk of |m| = 0..7
nu, dn = R[:, 0], R[1, 0] - R[0, 0]


def spectrum(f):
    cache = f + ".spec.txt"
    if not os.path.exists(cache):
        out = subprocess.run([H5X, f, "spectrum", "0", "1e9", "0.3", "3.6", "661"], capture_output=True, text=True,
                             check=True).stdout
        open(cache, "w").write(out)
    return np.loadtxt(cache)


def table(k, modes):
    """modes: dW/dk of m = 0..3 -> energies per mode, then odd/even per band"""
    dk = k[1] - k[0]
    r = [modes[:, m].sum() * dk for m in range(4)]
    for h in (1, 2, 3):
        b = (k > h - 0.5) & (k < h + 0.5)
        r += [(modes[b, 1] + modes[b, 3]).sum() * dk, (modes[b, 0] + modes[b, 2]).sum() * dk]
    return np.array(r)


ref = table(nu, R[:, 1:5])
names = ["m = 0", "m = 1", "m = 2", "m = 3"] + [f"{p} m, k {h - 0.5:.1f}-{h + 0.5:.1f}" for h in (1, 2, 3)
                                                for p in ("odd", "even")]
print(f"{'':18s} {'radint':>10s}" + "".join(f" {d:>18s}" for d in sys.argv[1:]))
print(f"{'(m >= 4 in radint:':18s} {R[:, 5:].sum() * dn:10.4e})")
res = []
for d in sys.argv[1:]:
    steps = sorted(f for f in os.listdir(d) if f.startswith("data_") and f.endswith(".h5"))
    S0, S1 = spectrum(os.path.join(d, steps[0])), spectrum(os.path.join(d, steps[-1]))
    modes = S1[:, 3:7].copy()
    modes[:, 0] -= S0[:, 3]
    res.append(table(S1[:, 0], modes))
for i, n in enumerate(names):
    print(f"{n:18s} {ref[i]:10.4e}" + "".join(f" {r[i]:10.4e} ({r[i] / ref[i] - 1:+5.1%})" for r in res))
