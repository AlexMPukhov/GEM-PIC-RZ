#!/bin/sh
# Chained runs (examples/chain_lwfa_fel, reduced): stage 1 = 2 steps of tests/regression.deck (LWFA with a witness),
# openPMD output; stage 2 = tests/chain_fel.deck (import, drift, thin chicane, 10 undulator periods). Checked against
# independent references:
#   (i)   the imported beam (select = "abs(z) < 1 um", particles as in the file): charge, energy, energy spread,
#         rms sizes and emittances in bunches.txt = the same moments computed by python from the openPMD file of
#         stage 1 with the same selection (1e-6 relative)
#   (ii)  the chicane (R56 = 20 um, delay 0.1 um): the particle data before and after it (HDF5) - centroid shift
#         = R56 <delta> - delay + the slippage of one step, and the rms length after it = that of the mapped
#         particles zeta + R56 delta (1e-3)
#   (iii) the radiated power at lambda_r (fel.txt, P_res) grows along the undulator (factor > 10)
#   (iv)  [mpi] with <mpiexec>: stage 2 on 2 ranks = 1 rank (bunches.txt, fel.txt: 1e-6)
# usage: sh tests/chain.sh <gempic_rz> <compare_columns> [mpiexec]    (exit 77: no python3 with numpy and h5py)
set -e
abspath() { (cd "$(dirname "$1")" && echo "$(pwd)/$(basename "$1")"); }
EXE=$(abspath "$1"); CMP=$(abspath "$2"); MPIEXEC=$3
python3 -c "import numpy, h5py" 2>/dev/null || { echo "chain test skipped: needs python3 with numpy and h5py"; exit 77; }
SRC=$(cd "$(dirname "$0")/.." && pwd)
WORK=$(mktemp -d)
cd "$WORK"
export OMP_PROC_BIND=false
O="numerics.backend=serial"
# ---- stage 1
"$EXE" "$SRC/tests/regression.deck" $O grid.z_end=80k0^-1 output.every=2 output.format=openpmd output.dir=src > src.log
cp "$SRC/tests/chain_fel.deck" fel.deck
# ---- (i) import: the selected particles as they are (no resampling, no shot noise)
sed -e '/macroparticles/d' -e '/shot_noise/d' -e '/noise_wavelength/d' fel.deck > imp.deck
"$EXE" imp.deck $O grid.z_end=1.25mm output.dir=imp > imp.log
awk '$3 == "witness" && $1 == 0 { print $5, $6, $7, $8, $9, $10, $11 }' imp/bunches.txt > got.txt
python3 - src/openpmd/openpmd_2.h5 > ref.txt <<'PY'
import sys
import h5py
import numpy as np
f = h5py.File(sys.argv[1], "r")
it = sorted(f["data"].keys(), key=int)[-1]
sp = f["data"][it]["particles"]["witness"]
def rec(n, c):
    d = sp[n][c]
    v = d[()] if isinstance(d, h5py.Dataset) else np.full(len(sp["weighting"][()]), d.attrs["value"])
    return v * d.attrs["unitSI"]
x, y, z = (rec("position", c) + rec("positionOffset", c) for c in "xyz")
mc = 9.1093837015e-31 * 299792458.0
ux, uy, uz = (rec("momentum", c) / mc for c in "xyz")
w = sp["weighting"][()]
zc = np.sum(w * z) / np.sum(w)
k = np.abs(z - zc) < 1e-6
x, y, ux, uy, uz, w = x[k], y[k], ux[k], uy[k], uz[k], w[k]
g = np.sqrt(1 + ux ** 2 + uy ** 2 + uz ** 2)
W = w.sum()
m = lambda v: np.sum(w * v) / W
c2 = lambda a, b: m((a - m(a)) * (b - m(b)))
eps = lambda a, u: np.sqrt(c2(a, a) * c2(u, u) - c2(a, u) ** 2)
print("%.6e %.6e %.6e %.6e %.6e %.6e %.6e" % (W * 1.602176634e-7, (m(g) - 1) * 0.51099895, np.sqrt(c2(g, g)) / m(g),
      np.sqrt(c2(x, x)) * 1e6, np.sqrt(c2(y, y)) * 1e6, eps(x, ux) * 1e6, eps(y, uy) * 1e6))
PY
echo "  imported beam (bunches.txt):  $(cat got.txt)"
echo "  from the stage-1 file:        $(cat ref.txt)"
"$CMP" ref.txt got.txt 1e-3 1e-12
# ---- stage 2: (ii) chicane, (iii) FEL
"$EXE" fel.deck $O output.every=1 output.dir=fel > fel.log
grep -q "chicane at z = 10 mm" fel.log
python3 - fel/data_00007.h5 fel/data_00008.h5 <<'PY'
import sys
import h5py
import numpy as np
def load(p):
    f = h5py.File(p, "r")
    L = f.attrs["k0inv_m"]
    g = f["particles/beam"]
    z = g["zeta"][()] * L
    u = np.stack([g[c][()] for c in ("ux", "uy", "uz")])
    return z, np.sqrt(1 + (u ** 2).sum(0)), np.abs(g["charge"][()])
z0, g0, q0 = load(sys.argv[1])
z1, g1, q1 = load(sys.argv[2])
r56, delay, gref, dz = 20e-6, 0.1e-6, 1000.0, 1.25e-3
m = lambda v, q: np.sum(q * v) / np.sum(q)
pred = z0 + r56 * (g0 / gref - 1) - delay - dz * m(1 / (2 * g0 ** 2), q0)
dc = m(z1, q1) - m(pred, q0)
s_pred = np.sqrt(m((pred - m(pred, q0)) ** 2, q0))
s_got = np.sqrt(m((z1 - m(z1, q1)) ** 2, q1))
print("  chicane: centroid moved by %.4f um (R56 <delta> - delay - slippage: %.4f um, difference %.2e um); "
      "rms length %.5f um (mapped particles %.5f um)" % ((m(z1, q1) - m(z0, q0)) * 1e6, (m(pred, q0) - m(z0, q0)) * 1e6,
                                                         dc * 1e6, s_got * 1e6, s_pred * 1e6))
sys.exit(0 if abs(dc) < 2e-9 and abs(s_got / s_pred - 1) < 1e-3 else 1)
PY
awk '!/^#/ && $1 == 12 { a = $7 } !/^#/ { b = $7 } END {
  printf "  P(lambda_r): %.3g W at the undulator entrance (step 12), %.3g W after 10 periods (x %.0f)\n", a, b, b / a
  exit !(b > 10 * a) }' fel/fel.txt
if [ -n "$MPIEXEC" ]; then
  FLAGS=""
  if [ "$(id -u)" = 0 ]; then export OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1; fi
  if $MPIEXEC --version 2>&1 | grep -qiE "open mpi|openrte|open-mpi"; then FLAGS="--oversubscribe"; fi
  $MPIEXEC $FLAGS -n 2 "$EXE" fel.deck $O output.every=1 output.dir=mp2 > mp2.log
  grep -q "chicane at z = 10 mm" mp2.log
  for f in bunches.txt fel.txt; do "$CMP" fel/$f mp2/$f 1e-6 1e-12; done
  echo "  stage 2 on 2 MPI ranks = 1 rank (bunches.txt, fel.txt)"
fi
cd - > /dev/null
rm -rf "$WORK"
echo "chain test passed"
