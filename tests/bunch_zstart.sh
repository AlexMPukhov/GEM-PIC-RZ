#!/bin/sh
# bunch.z_start: a bunch loaded when the run reaches z_start, with its co-moving self-field. Vacuum: a rotating
# ring (gamma 50, x offset: modes 0 and 1; beam_consistent: also the field of its azimuthal current) loaded at
# z = 0.1 mm. Reference: the same deck with the ring at z = 0 (the random numbers are drawn in the same order).
# Checked: no particles and zero fields before z_start; at the loading the state (fields, particles) equals the
# reference start, and 5 steps later the reference after 5 steps (to rounding: the window does not know z).
# usage: sh tests/bunch_zstart.sh <gempic_rz> [backend]     (exit 77: skipped, needs python3 with numpy and h5py)
set -e
abspath() { (cd "$(dirname "$1")" && echo "$(pwd)/$(basename "$1")"); }
EXE=$(abspath "$1")
BACKEND=${2:-serial}
python3 -c "import numpy, h5py" 2>/dev/null || { echo "z_start test skipped: needs python3 with numpy and h5py"; exit 77; }
WORK=$(mktemp -d)
cd "$WORK"
export OMP_PROC_BIND=false
cat > zs.deck <<'DECK'
constants { Rw = 1 um }
grid {
  lambda0 = 0.8 um; m_max = 1; radial_grid = uniform; dr = 0.1 um; r_max = 30 um
  dzeta = 0.1 um; window_length = 20 um; dz = 20 um; z_end = 0.2 mm; r_absorber = 25 um
}
numerics { beam_field = mean; beam_dz_max = 5 um; beam_consistent = true; theta_stiff = 1 }
bunch { name = ring; charge = 3 pC; gamma = 50; distribution_z = flattop; length = 1 um; zeta_center = 6 um
        distribution_r = function; profile_r = "gauss(r - Rw, 50 nm)"; r_extent = 1.5 um; x0 = 0.2 um
        ux_function = "-0.5 * y / Rw"; uy_function = "0.5 * x / Rw"; macroparticles = 20000; z_start = 0.1 mm }
output { dir = late; every = 1; maps = false; particles = true; format = hdf5 }
DECK
"$EXE" zs.deck numerics.backend=$BACKEND > late.log
grep -q "bunch 'ring' loaded at z = 0.1 mm" late.log || { echo "FAIL: no loading message"; cat late.log; exit 1; }
"$EXE" zs.deck numerics.backend=$BACKEND bunch.z_start=0mm grid.z_end=0.1mm output.dir=ref > ref.log
python3 - <<'PY'
import sys
import h5py
import numpy as np
ok = True
def parts(f):  # sorted by (zeta, x, y): the GPU appends particles in a nondeterministic order
    g = f["particles/beam"]
    p = np.stack([g[c][()] for c in ("x", "y", "zeta", "ux", "uy", "uz")], 1)
    return p[np.lexsort((p[:, 1], p[:, 0], p[:, 2]))]
with h5py.File("late/data_00004.h5") as A:  # z = 0.08 mm: before the loading
    n = len(A["particles/beam/x"][()])
    f = max(abs(A["fields"][k][()]).max() for k in A["fields"])
    print("before z_start: %d particles, max |field| %.1e" % (n, f))
    ok = ok and n == 0 and f == 0.0
for la, re in (("late/data_00005.h5", "ref/data_00000.h5"), ("late/data_00010.h5", "ref/data_00005.h5")):
    with h5py.File(la) as A, h5py.File(re) as B:
        worst = max(abs(A["fields"][k][()] - B["fields"][k][()]).max() / abs(B["fields"]["Er"][()]).max()
                    for k in B["fields"])
        pa, pb = parts(A), parts(B)
        dp = abs(pa - pb).max() / abs(pb).max() if pa.shape == pb.shape else 1.0
        print("%s vs %s: fields %.2e, particles (%d) %.2e (limits 1e-10)" % (la, re, worst, len(pa), dp))
        ok = ok and worst < 1e-10 and dp < 1e-10
sys.exit(0 if ok else 1)
PY
cd - > /dev/null
rm -rf "$WORK"
echo "z_start test passed"
