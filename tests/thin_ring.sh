#!/bin/sh
# Self-fields of a thin, intense, low-gamma ring (numerics.beam_consistent, numerics.theta_stiff; tests/thin_ring.deck):
#  1. static ring at gamma = 1e6 (exact discrete steady state at s = 0): E_r - B_phi stays 0 (the old scheme gave
#     6e-4 of max E_r after one step: the beam j_z of the E_z equation did not match the slice charge of R);
#  2. ring expanding with prescribed v_r = 3e-4 c (gamma = 1e6, no forces): E_r - B_phi = O(v_r^2) exactly; the old
#     scheme gave 1.7e-3 of max E_r (undamped odd-even mode of the stiff E_z constraint, non-conserving j_r);
#  3. 10 kA ring at rest, gamma 50, no lens: r percentiles after 0.5 mm against the shell model (tests/thin_ring.py);
#  4. the same ring rotating on matched orbits in a plasma lens: radius, u_r rms and <r u_phi> against the shell
#     model (old scheme: u_r rms 0.32, <r u_phi> +137 %). The shell model has no induced E_phi, so its <r u_phi> is
#     exactly constant; the code keeps it to 0.3 % (limit 0.5 %).
# Beyond ~0.5 mm the longitudinal space charge of the flat-top (head and tail +-10 % in gamma per mm) mismatches
# head and tail to the lens, which the shell model (no E_z) does not describe.
# usage: sh tests/thin_ring.sh <gempic_rz> [backend]     (exit 77: skipped, needs python3 with numpy and h5py)
set -e
abspath() { (cd "$(dirname "$1")" && echo "$(pwd)/$(basename "$1")"); }
EXE=$(abspath "$1")
BACKEND=${2:-serial}
SRC=$(cd "$(dirname "$0")/.." && pwd)
python3 -c "import numpy, h5py" 2>/dev/null || { echo "thin-ring test skipped: needs python3 with numpy and h5py"; exit 77; }
WORK=$(mktemp -d)
cd "$WORK"
DECK="$SRC/tests/thin_ring.deck"
export OMP_PROC_BIND=false
REST="plasma_lens.gradient=0T/m bunch.ux_function=0 bunch.uy_function=0"
H5="output.format=hdf5 output.particles=false output.every=1 grid.z_end=60um"

echo "--- 1. static ring, gamma = 1e6"
"$EXE" "$DECK" numerics.backend=$BACKEND output.dir=static bunch.gamma=1e6 bunch.charge=0.001pC $REST $H5 > static.log
python3 "$SRC/tests/thin_ring.py" static static 1e-9 1e-9

echo "--- 2. prescribed expansion v_r = 3e-4 c, gamma = 1e6"
"$EXE" "$DECK" numerics.backend=$BACKEND output.dir=expand bunch.gamma=1e6 bunch.charge=0.001pC \
  plasma_lens.gradient=0T/m "bunch.ux_function=300 * x / Rw" "bunch.uy_function=300 * y / Rw" $H5 > expand.log
python3 "$SRC/tests/thin_ring.py" static expand 2e-5 1e-4

echo "--- 3. 10 kA ring at rest, gamma 50"
"$EXE" "$DECK" numerics.backend=$BACKEND output.dir=rest $REST > rest.log
python3 "$SRC/tests/thin_ring.py" ring rest 33.4 1 0 0.5 0.06 1 1

echo "--- 4. 10 kA ring on matched orbits in a plasma lens"
"$EXE" "$DECK" numerics.backend=$BACKEND output.dir=lens > lens.log
python3 "$SRC/tests/thin_ring.py" ring lens 33.4 1 3.0e7 0.5 0.03 0.01 5e-3

cd - > /dev/null
rm -rf "$WORK"
echo "thin-ring test passed"
