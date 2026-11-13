#!/bin/sh
# Dielectric (block 'dielectric'): constant eps(r, z) in the field equations, against independent references
# (tests/dielectric.py):
#  1. eps = 1 through the dielectric slice solver = the vacuum solver (laser in a metal tube, tests/wall.deck);
#  2. dielectric wakefield accelerator (tests/dielectric_dwa.deck: vacuum channel 50 um, eps = 3.8 up to the metal
#     wall at 75 um, 100 pC, 10 um bunch at v = c, offset by 10 um): on-axis E_z behind the bunch against the
#     modal theory of the synchronous TM_0n modes (wavenumbers 0.5 %, amplitudes from the loss factors 0.3 / 1.7 /
#     0.2 %), the end of the filled TM01 wake against the group velocity of the GEM system (0.3 %; Maxwell's is
#     slower: the system drops d_s L, exact for the synchronous field only), and the on-axis E_x (m = 1) against the
#     synchronous hybrid modes HEM_11, HEM_12 (field matching of E_z, B_z, E_phi, B_phi: 0.6 / 0.3 %);
#  3. Cherenkov radiation of a v = c bunch in an open eps = 2 medium: the cone zeta = zeta_b - r sqrt(eps - 1)
#     (slope 0.8 %);
#  4. a laser in a uniform eps = 1.1 medium: the slippage (n_g - 1) z against the dispersion of the scheme
#     (2e-5; GEM continuum and Maxwell printed).
# usage: sh tests/dielectric.sh <gempic_rz> [backend]     (exit 77: skipped, needs python3 with numpy)
set -e
abspath() { (cd "$(dirname "$1")" && echo "$(pwd)/$(basename "$1")"); }
EXE=$(abspath "$1")
BACKEND=${2:-serial}
SRC=$(cd "$(dirname "$0")/.." && pwd)
python3 -c "import numpy" 2>/dev/null || { echo "dielectric test skipped: needs python3 with numpy"; exit 77; }
WORK=$(mktemp -d)
cd "$WORK"
export OMP_PROC_BIND=false
PY="python3 $SRC/tests/dielectric.py"

echo "--- 1. eps = 1: dielectric solver = vacuum solver"
"$EXE" "$SRC/tests/wall.deck" numerics.backend=$BACKEND output.dir=vac grid.z_end=0.2mm output.every=40 > vac.log
cp "$SRC/tests/wall.deck" one.deck
echo 'dielectric { epsilon = 1 }' >> one.deck
"$EXE" one.deck numerics.backend=$BACKEND output.dir=one grid.z_end=0.2mm output.every=40 > one.log
cmp vac/map_00040.txt one/map_00040.txt && echo "  map after 40 steps identical"

echo "--- 2. dielectric-lined tube (DWA), offset bunch"
"$EXE" "$SRC/tests/dielectric_dwa.deck" numerics.backend=$BACKEND grid.m_max=1 bunch.x0=10um > dwa.log
$PY dwa dwa/axis_00055.txt 50 75 3.8 100 10 2260 1000 0.8 0.01 0.05

echo "--- 3. Cherenkov cone in an open dielectric"
"$EXE" "$SRC/tests/dielectric_cherenkov.deck" numerics.backend=$BACKEND > cher.log
$PY cherenkov cherenkov/map_00022.txt 2 140 10 60 0.02

echo "--- 4. laser in a uniform dielectric"
"$EXE" "$SRC/tests/dielectric_laser.deck" numerics.backend=$BACKEND > laser.log
sed '/^dielectric/d' "$SRC/tests/dielectric_laser.deck" > vlaser.deck
"$EXE" vlaser.deck numerics.backend=$BACKEND output.dir=vlaser > vlaser.log
$PY slip laser/axis_00000.txt laser/axis_00400.txt vlaser/axis_00000.txt vlaser/axis_00400.txt 1.1 200 0.025 0.5 0.8 1e-3
echo "dielectric test passed"
