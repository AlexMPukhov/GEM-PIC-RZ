#!/bin/sh
# Metal wall (block 'wall'): a laser in a perfectly conducting circular waveguide (tests/wall.deck), against the
# modal theory of the guide (tests/wall.py: TE_1n and TM_1n modes with the dispersion relation of the scheme):
#  1. straight tube, a = 12 um: E_x to 0.8 %, E_z to 1.2 % rms of the maximum (no wall: 22 %; a wall at 12.5 um
#     against the 12 um reference: 23 %);
#  2. narrowing tube, a = 12 -> 9 um over 2 mm (adiabatic local modes): E_x 1.0 %, E_z 2.7 % (with E_z = 0 on the
#     sloped wall instead of E_z = -slope E_r the wall swallowed a third of the energy: peak 20 % low);
#  3. widening tube, 9 -> 12 um (w0 = 3.5 um): E_x 1.4 %, E_z 3.7 % (without the continuation of the field into the
#     opened annulus: E_z 5 times too large);
#  4. plasma in the tube: no plasma inside the metal, the density inside is the deck's.
# usage: sh tests/wall.sh <gempic_rz> [backend]     (exit 77: skipped, needs python3 with numpy)
set -e
abspath() { (cd "$(dirname "$1")" && echo "$(pwd)/$(basename "$1")"); }
EXE=$(abspath "$1")
BACKEND=${2:-serial}
SRC=$(cd "$(dirname "$0")/.." && pwd)
python3 -c "import numpy" 2>/dev/null || { echo "wall test skipped: needs python3 with numpy"; exit 77; }
WORK=$(mktemp -d)
cd "$WORK"
DECK="$SRC/tests/wall.deck"
export OMP_PROC_BIND=false

echo "--- 1. straight tube, a = 12 um"
"$EXE" "$DECK" numerics.backend=$BACKEND output.dir=straight > straight.log
python3 "$SRC/tests/wall.py" straight straight/map_00400.txt 0.8 2000 5 12 12 0.015 0.025

echo "--- 2. narrowing tube, 12 -> 9 um"
sed 's/^wall.*/wall { z_points = 0 2 mm; r_points = 12 9 um }/' "$DECK" > narrow.deck
"$EXE" narrow.deck numerics.backend=$BACKEND output.dir=narrow > narrow.log
python3 "$SRC/tests/wall.py" narrow narrow/map_00400.txt 0.8 2000 5 12 9 0.02 0.05

echo "--- 3. widening tube, 9 -> 12 um"
sed 's/^wall.*/wall { z_points = 0 2 mm; r_points = 9 12 um }/' "$DECK" > widen.deck
"$EXE" widen.deck numerics.backend=$BACKEND output.dir=widen laser.waist=3.5um > widen.log
python3 "$SRC/tests/wall.py" widen widen/map_00400.txt 0.8 2000 5 9 12 0.03 0.06

echo "--- 4. plasma in the tube"
cat > plasma.deck <<EOF
grid { lambda0 = 0.8 um; m_max = 0; radial_grid = uniform; dr = 0.2 um; r_max = 14 um
       dzeta = 0.1 um; window_length = 10 um; dz = 5 um; z_end = 10 um }
plasma { density = 1e18 cm^-3 }
wall { radius = 10.1 um }
output { dir = plasma; every = 1; format = text; map_stride = 10; particles = false }
EOF
"$EXE" plasma.deck numerics.backend=$BACKEND > plasma.log
awk '!/^#/ && NF > 4 { r = ($2 < 0) ? -$2 : $2
       if (r > 10.1) { if ($5 != 0) bad++ } else if (r < 9.5) { if ($5 < 0.98 || $5 > 1.02) bad++; n++ } }
     END { printf "density inside (r < 9.5 um, %d points) within 2 %% of 1, outside 0: %s\n", n, bad ? "FAIL" : "OK"
           exit bad ? 1 : 0 }' plasma/map_00002.txt
echo "wall test passed"
