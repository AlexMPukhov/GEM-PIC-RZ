#!/bin/sh
# Ion blocks with their own density profile (ion.profile = own):
#   1. tests/ionization_H.deck with the hydrogen given by its own profile (same as the
#      plasma's): identical results
#   2. tests/regression.deck with the plasma electrons and H+ ions given as a pre-ionized
#      own-profile species (free_electrons = false; the code loads Z electrons per atom):
#      identical results
#   3. tests/ionization_N.deck with the N5+ dopant only at 3 < z < 6 k0^-1: ionization
#      only in the steps 4 ... 6, as strong as in the uniform case (+- 5 %)
# usage: sh tests/dopant.sh <gempic_rz> <compare_columns> [backend]
set -e
abspath() { (cd "$(dirname "$1")" && echo "$(pwd)/$(basename "$1")"); }
EXE=$(abspath "$1"); CMP=$(abspath "$2"); BACKEND=${3:-serial}
SRC=$(cd "$(dirname "$0")/.." && pwd)
WORK=$(mktemp -d)
cd "$WORK"
O="numerics.backend=$BACKEND"
# serial / openmp: identical; cuda (atomic deposition in a nondeterministic order): equal to rounding
if [ "$BACKEND" = cuda ]; then T0="1e-6 1e-9"; T1="1e-6 1e-9"; else T0="0 0"; T1="1e-9 1e-12"; fi
# 1
H="$SRC/tests/ionization_H.deck"
sed 's/^ion    { element = H;  ionization = adk;  density = 1e-8 n_c;  mobile = false }/ion    { element = H;  ionization = adk;  profile = own;  density = 1e-8 n_c;  r_max = 50 k0^-1;  mobile = false }/' \
  "$H" > h.deck
grep -q "profile = own" h.deck
"$EXE" "$H" $O output.dir=h0 > h0.log
"$EXE" h.deck $O output.dir=h1 > h1.log
for f in ionization.txt axis_00001.txt; do "$CMP" h0/$f h1/$f $T0; done
# 2
R="$SRC/tests/regression.deck"
sed -e 's/^  density       = 0.002 n_c/  density       = 0.002 n_c\n  free_electrons = false/' \
    -e 's/^  element = H/  element = H\n  profile = own\n  density = 0.002 n_c\n  z_points = 0, 400 k0^-1\n  z_values = 0, 1\n  r_max = 150 k0^-1/' \
    "$R" > pre.deck
grep -q "profile = own" pre.deck
"$EXE" "$R" $O grid.z_end=80k0^-1 output.every=2 output.dir=p0 > p0.log
"$EXE" pre.deck $O grid.z_end=80k0^-1 output.every=2 output.dir=p1 > p1.log
for f in axis_00002.txt map_00002.txt bunches.txt; do "$CMP" p0/$f p1/$f $T1; done
# 3
N="$SRC/tests/ionization_N.deck"
sed 's/^ion    { element = N;  charge_state = 5;  ionization = adk }/ion    { element = N;  charge_state = 5;  ionization = adk;  profile = own;  density = 0.2e-8 n_c\n         z_points = 3, 3.01, 5.99, 6 k0^-1;  z_values = 0, 1, 1, 0;  r_max = 50 k0^-1 }/' \
  "$N" > n.deck
grep -q "profile = own" n.deck
"$EXE" "$N" $O grid.z_end=10k0^-1 output.dir=n0 > n0.log
"$EXE" n.deck $O grid.z_end=10k0^-1 output.dir=n1 > n1.log
paste n0/ionization.txt n1/ionization.txt | awk '!/^#/ {
    inside = ($1 >= 4 && $1 <= 6)
    if (!inside && $(NF/2 + 5) != 0) { print "  ionization outside the dopant at step " $1; bad = 1 }
    if (inside) { a += $5; b += $(NF/2 + 5) } }
  END { printf "  freed charge in the dopant steps: %.4e pC (uniform dopant: %.4e pC)\n", b, a
        exit !(!bad && b > 0.95 * a && b < 1.05 * a) }'
cd - > /dev/null
rm -rf "$WORK"
echo "dopant test passed ($BACKEND)"
