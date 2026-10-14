#!/bin/sh
# Transport elements (tests/transport.deck: cold 1 fC bunch, gamma = 1000, B rho = 1.70451 T m,
# from z = 0 to 0.25 m). The ratio of the rms sizes at the end to the start against the
# linear optics (tests/transport_expect.awk integrates x'' = -K(z) x):
#   quadrupole 100 T/m, 5 cm long at 1 cm:   x 0.38374, y 1.64531 (+- 0.5 %)
#   plasma lens 100 T/m, the same:           x and y 0.38374
#   solenoid 5 T, 10 cm with 1 cm fringes:   r 0.96466 (Larmor frame, K = (B / 2 B rho)^2)
#   chicane at 10 cm (R56 = 100 um, delay 10 um) on a 4 um flat-top bunch with the chirp
#   delta = -0.01 (zeta - zeta_c) / 2 um: length halved (rms 1.1532 -> 0.5766 um +- 0.5 %),
#   centroid moved back by the delay; the self-field left at the old position < 2 % of its energy
# usage: sh tests/transport.sh <gempic_rz> <h5_extract> [backend]
set -e
abspath() { (cd "$(dirname "$1")" && echo "$(pwd)/$(basename "$1")"); }
EXE=$(abspath "$1"); H5X=$(abspath "$2"); BACKEND=${3:-openmp}
SRC=$(cd "$(dirname "$0")/.." && pwd)
D="$SRC/tests/transport.deck"
WORK=$(mktemp -d)
cd "$WORK"
check() {  # <label> <value> <expected> <relative tolerance>
  if awk -v v="$2" -v e="$3" -v t="$4" 'BEGIN { d = (v - e) / e; if (d < 0) d = -d; exit !(d <= t) }'; then
    echo "  $1 = $2 (expected $3)  ok"
  else
    echo "  $1 = $2, expected $3 +- $4 (relative)  FAILED"; exit 1
  fi
}
ratio() {  # <dir> <column>: last / first value of the bunch b
  awk -v c="$2" '$3 == "b" { if (!n++) a = $c; b = $c } END { printf "%.5f", b / a }' "$1/bunches.txt"
}
rratio() {  # r: sqrt(sx^2 + sy^2)
  awk '$3 == "b" { r = sqrt($8 * $8 + $9 * $9); if (!n++) a = r; b = r } END { printf "%.5f", b / a }' "$1/bunches.txt"
}
EXP="$SRC/tests/transport_expect.awk"
K=$(awk 'BEGIN { printf "%.10g", 100 / 1.7045081717848325 }')
QX=$(awk -v kind=quad_x -v k="$K" -v z0=0.01 -v z1=0.06 -v zend=0.25 -f "$EXP")
QY=$(awk -v kind=quad_y -v k="$K" -v z0=0.01 -v z1=0.06 -v zend=0.25 -f "$EXP")
KS=$(awk 'BEGIN { printf "%.10g", 5 / (2 * 1.7045081717848325) }')
SR=$(awk -v kind=sol -v k="$KS" -v z0=0.01 -v z1=0.11 -v ramp=0.01 -v zend=0.25 -f "$EXP")
O="numerics.backend=$BACKEND"
cp "$D" q.deck; echo 'quadrupole { z_start = 1 cm; length = 5 cm; gradient = 100 T/m }' >> q.deck
"$EXE" q.deck $O output.dir=q > q.log
check "quadrupole: x" "$(ratio q 8)" "$QX" 0.005
check "quadrupole: y" "$(ratio q 9)" "$QY" 0.005
cp "$D" l.deck; echo 'plasma_lens { z_start = 1 cm; length = 5 cm; gradient = 100 T/m }' >> l.deck
"$EXE" l.deck $O output.dir=l > l.log
check "plasma lens: x" "$(ratio l 8)" "$QX" 0.005
check "plasma lens: y" "$(ratio l 9)" "$QX" 0.005
cp "$D" s.deck; echo 'solenoid { z_start = 1 cm; length = 10 cm; field = 5 T; ramp = 1 cm }' >> s.deck
"$EXE" s.deck $O output.dir=s > s.log
check "solenoid: r" "$(rratio s)" "$SR" 0.003
# chicane
sed -e 's/^bunch .*/bunch { name = b; charge = 1 pC; sigma_x = 100 um; emittance_x = 0 um; distribution_z = function; profile_z = "1"; zeta_min = 28 um; zeta_max = 32 um/' \
    -e 's/^        sigma_z = 2 um.*/        gamma_function = "1000 * (1 - 0.01 * (zeta - 30 um) \/ (2 um))"; macroparticles = 20000 }/' "$D" > c.deck
echo 'chicane { z = 10 cm; r56 = 100 um; delay = 10 um; energy = 510.49 MeV }' >> c.deck
grep -q gamma_function c.deck
"$EXE" c.deck $O output.dir=c output.format=text,hdf5 output.every=1 grid.z_end=0.11m > c.log
grep -q "beam moved by -10 um" c.log
SZ=$(awk '$3 == "b" { s = $13 } END { print s }' c/bunches.txt)
check "chicane: rms length" "$SZ" 0.5766 0.005
ZC=$(awk '$3 == "b" { if (!n++) a = $12; b = $12 } END { printf "%.3f", a - b }' c/bunches.txt)
check "chicane: centroid moved back by" "$ZC" 10.055 0.002   # 10 um delay + slippage z / (2 gamma^2)
W0=$("$H5X" c/data_00020.h5 energy 26 36 | awk '!/^#/ { print $3 }')
W1=$("$H5X" c/data_00021.h5 energy 26 36 | awk '!/^#/ { print $3 }')
awk -v a="$W0" -v b="$W1" 'BEGIN { printf "  chicane: field left at the old position %.2g of its energy\n", b / a; exit !(b / a < 0.02) }'
cd - > /dev/null
rm -rf "$WORK"
echo "transport test passed ($BACKEND)"
