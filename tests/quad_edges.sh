#!/bin/sh
# Element edges between the sub-steps of the beam push (tests/transport.deck: cold 1 fC bunch, gamma = 1000,
# B rho = 1.70451 T m). The push averages the element field over each sub-step, so an edge anywhere inside a
# sub-step gives the exact integrated kick. Edges here lie off the step grid, sub-steps are 5 mm (= dz):
#   FODO: 6 quadrupoles of 4.37 cm, +-100 T/m, back to back from z = 1.13 cm, drift to 0.30 m;
#   hard-edge solenoid (ramp = 0): 5 T, 10.37 cm from z = 1.13 cm, drift to 0.25 m (Larmor frame,
#   K = (B / 2 B rho)^2 inside: the edge fields B_r are part of it).
# Reference: tests/transport_expect.awk (RK4 of x'' = -K(z) x, the ratio of the rms sizes of a cold beam).
# With the field sampled at points (before) the FODO was off by several per cent and the hard-edge solenoid
# did not focus at all; now both agree to < 0.5 %.
# usage: sh tests/quad_edges.sh <gempic_rz> [backend]
set -e
abspath() { (cd "$(dirname "$1")" && echo "$(pwd)/$(basename "$1")"); }
EXE=$(abspath "$1"); BACKEND=${2:-openmp}
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
BRHO=1.7045081717848325
O="numerics.backend=$BACKEND numerics.beam_dz_max=5mm"
# FODO with edges off the 5 mm grid
K=$(awk -v b=$BRHO 'BEGIN { printf "%.10g", 100 / b }')
cp "$D" f.deck
ELEMS=""
z=0.0113; s=1
for i in 1 2 3 4 5 6; do
  z1=$(awk -v z=$z 'BEGIN { printf "%.4f", z + 0.0437 }')
  echo "quadrupole { z_start = $z m; length = 0.0437 m; gradient = $((s * 100)) T/m }" >> f.deck
  ELEMS="$ELEMS$z $z1 $(awk -v k=$K -v s=$s 'BEGIN { printf "%.10g", s * k }');"
  z=$z1; s=$((-s))
done
QX=$(awk -v kind=list -v elems="$ELEMS" -v zend=0.30 -f "$EXP")
QY=$(awk -v kind=list -v sign=-1 -v elems="$ELEMS" -v zend=0.30 -f "$EXP")
"$EXE" f.deck $O grid.z_end=0.30m output.dir=f > f.log
check "FODO, edges off the grid: x" "$(ratio f 8)" "$QX" 0.005
check "FODO, edges off the grid: y" "$(ratio f 9)" "$QY" 0.005
# hard-edge solenoid
KS=$(awk -v b=$BRHO 'BEGIN { printf "%.10g", 5 / (2 * b) }')
SR=$(awk -v kind=sol -v k="$KS" -v z0=0.0113 -v z1=0.115 -v ramp=0 -v zend=0.25 -f "$EXP")
cp "$D" s.deck; echo 'solenoid { z_start = 1.13 cm; length = 10.37 cm; field = 5 T }' >> s.deck
"$EXE" s.deck $O output.dir=s > s.log
check "hard-edge solenoid: r" "$(rratio s)" "$SR" 0.005
cd - > /dev/null
rm -rf "$WORK"
echo "quad_edges test passed ($BACKEND)"
