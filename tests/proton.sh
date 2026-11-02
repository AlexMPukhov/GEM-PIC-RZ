#!/bin/sh
# Proton bunches (bunch.species = proton): mass m_p, charge +e.
# tests/transport.deck with a cold 1 fC proton bunch of 400 GeV kinetic energy (gamma = 427.3156,
# B rho = p / e = 1337.3825 T m). The gradients are scaled by B rho_p / B rho_e (electrons of gamma = 1000,
# B rho = 1.70451 T m), so K = G / B rho is the one of tests/transport.sh and the size ratios are those of the
# linear optics (tests/transport_expect.awk), but with the roles of the planes swapped by the positive charge:
#   quadrupole +78460 T/m (focuses electrons in x):  protons x 1.64531, y 0.38374 (+- 0.5 %)
#   plasma lens -78460 T/m (defocuses electrons):     protons x and y 0.38374
# and the energy in bunches.txt is the kinetic energy (gamma - 1) m_p c^2 = 400 GeV.
# usage: sh tests/proton.sh <gempic_rz> [backend]
set -e
abspath() { (cd "$(dirname "$1")" && echo "$(pwd)/$(basename "$1")"); }
EXE=$(abspath "$1"); BACKEND=${2:-openmp}
SRC=$(cd "$(dirname "$0")/.." && pwd)
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
EXP="$SRC/tests/transport_expect.awk"
BRE=1.7045081717848325
# B rho of the protons from m_p c^2 = 938.27208816 MeV and 400 GeV kinetic energy
BRP=$(awk 'BEGIN { m = 938.27208816; E = 400000 + m; printf "%.10g", sqrt(E * E - m * m) * 1e6 / 299792458 }')
G=$(awk -v a="$BRP" -v b="$BRE" 'BEGIN { printf "%.8g", 100 * a / b }')
K=$(awk -v b="$BRE" 'BEGIN { printf "%.10g", 100 / b }')
QX=$(awk -v kind=quad_x -v k="$K" -v z0=0.01 -v z1=0.06 -v zend=0.25 -f "$EXP")
QY=$(awk -v kind=quad_y -v k="$K" -v z0=0.01 -v z1=0.06 -v zend=0.25 -f "$EXP")
sed 's/^bunch { name = b; charge = 1 fC; gamma = 1000;/bunch { name = b; species = proton; charge = 1 fC; energy = 400 GeV;/' \
  "$SRC/tests/transport.deck" > p.deck
grep -q "species = proton" p.deck
O="numerics.backend=$BACKEND"
cp p.deck q.deck; echo "quadrupole { z_start = 1 cm; length = 5 cm; gradient = $G T/m }" >> q.deck
"$EXE" q.deck $O output.dir=q > q.log
grep -q "protons, Q = 0.001 pC, E = 4e+05 MeV (gamma 427.32" q.log
check "energy [MeV]" "$(awk '$3 == "b" { e = $6 } END { print e }' q/bunches.txt)" 400000 1e-6
check "quadrupole: x" "$(ratio q 8)" "$QY" 0.005
check "quadrupole: y" "$(ratio q 9)" "$QX" 0.005
cp p.deck l.deck; echo "plasma_lens { z_start = 1 cm; length = 5 cm; gradient = -$G T/m }" >> l.deck
"$EXE" l.deck $O output.dir=l > l.log
check "plasma lens: x" "$(ratio l 8)" "$QX" 0.005
check "plasma lens: y" "$(ratio l 9)" "$QX" 0.005
cd - > /dev/null
rm -rf "$WORK"
echo "proton test passed ($BACKEND)"
