#!/bin/sh
# Protons and ions from openPMD files (bunch.file): charge and mass of the species from the file's records.
# tests/ion_beam.py writes an openPMD file (any code could have written it) with the species
#   c6     C6+ (charge 6 e, mass 12 u), cold, gamma = 100, 1 fC; charge and mass as constant records
#   c6w    C6+, warm (energy spread, emittance, unequal weights), 50 pC; charge, mass and momentum macroWeighted
#   p      protons, 400 GeV kinetic energy
#   mixed  two charge states in one species (refused);  neutral  charge 0 (refused)
# Checked against independent references:
#   (i)   c6w: Q, kinetic energy per ion, dE/E, rms sizes and normalized emittances (u = p / (m c), m = 12 u) in
#         bunches.txt = the moments computed by python from the file (1e-3, the printed digits)
#   (ii)  c6w written by the run as openPMD and as native HDF5 and read back: the same bunches.txt row (charge Z e,
#         weighting = ions, mass; native: bunch_masses, bunch_charges)
#   (iii) c6 through the quadrupole of tests/transport.sh with the gradient scaled by B rho = p / (6 e): the size
#         ratios of the linear optics (tests/transport_expect.awk, planes swapped by the positive charge), +- 0.5 %
#   (iv)  p: summary "protons", kinetic energy 400 GeV in bunches.txt
# usage: sh tests/ion_import.sh <gempic_rz> <compare_columns>     (exit 77: no python3 with numpy and h5py)
set -e
abspath() { (cd "$(dirname "$1")" && echo "$(pwd)/$(basename "$1")"); }
EXE=$(abspath "$1"); CMP=$(abspath "$2")
python3 -c "import numpy, h5py" 2>/dev/null || { echo "ion import test skipped: needs python3 with numpy and h5py"; exit 77; }
SRC=$(cd "$(dirname "$0")/.." && pwd)
WORK=$(mktemp -d)
cd "$WORK"
export OMP_PROC_BIND=false
O="numerics.backend=serial"
python3 "$SRC/tests/ion_beam.py" ions.h5 > ref.txt
check() {  # <label> <value> <expected> <relative tolerance>
  if awk -v v="$2" -v e="$3" -v t="$4" 'BEGIN { d = (v - e) / e; if (d < 0) d = -d; exit !(d <= t) }'; then
    echo "  $1 = $2 (expected $3)  ok"
  else
    echo "  $1 = $2, expected $3 +- $4 (relative)  FAILED"; exit 1
  fi
}
# the transport deck with a bunch from the file
sed -e '/^bunch {/,/macroparticles = 20000 }/d' "$SRC/tests/transport.deck" > base.deck
grep -q "^bunch" base.deck && { echo "bunch block of tests/transport.deck not removed"; exit 1; }
printf 'bunch { name = b; file = %s/ions.h5; zeta_center = 25 um }\n' "$WORK" >> base.deck
row() { awk -v s="$2" '$3 == "b" && $1 == s { print $4, $5, $6, $7, $8, $9, $10, $11, $12, $13 }' "$1"; }
# ---- (i) moments of the warm C6+ beam
"$EXE" base.deck $O bunch.file_species=c6w grid.z_end=5mm output.every=1 output.particles=true \
  output.format=text,hdf5,openpmd output.dir=w > w.log
grep -q "ions (charge +6 e, mass 12 u)" w.log || { grep "^bunch\[0\]" w.log; echo "summary: species not named"; exit 1; }
row w/bunches.txt 0 | awk '{ print $2, $3, $4, $5, $6, $7, $8 }' > got.txt
grep "^c6w" ref.txt | cut -d' ' -f2- > refw.txt
echo "  imported C6+ (bunches.txt): $(cat got.txt)"
echo "  from the file (python):     $(cat refw.txt)"
"$CMP" refw.txt got.txt 1e-3 1e-12
# ---- (ii) round trip through the openPMD and native output of the run
row w/bunches.txt 1 > ref1.txt
ZC=$(awk '{ print $9 }' ref1.txt)
for f in w/openpmd/openpmd_1.h5 w/data_00001.h5; do
  "$EXE" base.deck $O bunch.file=$f bunch.file_species=b bunch.zeta_center=${ZC}um grid.z_end=1k0^-1 \
    output.dir=rt > rt.log
  grep -q "ions (charge +6 e, mass 12 u)" rt.log || { grep "^bunch\[0\]" rt.log; echo "$f: species lost"; exit 1; }
  row rt/bunches.txt 0 > got1.txt
  "$CMP" ref1.txt got1.txt 1e-12 1e-12
done
echo "  written and read back (openPMD, native): same moments"
# ---- (iii) quadrupole on the cold C6+ beam, gradient scaled by B rho
BRE=1.7045081717848325
BRI=$(grep "^brho" ref.txt | cut -d' ' -f2)
G=$(awk -v a="$BRI" -v b="$BRE" 'BEGIN { printf "%.10g", 100 * a / b }')
K=$(awk -v b="$BRE" 'BEGIN { printf "%.10g", 100 / b }')
EXP="$SRC/tests/transport_expect.awk"
QX=$(awk -v kind=quad_x -v k="$K" -v z0=0.01 -v z1=0.06 -v zend=0.25 -f "$EXP")
QY=$(awk -v kind=quad_y -v k="$K" -v z0=0.01 -v z1=0.06 -v zend=0.25 -f "$EXP")
cp base.deck q.deck; echo "quadrupole { z_start = 1 cm; length = 5 cm; gradient = $G T/m }" >> q.deck
"$EXE" q.deck $O bunch.file_species=c6 output.dir=q > q.log
ratio() { awk -v c="$2" '$3 == "b" { if (!n++) a = $c; b = $c } END { printf "%.5f", b / a }' "$1/bunches.txt"; }
check "C6+ energy [MeV]" "$(awk '$3 == "b" { print $6; exit }' q/bunches.txt)" "$(grep "^ekin" ref.txt | cut -d' ' -f2)" 1e-5
check "C6+ quadrupole ($G T/m): x" "$(ratio q 8)" "$QY" 0.005
check "C6+ quadrupole ($G T/m): y" "$(ratio q 9)" "$QX" 0.005
# ---- (iv) protons
"$EXE" base.deck $O bunch.file_species=p grid.z_end=1k0^-1 output.dir=p > p.log
grep -q "'b': protons" p.log || { grep "^bunch\[0\]" p.log; echo "summary: not protons"; exit 1; }
check "proton energy [MeV]" "$(awk '$3 == "b" { print $6; exit }' p/bunches.txt)" 400000 1e-6
# ---- refused species
"$EXE" base.deck $O bunch.file_species=mixed grid.z_end=1k0^-1 output.dir=x > x.log 2>&1 && exit 1
grep -q "the charge differs between particles" x.log
"$EXE" base.deck $O bunch.file_species=neutral grid.z_end=1k0^-1 output.dir=x > x.log 2>&1 && exit 1
grep -q "is not a non-zero multiple of e" x.log
echo "  mixed charge states and neutral particles refused"
cd - > /dev/null
rm -rf "$WORK"
echo "ion import test passed"
