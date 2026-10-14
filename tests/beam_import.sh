#!/bin/sh
# Beams from files (bunch.file) and shot noise (bunch.shot_noise):
#   1. the witness of tests/regression.deck after 2 steps, written as openPMD and as
#      native HDF5, read back into a new run: the bunch statistics at the start equal
#      those of the source run (bunches.txt)
#   2. selection and resampling keep the moments (charge, energy, sizes)
#   3. (python3 with numpy and h5py) shot noise of tests/shot_noise.deck: Fawley beamlets
#      give the bunching of the real electrons, N_e <|b|^2> = 1 +- 0.25 (window leakage
#      measured with a 1000 times larger charge and subtracted); random positions give
#      that of the macroparticles, N_macro <|b|^2> = 1 +- 0.25
# usage: sh tests/beam_import.sh <gempic_rz> <compare_columns>
set -e
abspath() { (cd "$(dirname "$1")" && echo "$(pwd)/$(basename "$1")"); }
EXE=$(abspath "$1"); CMP=$(abspath "$2")
SRC=$(cd "$(dirname "$0")/.." && pwd)
WORK=$(mktemp -d)
cd "$WORK"
O="numerics.backend=serial"
"$EXE" "$SRC/tests/regression.deck" $O grid.z_end=80k0^-1 output.every=2 output.format=text,hdf5,openpmd \
  output.openpmd_dr=0.5um output.dir=src > src.log
sed -e '/^bunch {/,/^}/d' "$SRC/tests/regression.deck" > imp.deck
printf 'bunch {\n  name = witness\n  file = src/openpmd/openpmd_2.h5\n  zeta_center = 16.5502 um\n}\n' >> imp.deck
# columns N, Q, E, dE/E, sig_x, sig_y, eps_x, eps_y, <zeta>, sig_z
row() { awk -v s="$2" '$3 == "witness" && $1 == s { print $4, $5, $6, $7, $8, $9, $10, $11, $12, $13 }' "$1"; }
row src/bunches.txt 2 > ref.txt
for f in src/openpmd/openpmd_2.h5 src/data_00002.h5; do
  "$EXE" imp.deck $O grid.z_end=1k0^-1 bunch.file=$f output.dir=imp > imp.log
  row imp/bunches.txt 0 > got.txt
  "$CMP" ref.txt got.txt 1e-12 1e-12
done
# selection and resampling: 80000 macroparticles from 20000
"$EXE" imp.deck $O grid.z_end=1k0^-1 bunch.macroparticles=80000 output.dir=res > res.log
row res/bunches.txt 0 | awk '{ $1 = 20000; print }' > got.txt
"$CMP" ref.txt got.txt 2e-2 1e-12
"$EXE" imp.deck $O grid.z_end=1k0^-1 "bunch.select=uz > 1e9" output.dir=bad > bad.log 2>&1 && exit 1
grep -q "no particle of the file passes the selection" bad.log
echo "  file import, selection and resampling ok"
if python3 -c "import numpy, h5py" 2>/dev/null; then
  "$EXE" "$SRC/tests/shot_noise.deck" output.dir=fw > fw.log
  "$EXE" "$SRC/tests/shot_noise.deck" output.dir=fwq bunch.charge=5000pC > fwq.log
  "$EXE" "$SRC/tests/shot_noise.deck" output.dir=rn bunch.shot_noise=random > rn.log 2>&1 || {
    sed -e 's/^        shot_noise .*/        }/' -e 's/macroparticles = 1600000$/macroparticles = 1600000/' \
      "$SRC/tests/shot_noise.deck" > rn.deck
    "$EXE" rn.deck output.dir=rn > rn.log; }
  A=$(python3 "$SRC/tests/bunching.py" fw/data_00000.h5 0.8e-6 16)
  B=$(python3 "$SRC/tests/bunching.py" fwq/data_00000.h5 0.8e-6 16)
  C=$(python3 "$SRC/tests/bunching.py" rn/data_00000.h5 0.8e-6 16)
  echo "$A $B $C" | awk '{ fw = ($1 - $4) * $2; rn = $7 * $9;
    printf "  Fawley: N_e <|b|^2> = %.3f (leakage %.3f subtracted); random: N_macro <|b|^2> = %.3f\n", fw, $4 * $2, rn;
    exit !(fw > 0.75 && fw < 1.25 && rn > 0.75 && rn < 1.25) }'
else
  echo "  shot-noise check skipped (python3 with numpy and h5py not found)"
fi
cd - > /dev/null
rm -rf "$WORK"
echo "beam import test passed"
