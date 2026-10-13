#!/bin/sh
# Physics tests: bunch self-field in vacuum (stretched and free-hand radial grid, free-hand zeta grid,
# gamma = 100: the co-moving field keeps the bunch energy),
# ion collapse, field ionization
# of H (one level) and N5+ -> N7+ (two levels); about 15 s in total.
# usage: sh tests/physics.sh <gempic_rz> <check_physics> [block.key=value ...]
#        (e.g. numerics.backend=openmp: the same tests for the parallel engine)
set -e
abspath() { (cd "$(dirname "$1")" && echo "$(pwd)/$(basename "$1")"); }
EXE=$(abspath "${1:-build/gempic_rz}")
CHK=$(abspath "${2:-build/check_physics}")
shift 2 2>/dev/null || true
SRC=$(cd "$(dirname "$0")/.." && pwd)
WORK=$(mktemp -d)
cd "$WORK"
"$EXE" "$SRC/tests/bunch_vacuum.deck" "$@" > bunch.log
"$CHK" bunchfield bunch_vacuum_out/map_00001.txt
"$EXE" "$SRC/tests/bunch_vacuum_table.deck" "$@" > bunch_table.log
"$CHK" bunchfield bunch_vacuum_table_out/map_00001.txt
"$EXE" "$SRC/tests/bunch_vacuum_zeta.deck" "$@" > bunch_zeta.log
"$CHK" bunchfield bunch_vacuum_zeta_out/map_00001.txt
"$EXE" "$SRC/tests/bunch_slow.deck" "$@" > bunch_slow.log
awk '!/^#/ && $3 == "b" { if (!n++) e0 = $6; e1 = $6 } END {
  d = (e1 - e0) / e0; ok = (d < 1e-4 && d > -1e-4)
  printf "%s bunch with gamma = 100 in vacuum: dE/E = %.2e over 4 cm (tolerance 1e-4)\n", ok ? "OK  " : "FAIL", d
  exit !ok }' bunch_slow_out/bunches.txt
"$EXE" "$SRC/tests/ion_collapse.deck" "$@" > ions.log
"$CHK" ioncollapse ion_collapse_out/axis_00001.txt
"$EXE" "$SRC/tests/ionization_H.deck" "$@" > ionH.log
"$CHK" ionH ionization_H_out/ionization.txt
"$EXE" "$SRC/tests/ionization_N.deck" "$@" > ionN.log
"$CHK" ionN ionization_N_out/ionization.txt
cd - > /dev/null
rm -rf "$WORK"
echo "physics tests passed"
