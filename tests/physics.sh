#!/bin/sh
# Physics tests: bunch self-field in vacuum (stretched and free-hand radial grid, free-hand zeta grid),
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
"$EXE" "$SRC/tests/ion_collapse.deck" "$@" > ions.log
"$CHK" ioncollapse ion_collapse_out/axis_00001.txt
"$EXE" "$SRC/tests/ionization_H.deck" "$@" > ionH.log
"$CHK" ionH ionization_H_out/ionization.txt
"$EXE" "$SRC/tests/ionization_N.deck" "$@" > ionN.log
"$CHK" ionN ionization_N_out/ionization.txt
cd - > /dev/null
rm -rf "$WORK"
echo "physics tests passed"
