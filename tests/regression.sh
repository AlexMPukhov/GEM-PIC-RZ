#!/bin/sh
# Regression test: run tests/regression.deck and compare with the reference
# output in tests/reference/ (serial backend, -O3 without -march=native).
# The reference was regenerated when the beam current filter
# (numerics.beam_filter = 1) and the zero-gradient closure at the rear end of
# the window were introduced, and again with the self-fields of bunches with
# v < c (gamma = 1000 here: changes of up to 2e-3 near the rear of the window).
# usage: sh tests/regression.sh <gempic_rz> <compare_columns> [block.key=value ...]
#        (e.g. numerics.backend=cuda: the parallel engine against the same reference)
set -e
abspath() { (cd "$(dirname "$1")" && echo "$(pwd)/$(basename "$1")"); }
EXE=$(abspath "${1:-build/gempic_rz}")
CMP=$(abspath "${2:-build/compare_columns}")
shift 2 2>/dev/null || true
SRC=$(cd "$(dirname "$0")/.." && pwd)
WORK=$(mktemp -d)
cd "$WORK"
"$EXE" "$SRC/tests/regression.deck" "$@" > run.log
# tolerances allow for FMA / compiler differences (bitwise identical with -O3, no -march=native)
# and for the order of the atomic sums of the GPU engine
"$CMP" "$SRC/tests/reference/axis_00005.txt" regression_out/axis_00005.txt 1e-6 1e-9
"$CMP" "$SRC/tests/reference/bunches.txt"    regression_out/bunches.txt    1e-6 1e-9
cd - > /dev/null
rm -rf "$WORK"
echo "regression test passed"
