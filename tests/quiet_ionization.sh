#!/bin/sh
# Quiet ionization (ion block key split = K): stratified thresholds of the K
# sibling ions of a loading slot (src/rng.hpp quietDepth).
#   1. Hydrogen (tests/ionization_H.deck, one level), 6 seeds each, the same
#      number of ion macroparticles: Monte Carlo (split = 1, ppc 2 x 32) and
#      quiet (split = 16, ppc 2 x 2). Both means within 5 % of the ADK reference
#      (check_physics ionH); the seed-to-seed scatter of the ionized fraction
#      of the quiet runs below 1/3 of the Monte Carlo scatter.
#   2. N5+ -> N7+ (tests/ionization_N.deck, two sequential levels) with
#      split = 8 and ppc 4 x 4: charge-state fractions vs the rate equations
#      (check_physics ionN tolerances): the stratification keeps the mean rates.
# usage: sh tests/quiet_ionization.sh <gempic_rz> <check_physics> [block.key=value ...]
set -e
abspath() { (cd "$(dirname "$1")" && echo "$(pwd)/$(basename "$1")"); }
EXE=$(abspath "${1:-build/gempic_rz}")
CHK=$(abspath "${2:-build/check_physics}")
shift 2 2>/dev/null || true
SRC=$(cd "$(dirname "$0")/.." && pwd)
WORK=$(mktemp -d)
cd "$WORK"
H="$SRC/tests/ionization_H.deck"
for seed in 1 2 3 4 5 6; do
  "$EXE" "$H" "$@" plasma.ppc_r=2 plasma.ppc_phi=32 numerics.seed=$seed output.dir=mc$seed > mc$seed.log
  "$CHK" ionH mc$seed/ionization.txt >> mc.txt || true
  "$EXE" "$H" "$@" plasma.ppc_r=2 plasma.ppc_phi=2 ion.split=16 numerics.seed=$seed output.dir=q$seed > q$seed.log
  "$CHK" ionH q$seed/ionization.txt >> q.txt || true
done
grep -q "quiet ionization: 16 ion macroparticles" q1.log || { echo "FAIL split = 16 not in the summary"; exit 1; }
# lines: "<OK|FAIL> H ionized fraction (disk r < 50): <f>, ADK reference <ref> (tolerance 5%)"
stats() { awk '{ f = $7; sub(",", "", f); r = $10; n++; s += f; s2 += f * f }
  END { m = s / n; printf "%.6g %.6g %.6g", m, sqrt((s2 - n * m * m) / (n - 1)), r }' "$1"; }
set -- $(stats mc.txt) $(stats q.txt)
awk -v mm=$1 -v ms=$2 -v ref=$3 -v qm=$4 -v qs=$5 'BEGIN {
  okm = (mm / ref - 1 < 0.05 && mm / ref - 1 > -0.05); okq = (qm / ref - 1 < 0.05 && qm / ref - 1 > -0.05)
  oks = (qs < ms / 3)
  printf "%s Monte Carlo (split 1, 2 x 32 per cell): ionized fraction %.5f +- %.2e (6 seeds), ADK reference %.5f\n",
    okm ? "OK  " : "FAIL", mm, ms, ref
  printf "%s quiet (split 16, 2 x 2 per cell):       ionized fraction %.5f +- %.2e (6 seeds)\n", okq ? "OK  " : "FAIL", qm, qs
  printf "%s scatter quiet / Monte Carlo = %.3f (must be < 1/3)\n", oks ? "OK  " : "FAIL", qs / ms
  exit !(okm && okq && oks) }'
"$EXE" "$SRC/tests/ionization_N.deck" "$@" plasma.ppc_r=4 plasma.ppc_phi=4 ion.split=8 output.dir=qn > qn.log
"$CHK" ionN qn/ionization.txt
cd - > /dev/null
rm -rf "$WORK"
echo "quiet ionization tests passed"
