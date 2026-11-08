#!/bin/sh
# Output volume and bounded memory (a DWA tube, tests/dielectric_dwa.deck with a shorter window and 400 macroparticles
# written per output step).
#  1. The volume report of the first output step ("output volume: X GB per output step in '<dir>', N more output
#     steps") against the files actually written: X = the bytes of the files of the last output step (axis_, map_,
#     beam_NNNNN.txt) within 1 %, N = the number of output steps after the first.
#  2. Bounded memory: the peak resident memory ("peak memory: X GB", getrusage) of a run of 8 x the steps (160) is
#     that of the short run (20) within 5 % + 20 MB (scatter of the peak under load ~10 MB; a list or buffer that grows
#     with the step count by more than ~0.15 MB per step would show here).
#  3. Linux, /dev/shm on tmpfs: an output directory there is reported (the files take RAM).
# usage: sh tests/output_volume.sh <gempic_rz> [backend]
set -e
abspath() { (cd "$(dirname "$1")" && echo "$(pwd)/$(basename "$1")"); }
EXE=$(abspath "$1"); BACKEND=${2:-serial}
SRC=$(cd "$(dirname "$0")/.." && pwd)
WORK=$(mktemp -d)
cd "$WORK"
D="$SRC/tests/dielectric_dwa.deck"
G="numerics.backend=$BACKEND grid.window_length=400um bunch.zeta_center=360um bunch.macroparticles=400"

# 1. volume report against the files
"$EXE" "$D" $G grid.z_end=0.4mm output.every=5 output.particles=true output.dir=v > v.log
X=$(sed -n "s/^output volume: \([0-9.e+-]*\) GB per output step.*/\1/p" v.log)
N=$(sed -n "s/^output volume: .*', \([0-9]*\) more output steps.*/\1/p" v.log)
[ -n "$X" ] && [ -n "$N" ] || { echo "no output volume report in the console"; cat v.log; exit 1; }
LAST=$(ls v/axis_*.txt | tail -1 | sed 's/.*axis_\([0-9]*\)\.txt/\1/')
BYTES=$(cat v/axis_$LAST.txt v/map_$LAST.txt v/beam_$LAST.txt | wc -c)
FILES=$(ls v/axis_*.txt | wc -l)
echo "  reported $X GB per output step, measured $BYTES bytes (step $LAST); $N more output steps, written $FILES"
awk -v x="$X" -v b="$BYTES" 'BEGIN { d = x * 1e9 / b - 1; if (d < 0) d = -d; exit !(d < 0.01) }' ||
  { echo "  volume per output step: FAILED"; exit 1; }
[ "$N" -eq $((FILES - 1)) ] || { echo "  number of output steps: FAILED"; exit 1; }

# 2. peak memory of 20 and 160 steps
"$EXE" "$D" $G grid.z_end=0.4mm output.every=0 output.dir=m1 > m1.log
"$EXE" "$D" $G grid.z_end=3.2mm output.every=0 output.dir=m4 > m4.log
P1=$(sed -n 's/^peak memory: \([0-9.e+-]*\) GB.*/\1/p' m1.log)
P4=$(sed -n 's/^peak memory: \([0-9.e+-]*\) GB.*/\1/p' m4.log)
if [ -z "$P1" ] || [ -z "$P4" ]; then
  echo "  no peak memory reported (not a unix system): memory check skipped"
else
  S1=$(grep -c '^step' m1.log); S4=$(grep -c '^step' m4.log)
  echo "  peak memory: $P1 GB after $((S1 - 1)) steps, $P4 GB after $((S4 - 1)) steps"
  awk -v a="$P1" -v b="$P4" 'BEGIN { exit !(b <= a * 1.05 + 0.02) }' || { echo "  memory grows with the steps: FAILED"; exit 1; }
fi

# 3. output on tmpfs
if [ "$(uname)" = Linux ] && [ "$(stat -f -c %T /dev/shm 2>/dev/null)" = tmpfs ]; then
  T=/dev/shm/gempic_volume_$$
  "$EXE" "$D" $G grid.z_end=0.04mm output.every=0 output.dir=$T > t.log
  rm -rf "$T"
  grep -q "^WARNING: output.dir '$T' is on tmpfs" t.log || { echo "  no tmpfs warning: FAILED"; cat t.log; exit 1; }
  echo "  output directory on tmpfs reported"
fi
echo "output_volume: ok"
cd - > /dev/null
rm -rf "$WORK"
