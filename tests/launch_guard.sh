#!/bin/sh
# Launch guard: when the launcher reports several tasks but MPI sees one rank (srun with an MPI library without
# Slurm/PMIx support, or a build without MPI), every task would run the whole deck and overwrite the files of the
# others. The run must stop at once with exit code 1 and a message (printed by the task with launcher index 0 only).
#  1. SLURM_STEP_NUM_TASKS = 3, SLURM_PROCID = 0: exit 1, message on stderr, no output directory.
#  2. SLURM_PROCID = 1 of the same step: exit 1, silent.
#  3. SLURM_STEP_NUM_TASKS = 1 (srun -n 1): the run goes ahead.
# usage: sh tests/launch_guard.sh <gempic_rz>
set -e
abspath() { (cd "$(dirname "$1")" && echo "$(pwd)/$(basename "$1")"); }
EXE=$(abspath "$1")
SRC=$(cd "$(dirname "$0")/.." && pwd)
WORK=$(mktemp -d)
cd "$WORK"
unset OMPI_COMM_WORLD_SIZE OMPI_COMM_WORLD_RANK PMI_SIZE PMI_RANK
D="$SRC/tests/regression.deck"
G="grid.z_end=40k0^-1 output.every=0"

rc=0
SLURM_STEP_NUM_TASKS=3 SLURM_PROCID=0 "$EXE" "$D" $G output.dir=g0 > g0.out 2> g0.err || rc=$?
[ "$rc" -eq 1 ] || { echo "  3 tasks, index 0: exit code $rc, expected 1: FAILED"; cat g0.out g0.err; exit 1; }
grep -q "the launcher started 3 tasks, but MPI sees only one rank" g0.err ||
  { echo "  3 tasks, index 0: no message: FAILED"; cat g0.err; exit 1; }
[ ! -e g0 ] || { echo "  3 tasks: output directory written: FAILED"; exit 1; }
echo "  3 tasks, index 0: stopped with exit code 1 and the message"

rc=0
SLURM_STEP_NUM_TASKS=3 SLURM_PROCID=1 "$EXE" "$D" $G output.dir=g1 > g1.out 2> g1.err || rc=$?
[ "$rc" -eq 1 ] && [ ! -s g1.err ] || { echo "  3 tasks, index 1: exit $rc, stderr '$(cat g1.err)': FAILED"; exit 1; }
echo "  3 tasks, index 1: stopped silently"

SLURM_STEP_NUM_TASKS=1 SLURM_PROCID=0 "$EXE" "$D" $G output.dir=g2 > g2.out 2> g2.err ||
  { echo "  1 task: the run failed: FAILED"; cat g2.out g2.err; exit 1; }
grep -q "^done" g2.out || { echo "  1 task: no 'done' line: FAILED"; cat g2.out; exit 1; }
echo "  1 task: runs"
rm -rf "$WORK"
echo "launch guard: passed"
