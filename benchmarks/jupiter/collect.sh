#!/bin/sh
# Table of the scaling series from $RUNS/scaling.txt: time per step, speed-up against one GPU (the pipelined
# step on one rank), parallel efficiency, and the ideal efficiency of the pipeline N / (N + P - 1) for N steps.
. "$(dirname "$0")/env.sh"
T=$RUNS/scaling.txt
[ -f "$T" ] || { echo "no $T yet"; exit 1; }
awk '
  { ngpu = $1; label = $2; ranks = $3; steps = $4; sec = $5
    if (steps == "NA" || sec == "NA" || steps == 0) { bad[NR] = $0; next }
    base = label; sub(/_.*/, "", base)
    key[NR] = label; n[NR] = ngpu; st[NR] = steps; s[NR] = sec; r[NR] = ranks; b[NR] = base
    if (ngpu == 1 && label == base) one[base] = sec / steps }
  END {
    printf "%-12s %6s %8s %6s %10s %9s %9s %9s\n", "run", "GPUs", "ranks", "steps", "s/step", "speed-up", "effic.", "pipeline"
    for (i = 1; i <= NR; ++i) {
      if (!(i in key)) continue
      ps = s[i] / st[i]; P = (n[i] == "cpu") ? 1 : n[i]
      sp = (b[i] in one) ? one[b[i]] / ps : 0
      eff = (n[i] == "cpu" || !(b[i] in one)) ? "" : sprintf("%.0f %%", 100 * sp / P)
      pipe = (n[i] == "cpu") ? "" : sprintf("%.0f %%", 100 * st[i] / (st[i] + P - 1))
      printf "%-12s %6s %8s %6d %10.3f %9s %9s %9s\n", key[i], n[i], r[i], st[i], ps, sp ? sprintf("%.2f", sp) : "", eff, pipe
    }
    for (i in bad) print "incomplete: " bad[i]
  }' "$T"
echo "(speed-up and efficiency against the pipelined step on 1 GPU; *_engine: the single-GPU engine without the"
echo " pipeline; *_cpu: one node without GPUs; 'pipeline': the ideal efficiency N / (N + P - 1) of the pipeline)"
