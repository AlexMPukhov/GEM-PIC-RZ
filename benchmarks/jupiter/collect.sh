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
# column 9 (newer runs): s per step after the pipeline has filled (rank 0)
awk 'NF >= 9 && $9 != "NA" { if (!h++) print "\nafter the pipeline filled (s/step):"; printf "  %-12s %6s GPUs %10s\n", $2, $1, $9 }' "$T"
W=$RUNS/weak.txt
if [ -f "$W" ]; then
  echo
  echo "weak scaling (1200 slices per GPU; efficiency = time per step after the fill on 1 GPU / on P GPUs):"
  awk '{ if ($6 != "NA" && $1 == 1) one = $6; n[NR] = $1; st[NR] = $4; s[NR] = $5; p[NR] = $6 }
    END { printf "  %6s %6s %10s %16s %8s\n", "GPUs", "steps", "wall s", "s/step (filled)", "effic."
          for (i = 1; i <= NR; ++i)
            printf "  %6s %6s %10s %16s %8s\n", n[i], st[i], s[i], p[i], (one && p[i] != "NA") ? sprintf("%.0f %%", 100 * one / p[i]) : "" }' "$W"
fi
