#!/bin/sh
# GEM-PIC-RZ on a GPU: correctness and speed. Run from benchmarks/gpu after building build-gpu (README.md).
#   sh gpu_check.sh [path to build-gpu]       (default ../../build-gpu)
# 1. ctest of the GPU build (regression_cuda, physics_cuda, restart_cuda, ...)
# 2. the LWFA benchmark (../lwfa_fbpic/gem.deck, 100 um) on the GPU, compared with FBPIC by compare.py
#    as the CPU run was (CPU: on-axis E_z of the first period to 1.6 % rms)
# 3. speed: lwfa_heavy.deck (200 steps, ~4e4 plasma particles per slice) with GEMPIC_PROFILE=1
# 4. speed, beam-dominated: the seeded planar FEL (../harmonics/fel.deck: 640 000 beam particles, 3 modes,
#    18 372 slices; 5.6 GB), first 96 steps only
set -e
B=$(cd "${1:-../../build-gpu}" && pwd)
X=$B/gempic_rz
nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv,noheader || true
echo "== 1. ctest"
(cd "$B" && ctest 2>&1 | tail -3)
echo "== 2. LWFA benchmark vs FBPIC"
(cd ../lwfa_fbpic && /usr/bin/time -f "wall %e s" "$X" gem.deck numerics.backend=cuda output.dir=gem_out_gpu > gem_out_gpu.log 2>&1 \
  && grep "^done" gem_out_gpu.log && python3 compare.py gem_out_gpu | tail -6 && python3 compare.py gem_out | tail -6)
echo "== 3. speed, plasma-dominated (lwfa_heavy.deck, 200 steps)"
GEMPIC_PROFILE=1 "$X" lwfa_heavy.deck numerics.backend=cuda > heavy.log 2>&1; grep "^done" heavy.log; tail -25 heavy.log | grep -iv "^step\|^ *b " || true
echo "== 4. speed, beam-dominated (planar FEL, 96 steps)"
(cd ../harmonics && GEMPIC_PROFILE=1 "$X" fel.deck numerics.backend=cuda grid.z_end=30mm output.every=0 output.dir=fel_gpu > fel_gpu.log 2>&1; \
  grep "^done" fel_gpu.log; awk '!/^#/ && $1 == 96 && $3 == 1 && $4 == 1 {print "GPU: P1 in slice 14 at step 96:", $19, "W"}' fel_gpu/fel_harmonics.txt;
  echo "CPU: P1 in slice 14 at step 96: 9.979e+06 W")
