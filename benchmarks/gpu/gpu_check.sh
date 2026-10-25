#!/bin/sh
# GEM-PIC-RZ on a GPU: correctness and speed. Run from benchmarks/gpu after building build-gpu (README.md).
#   sh gpu_check.sh [path to build-gpu]       (default ../../build-gpu)
# 1. ctest of the GPU build (regression_cuda, physics_cuda, restart_cuda, ...)
# 2. the LWFA benchmark (../lwfa_fbpic/gem.deck, 100 um) on the GPU, compared with FBPIC by compare.py
#    as the CPU run was (CPU: on-axis E_z of the first period to 1.6 % rms)
# 3. speed: lwfa_heavy.deck (200 steps, ~4e4 plasma particles per slice) with GEMPIC_PROFILE=1
# 4. speed, beam-dominated: the seeded planar FEL (../harmonics/fel.deck: 640 000 beam particles, 3 modes,
#    18 372 slices; 5.6 GB), first 30 mm only (48 steps at dz = lambda_u / 16)
# Portable: no /usr/bin/time (missing on some clusters), python3 only for the comparison plots.
B=$(cd "${1:-../../build-gpu}" && pwd) || exit 1
X=$B/gempic_rz
nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv,noheader || true
echo "== 1. ctest"
(cd "$B" && ctest 2>&1 | tail -15)   # summary with the names of failed tests
echo "== 2. LWFA benchmark vs FBPIC"
(cd ../lwfa_fbpic
 t0=$(date +%s)
 "$X" gem.deck numerics.backend=cuda output.dir=gem_out_gpu > gem_out_gpu.log 2>&1
 echo "exit code $?, wall $(( $(date +%s) - t0 )) s"; grep "^done" gem_out_gpu.log
 if python3 -c "import numpy, matplotlib" 2>/dev/null; then
   python3 compare.py gem_out_gpu | tail -6; echo "-- CPU run for comparison:"; python3 compare.py gem_out | tail -6
 else
   echo "python3 with numpy + matplotlib not found: copy benchmarks/lwfa_fbpic/gem_out_gpu back and run compare.py there"
 fi)
echo "== 3. speed, plasma-dominated (lwfa_heavy.deck, 200 steps)"
GEMPIC_PROFILE=1 "$X" lwfa_heavy.deck numerics.backend=cuda > heavy.log 2>&1; grep "^done" heavy.log; tail -25 heavy.log | grep -iv "^step\|^ *b " || true
echo "== 4. speed, beam-dominated (planar FEL, 30 mm = 48 steps)"
(cd ../harmonics && GEMPIC_PROFILE=1 "$X" fel.deck numerics.backend=cuda grid.z_end=30mm output.every=0 output.dir=fel_gpu > fel_gpu.log 2>&1; \
  grep "^done" fel_gpu.log; awk '!/^#/ && $1 == 48 && $3 == 1 && $4 == 1 {print "GPU: P1 in slice 14 at step 48:", $19, "W"}' fel_gpu/fel_harmonics.txt;
  echo "CPU: P1 in slice 14 at step 48: 9.96529e+06 W")
