#!/bin/sh
# resolution scan of emission.deck
export OMP_PROC_BIND=false OMP_NUM_THREADS=2
X=../../build/gempic_rz
$X emission.deck numerics.backend=openmp grid.dzeta=12.5nm output.dir=em_dzeta2 > em_dzeta2.log 2>&1
$X emission.deck numerics.backend=openmp grid.dzeta=12.5nm "grid.dr_points=0.5, 0.5, 6, 16, 25 um" output.dir=em_dzeta2_drc > em_dzeta2_drc.log 2>&1
$X emission.deck numerics.backend=openmp grid.dzeta=12.5nm "grid.dr_points=0.5, 0.5, 1.5, 4, 25 um" output.dir=em_dzeta2_drf > em_dzeta2_drf.log 2>&1
python3 compare_emission.py em_out em_dzeta2 em_dzeta2_drc em_dzeta2_drf > conv.txt 2>&1
