#!/bin/sh
# Convergence scan for the flying-focus decks (run on the cluster; adapt the launcher).
# Each run changes one numerical parameter of ff_plasma.deck over the first 2 mm. Compare with
#   python3 ../check_axis.py <dir> 2 0.05 50 95      (laser on the axis vs the ideal pattern)
#   python3 ../check_wake.py <dir> 95                (wake amplitude and position)
X=${GEMPIC:-gempic_rz}
B="ff_plasma.deck grid.z_end=2mm"
$X $B output.dir=conv_ref
$X $B grid.dr=0.25um               output.dir=conv_dr
$X $B grid.dzeta=0.0167um          output.dir=conv_dzeta
$X $B grid.dz=1.25um               output.dir=conv_dz
$X $B grid.r_max=200um grid.r_absorber=160um output.dir=conv_rmax
$X $B plasma.ppc_r=4 plasma.ppc_phi=16 output.dir=conv_ppc
