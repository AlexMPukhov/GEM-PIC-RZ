#!/bin/sh
# numerics.relaxation: the passes per slice at a coarse zeta step. A 50 pC driver in a uniform plasma, m_max = 1,
# theta 0.6, dz = 200 c/omega_p, dzeta = 1/16 c/omega_p (theta dz dzeta k_p^2 = 7.5, 1.25 times the limit), two steps:
#   - relaxation = 1 in the deck (4 passes): kept, the code warns, and the passes diverge at the axis (E_x > 1e-2);
#   - no relaxation in the deck: the code switches on relaxation 0.4 with 4 passes and says so; E_x on the axis stays
#     at the noise (< 1e-4) and max |E_z| on the axis agrees with dzeta = 1/32 c/omega_p within 5 %;
#   - dzeta = 1/32: nothing switched on, no warning;
#   - serial and the parallel backend give the same file.
# usage: sh tests/pass_relaxation.sh <gempic_rz> [backend]
set -e
abspath() { (cd "$(dirname "$1")" && echo "$(pwd)/$(basename "$1")"); }
EXE=$(abspath "$1"); BACKEND=${2:-openmp}
WORK=$(mktemp -d)
cd "$WORK"
cat > d.deck <<'EOF'
grid {
  lambda0       = 63.0986 um
  m_max         = 1
  radial_grid   = uniform
  dr            = 20.0849 um
  r_max         = 1606.79 um
  dzeta         = 12.5531 um
  window_length = 4016.97 um
  dz            = 40169.7 um
  z_end         = 50000 um
}
plasma { density = 7.0e14 cm^-3; r_max = 1205.09 um; ppc_r = 2; ppc_phi = 4 }
bunch {
  name = driver; charge = 50 pC; energy = 10 GeV; sigma_x = 200.849 um; emittance_x = 1 um
  sigma_z = 200.849 um; zeta_center = 3200 um; macroparticles = 10000
}
numerics { theta = 0.6; iterations = 4 }
output { dir = out; every = 2; format = text; axis = true; maps = false; particles = false }
EOF
O="numerics.backend=$BACKEND"
"$EXE" d.deck $O output.dir=a numerics.relaxation=1 > a.log
"$EXE" d.deck $O output.dir=b > b.log
"$EXE" d.deck $O output.dir=f grid.dzeta=6.27655um > f.log
"$EXE" d.deck numerics.backend=serial output.dir=s > s.log
grep -q "with numerics.relaxation = 1 and 4 passes they may diverge" a.log || { echo "  no warning with relaxation 1  FAILED"; exit 1; }
grep -q "4 iterations/slice (predictor zeta, relaxation 0.4 automatic)" b.log || { echo "  relaxation not switched on  FAILED"; exit 1; }
grep -q "The code switched on numerics.relaxation = 0.4" b.log || { echo "  no explanation  FAILED"; exit 1; }
if grep -q "limit 6 of the passes" f.log || grep -q "relaxation 0.4" f.log; then echo "  relaxation at dzeta / 2  FAILED"; exit 1; fi
mx() { awk -v c="$2" '!/^#/ { v = $c < 0 ? -$c : $c; if (v > m) m = v } END { printf "%.6e", m }' "$1"; }
EXA=$(mx a/axis_00002.txt 2); EXB=$(mx b/axis_00002.txt 2)
EZB=$(mx b/axis_00002.txt 3); EZF=$(mx f/axis_00002.txt 3)
echo "  max |E_x| on the axis: relaxation 1: $EXA, automatic: $EXB; max |E_z|: $EZB, dzeta / 2: $EZF"
awk -v x="$EXA" 'BEGIN { exit !(x > 1e-2) }' || { echo "  passes without relaxation did not diverge  FAILED"; exit 1; }
awk -v x="$EXB" 'BEGIN { exit !(x < 1e-4) }' || { echo "  passes with automatic relaxation diverged  FAILED"; exit 1; }
awk -v a="$EZB" -v b="$EZF" 'BEGIN { d = a / b - 1; if (d < 0) d = -d; exit !(d < 0.05) }' ||
  { echo "  wake differs from dzeta / 2 by more than 5 %  FAILED"; exit 1; }
cmp -s b/axis_00002.txt s/axis_00002.txt || { echo "  serial and $BACKEND differ  FAILED"; exit 1; }
echo "  ok"
cd /
rm -rf "$WORK"
