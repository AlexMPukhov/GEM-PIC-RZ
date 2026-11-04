#!/bin/sh
# bunch.azimuthal_copies: every sampled particle with K copies turned by 2 pi j / K about the bunch axis.
# A centred 50 pC driver in a uniform plasma (8 particles per ring), m_max = 1, two steps:
#   - K = 4: E_x on the axis (mode 1) stays at round-off (< 1e-9 of max |E_z|); with random loading it is the
#     macroparticle noise (> 1e-4 of max |E_z|);
#   - the wake (mode 0) agrees with the random loading (max |E_z| within 2 %), and so do the charge and sigma_x;
#   - with x0 = 30 um, y0 = -20 um the centroid of the particles is (x0, y0) to 1e-6 um.
# usage: sh tests/azimuthal_copies.sh <gempic_rz> [backend]
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
  dzeta         = 3.13826 um
  window_length = 4016.97 um
  dz            = 40169.7 um
  z_end         = 100000 um
}
plasma { density = 7.0e14 cm^-3; r_max = 1205.09 um; ppc_r = 2; ppc_phi = 8 }
bunch {
  name = driver; charge = 50 pC; energy = 10 GeV; sigma_x = 200.849 um; emittance_x = 1 um
  sigma_z = 200.849 um; zeta_center = 3200 um; macroparticles = 40000
}
numerics { theta = 0.6; iterations = 4 }
output { dir = out; every = 2; format = text; axis = true; maps = false; particles = false }
EOF
O="numerics.backend=$BACKEND"
"$EXE" d.deck $O output.dir=a > a.log
"$EXE" d.deck $O output.dir=b bunch.azimuthal_copies=4 > b.log
grep -q "azimuthal copies: 4" b.log || { echo "  no azimuthal copies reported  FAILED"; exit 1; }
mx() { awk -v c="$2" '!/^#/ { v = $c < 0 ? -$c : $c; if (v > m) m = v } END { printf "%.6e", m }' "$1"; }
EXA=$(mx a/axis_00002.txt 2); EZA=$(mx a/axis_00002.txt 3)
EXB=$(mx b/axis_00002.txt 2); EZB=$(mx b/axis_00002.txt 3)
echo "  max |E_x| on the axis: random $EXA, K = 4 $EXB; max |E_z|: $EZA, $EZB"
awk -v x="$EXA" -v z="$EZA" 'BEGIN { exit !(x > 1e-4 * z) }' || { echo "  random loading without m = 1 noise?  FAILED"; exit 1; }
awk -v x="$EXB" -v z="$EZB" 'BEGIN { exit !(x < 1e-9 * z) }' || { echo "  m = 1 not at round-off with K = 4  FAILED"; exit 1; }
awk -v a="$EZA" -v b="$EZB" 'BEGIN { d = b / a - 1; if (d < 0) d = -d; exit !(d < 0.02) }' ||
  { echo "  wake changed  FAILED"; exit 1; }
row() { awk '$1 == 0 && $3 == "driver"' "$1/bunches.txt"; }
QA=$(row a | awk '{print $5}'); QB=$(row b | awk '{print $5}')
SA=$(row a | awk '{print $8}'); SB=$(row b | awk '{print $8}')
NB=$(row b | awk '{print $4}')
echo "  macroparticles $NB, charge $QA / $QB pC, sigma_x $SA / $SB um"
[ "$NB" = 40000 ] || { echo "  macroparticles $NB, expected 40000  FAILED"; exit 1; }
awk -v a="$QA" -v b="$QB" 'BEGIN { d = b / a - 1; if (d < 0) d = -d; exit !(d < 1e-9) }' || { echo "  charge changed  FAILED"; exit 1; }
awk -v a="$SA" -v b="$SB" 'BEGIN { d = b / a - 1; if (d < 0) d = -d; exit !(d < 0.02) }' || { echo "  sigma_x changed  FAILED"; exit 1; }
"$EXE" d.deck $O output.dir=c bunch.azimuthal_copies=4 bunch.x0=30um bunch.y0=-20um grid.z_end=1um \
  output.particles=true output.every=1 > c.log
F=$(ls c/beam_*.txt | head -n 1)
C=$(awk '!/^#/ && $8 == 0 { x += $1; y += $2; n++ } END { printf "%.9f %.9f", x / n, y / n }' "$F")
echo "  centroid with x0 = 30 um, y0 = -20 um: $C um"
echo "$C" | awk '{ dx = $1 - 30; dy = $2 + 20; exit !(dx * dx + dy * dy < 1e-12) }' || { echo "  centroid  FAILED"; exit 1; }
echo "  ok"
cd /
rm -rf "$WORK"
