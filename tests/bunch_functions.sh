#!/bin/sh
# bunch.distribution_r = function with profile_r, ux_function, uy_function: a helical ring (r = 100 um, rms width
# 1 um, rigid rotation u_phi = K r / R0, K = 0.5, gamma = 50) in a plasma lens matched to its circular orbits
# (G = K^2 m c / (e gamma R0^2) = 852 T/m). Checked: the sampled ring (mean radius, width, u_phi, gamma exact:
# the momentum functions keep the energy) and that it stays on its orbits over 1 cm (1/6 of a betatron period):
# mean radius within 0.1 %, width within 2 %. A ring without rotation would fall to 0.54 R0.
# usage: sh tests/bunch_functions.sh <gempic_rz> [backend]
set -e
abspath() { (cd "$(dirname "$1")" && echo "$(pwd)/$(basename "$1")"); }
EXE=$(abspath "$1"); BACKEND=${2:-openmp}
WORK=$(mktemp -d)
cd "$WORK"
cat > ring.deck <<'DECK'
grid {
  lambda0 = 1 um; m_max = 1; radial_grid = uniform; dr = 5 um; r_max = 400 um
  dzeta = 0.5 um; window_length = 40 um; dz = 1 mm; z_end = 10 mm; r_absorber = 500 um
}
constants { R0 = 100 um; K = 0.5 }
bunch {
  charge = 1 fC; gamma = 50; distribution_z = flattop; length = 10 um; zeta_center = 20 um
  distribution_r = function; profile_r = "gauss(r - R0, 1 um)"; r_extent = 120 um
  ux_function = "-K * y / R0"; uy_function = "K * x / R0"
  macroparticles = 20000
}
plasma_lens { z_start = 0 mm; length = 20 mm; gradient = 852 T/m }
output { dir = o; every = 10; maps = false }
numerics { beam_dz_max = 0.02 mm }
DECK
"$EXE" ring.deck numerics.backend="$BACKEND" > run.log
grep -q "profile_r(x, y, r, zeta) = gauss" run.log || { echo "  profile_r not reported  FAILED"; exit 1; }
mom() {  # mean r, rms width of r, mean u_phi, max |gamma - 50|
  awk '$1 !~ /^#/ { r = sqrt($1*$1 + $2*$2); n++; s += r; s2 += r*r; up += ($1*$5 - $2*$4) / r
                    g = sqrt(1 + $4*$4 + $5*$5 + $6*$6) - 50; if (g < 0) g = -g; if (g > gm) gm = g }
       END { m = s / n; printf "%.5f %.5f %.6f %.3g\n", m, sqrt(s2 / n - m*m), up / n, gm }' "$1"
}
A=$(mom o/beam_00000.txt); B=$(mom o/beam_00010.txt)
echo "  start: r = $(echo $A | cut -d' ' -f1) um, width $(echo $A | cut -d' ' -f2) um, u_phi $(echo $A | cut -d' ' -f3), |dgamma| $(echo $A | cut -d' ' -f4)"
echo "  1 cm : r = $(echo $B | cut -d' ' -f1) um, width $(echo $B | cut -d' ' -f2) um, u_phi $(echo $B | cut -d' ' -f3)"
echo "$A" | awk '{ exit !($1 > 99.9 && $1 < 100.1 && $2 > 0.98 && $2 < 1.02 && $3 > 0.499 && $3 < 0.501 && $4 < 1e-4) }' ||
  { echo "  sampled ring wrong  FAILED"; exit 1; }
echo "$B" | awk '{ exit !($1 > 99.9 && $1 < 100.1 && $2 > 0.98 && $2 < 1.02) }' || { echo "  ring left its orbits  FAILED"; exit 1; }
echo "  ok"
cd /
rm -rf "$WORK"
