#!/bin/sh
# Profiles given as functions of the deck (function parser) reproduce the
# built-in profiles of tests/regression.deck:
#   plasma   z_points / z_values ramp          -> density_function(z)
#   laser    waist, length_1e (Gaussian focus) -> profile_function(r), envelope_function(zeta)
#   bunch    flat-top, gamma                   -> distribution_z = function, profile_z, gamma_function
#   numbers  given through the block 'constants'
# Both runs draw the same random numbers: the results agree to rounding.
# Undulator taper K_function(z) = 0.5 + 0.3 z / L in tests/undulator.deck: the
# slippage of the bunch int (1 + K(z)^2 g(z)^2 / 2) / (2 gamma^2) dz (g: ramps)
# is 8.784 um (+- 1 %; untapered: 8.235 um).
# usage: sh tests/functions.sh <gempic_rz> <compare_columns> [backend]
set -e
abspath() { (cd "$(dirname "$1")" && echo "$(pwd)/$(basename "$1")"); }
EXE=$(abspath "$1"); CMP=$(abspath "$2"); BACKEND=${3:-serial}
SRC=$(cd "$(dirname "$0")/.." && pwd)
WORK=$(mktemp -d)
cd "$WORK"
# reference: the regression deck with a flat-top bunch
sed -e 's/^  sigma_z .*/  distribution_z = flattop\n  length         = 8 k0^-1/' "$SRC/tests/regression.deck" > ref.deck
# the same with functions
sed -e 's/^  z_points .*/  density_function = "min(max(z \/ ramp, 0), 1)"/' -e '/^  z_values /d' \
    -e 's/^  waist .*/  profile_function  = "exp(-(r \/ w0)^2)"/' \
    -e 's/^  length_1e .*/  envelope_function = "gauss(zeta - zl, len \/ sqrt(2))"/' \
    -e 's/^  zeta_center *= 240 k0^-1/  zeta_center = zl/' \
    -e 's/^  gamma .*/  gamma_function = "1000"/' \
    -e 's/^  sigma_z .*/  distribution_z = function\n  profile_z      = "rect(zeta - 130 k0^-1, 8 k0^-1)"\n  zeta_min       = 126 k0^-1\n  zeta_max       = 134 k0^-1/' \
    -e '/^  zeta_center *= 130 k0^-1/d' \
    "$SRC/tests/regression.deck" > fn.deck
cat >> fn.deck <<'END'
constants {
  ramp = 400 k0^-1
  w0   = 63 k0^-1
  zl   = 240 k0^-1
  len  = 25 k0^-1
}
END
grep -q density_function fn.deck && grep -q envelope_function fn.deck && grep -q profile_z fn.deck
O="numerics.backend=$BACKEND grid.z_end=80k0^-1 output.every=1 output.map_stride=8"
"$EXE" ref.deck $O output.dir=ref > ref.log
"$EXE" fn.deck $O output.dir=fn > fn.log
for f in axis_00000.txt axis_00002.txt map_00002.txt bunches.txt; do "$CMP" ref/$f fn/$f 1e-7 1e-10; done
# undulator taper
sed 's/K = 0.5;/K_function = "0.5 + 0.3 * z \/ (142.222 mm)";/' "$SRC/tests/undulator.deck" > taper.deck
grep -q K_function taper.deck
"$EXE" taper.deck numerics.backend=$BACKEND bunch.macroparticles=5000 output.format=text output.maps=false \
  output.axis=false output.particles=false output.dir=taper > taper.log
S=$(awk '!/^#/ && $3 == "b" { if (!seen) z0 = $12; seen = 1; z1 = $12 } END { printf "%.3f", z0 - z1 }' taper/bunches.txt)
if awk -v s="$S" 'BEGIN { d = s - 8.784; if (d < 0) d = -d; exit !(d <= 0.088) }'; then
  echo "  slippage in the tapered undulator $S um (expected 8.784 +- 0.088)  ok"
else
  echo "  slippage in the tapered undulator $S um, expected 8.784 +- 0.088  FAILED"; exit 1
fi
cd - > /dev/null
[ -n "$KEEP" ] && echo "$WORK" || rm -rf "$WORK"
echo "functions test passed ($BACKEND)"
