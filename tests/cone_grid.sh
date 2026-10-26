#!/bin/sh
# --check warning for undulator harmonic cones that the radial grid cannot carry (config.cpp, coneGridWarnings).
# At h omega_1 the harmonic h + 1 is emitted on the cone theta^2 = (1/h)(1 + K^2/2)/gamma^2, k_perp = h k_r theta.
# For benchmarks/harmonics/fel.deck (K = 1.5, gamma = 1000, period 1 cm, lambda_r = 10.625 nm) at 2 omega_1:
# theta = 1.03 mrad, k_perp = 1.22 /um - the 5.1 um structure of the m = 2 field near the beam, whose wavenumber
# was measured in the full run as 1.31 /um on the dr = 1 um grid (= (2/dr) asin(k_perp dr/2)) with K = -k_u.
# Checks:
#  1. theta and k_perp of the warnings for h = 1, 2, 3 against the same formula evaluated here (awk) from the deck;
#  2. the deck's table grid (dr = 1 um to 40 um, then stretching): the h = 2 cone is turned back where the grid
#     stretches (r = 40 um) and thrown back at r = 51.5 um, where the step first exceeds 2/k_perp = 1.64 um;
#  3. dr = 1 um out to 100 um moves the cutoff of the h = 2 cone beyond 100 um;
#  4. a uniform grid of 0.5 um (k_perp dr = 0.61 for h = 2, 0.75 for h = 3) gives no cone warning.
# usage: sh tests/cone_grid.sh <gempic_rz>
set -e
EXE=$1
SRC=$(cd "$(dirname "$0")/.." && pwd)
DECK=$SRC/benchmarks/harmonics/fel.deck
WORK=$(mktemp -d)
cd "$WORK"
fail=0
"$EXE" "$DECK" --check > base.txt 2>&1 || true
grep 'on a cone' base.txt > cones.txt || { echo "FAIL: no cone warning for fel.deck"; cat base.txt; exit 1; }
# 1. theta [mrad] and k_perp [/um] from K, gamma, period, lambda_r (deck values), h = 1, 2, 3
for h in 1 2 3; do
  line=$(grep "at $h omega_1 on a cone" cones.txt) || { echo "FAIL: no warning for h = $h"; fail=1; continue; }
  th=$(echo "$line" | sed 's/.*theta = \([0-9.e+-]*\) mrad.*/\1/')
  kp=$(echo "$line" | sed 's/.*k_perp = \([0-9.e+-]*\) \/um.*/\1/')
  awk -v h=$h -v th="$th" -v kp="$kp" 'BEGIN {
    K = 1.5; g = 1000; lu = 1e4; pi = 3.14159265358979
    lr = lu * (1 + K * K / 2) / (2 * g * g)                     # um
    t = sqrt((1 / h) * (1 + K * K / 2)) / g                        # rad
    k = h * 2 * pi / lr * t                                        # /um
    ok = (sqrt((th / (1e3 * t) - 1)^2) < 1e-3) && (sqrt((kp / k - 1)^2) < 1e-3)
    printf "h = %d: theta %.4g mrad (expected %.4g), k_perp %.4g /um (expected %.4g, k_perp^2 = 2 h k_r k_u: %.4g)  %s\n",
           h, th, t * 1e3, kp, k, sqrt(2 * h * (2 * pi / lr) * (2 * pi / lu)), ok ? "ok" : "FAIL"
    exit ok ? 0 : 1 }' || fail=1
done
# 2. reflection radii of the h = 2 cone on the deck grid
line=$(grep "at 2 omega_1 on a cone" cones.txt)
echo "$line" | grep -q 'stretches at r = 40 um' && echo "h = 2: turned back where the grid stretches (40 um)  ok" ||
  { echo "FAIL: stretching point"; echo "$line"; fail=1; }
echo "$line" | grep -q 'at r = 51.5[0-9]* um the radial step 1.6[0-9]* um exceeds 2/k_perp = 1.64 um' &&
  echo "h = 2: thrown back at r = 51.5 um (dr > 1.64 um)  ok" || { echo "FAIL: cutoff radius"; echo "$line"; fail=1; }
# 3. dr = 1 um out to 100 um
"$EXE" "$DECK" grid.r_points=0,100,400,1000um --check 2>&1 | grep "at 2 omega_1 on a cone" > g100.txt || true
r=$(sed 's/.*the grid stretches at r = \([0-9.]*\) um.*/\1/' g100.txt)
awk -v r="$r" 'BEGIN { ok = (r + 0 >= 100); printf "dr = 1 um to 100 um: the h = 2 cone is turned back at r = %s um  %s\n", r, ok ? "ok" : "FAIL"; exit ok ? 0 : 1 }' || fail=1
# 4. uniform 0.5 um: no warning
if "$EXE" "$DECK" grid.dr_points=0.5,0.5,0.5,0.5um --check 2>&1 | grep -q 'on a cone'; then
  echo "FAIL: cone warning on a uniform 0.5 um grid"; fail=1
else
  echo "uniform 0.5 um grid: no cone warning  ok"
fi
cd - > /dev/null
rm -rf "$WORK"
[ $fail -eq 0 ] && echo "cone-grid test passed" || { echo "cone-grid test FAILED"; exit 1; }
