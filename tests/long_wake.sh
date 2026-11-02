#!/bin/sh
# Long linear plasma wake (numerics.predictor). tests/long_wake.deck: a weak Gaussian electron driver
# (n_b/n0 = 1.05e-5, k_p sigma_r = k_p sigma_z = 1) in a uniform plasma, window 60 c/omega_p, 8 steps of 100 c/omega_p.
# The on-axis E_z behind the driver is fitted with a cos(k_p zeta) + b sin(k_p zeta) over 5 ... 45 c/omega_p and
# compared with the linear wake (Keinigs & Jones, Phys. Fluids 30, 252 (1987)):
#   E_z / E0 = (n_b/n0) sqrt(2 pi) k_p sigma_z exp(-k_p^2 sigma_z^2 / 2) int_0^inf k_p^2 r K0(k_p r) exp(-r^2 / 2 sigma_r^2) dr,
# E0 = m c omega_p / e = (1/20) m c omega0 / e here (lambda0 = lambda_p / 20).
#   1) predictor = zeta (default): amplitude at step 8 within 3 % of the theory, steps 6 -> 8 within 3 %
#   2) predictor = level (the old one): the wake grows from step to step (> 10 x the theory at step 8)
# usage: sh tests/long_wake.sh <gempic_rz> [backend]
set -e
abspath() { (cd "$(dirname "$1")" && echo "$(pwd)/$(basename "$1")"); }
EXE=$(abspath "$1"); BACKEND=${2:-openmp}
SRC=$(cd "$(dirname "$0")/.." && pwd)
WORK=$(mktemp -d)
cd "$WORK"
# linear theory in m c omega0 / e: K0(r) = int_0^inf exp(-r cosh t) dt
THEORY=$(awk 'BEGIN {
  pi = 3.141592653589793; s = 0; dr = 0.002
  for (r = dr / 2; r < 12; r += dr) {
    k0 = 0; dt = 0.002
    for (t = dt / 2; t < 12; t += dt) { x = r * (exp(t) + exp(-t)) / 2; if (x > 700) break; k0 += exp(-x) * dt }
    s += r * k0 * exp(-r * r / 2) * dr
  }
  N = 0.15e-12 / 1.602176634e-19; sig = 200.849e-4                       # electrons; sigma in cm
  nb = N / ((2 * pi) ^ 1.5 * sig ^ 3) / 7.0e14
  printf "%.6e", nb * sqrt(2 * pi) * exp(-0.5) * s / 20 }')
amp() {  # <axis file>: amplitude of the k_p component of E_z on 5 ... 45 c/omega_p
  awk '!/^#/ { z = $1 / 200.8486; if (z < 5 || z > 45) next; c = cos(z); s = sin(z)
               cc += c * c; ss += s * s; cs += c * s; ce += c * $3; se += s * $3 }
       END { d = cc * ss - cs * cs; a = (ce * ss - se * cs) / d; b = (se * cc - ce * cs) / d; printf "%.6e", sqrt(a * a + b * b) }' "$1"
}
check() {  # <label> <value> <expected> <relative tolerance>
  if awk -v v="$2" -v e="$3" -v t="$4" 'BEGIN { d = (v - e) / e; if (d < 0) d = -d; exit !(d <= t) }'; then
    echo "  $1 = $2 (expected $3)  ok"
  else
    echo "  $1 = $2, expected $3 +- $4 (relative)  FAILED"; exit 1
  fi
}
"$EXE" "$SRC/tests/long_wake.deck" numerics.backend=$BACKEND > zeta.log
A6=$(amp wake/axis_00006.txt); A8=$(amp wake/axis_00008.txt)
check "zeta: E_z amplitude at step 8" "$A8" "$THEORY" 0.03
check "zeta: step 8 / step 6" "$(awk -v a="$A8" -v b="$A6" 'BEGIN { printf "%.4f", a / b }')" 1 0.03
"$EXE" "$SRC/tests/long_wake.deck" numerics.backend=$BACKEND numerics.predictor=level output.dir=lev > lev.log
AL=$(amp lev/axis_00008.txt)
if awk -v a="$AL" -v t="$THEORY" 'BEGIN { exit !(a > 10 * t) }'; then
  echo "  level: amplitude at step 8 = $AL (> 10 x theory: the instability is seen)  ok"
else
  echo "  level: amplitude at step 8 = $AL, expected > 10 x $THEORY  FAILED"; exit 1
fi
cd - > /dev/null
rm -rf "$WORK"
echo "long wake test passed ($BACKEND)"
