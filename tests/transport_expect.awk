# expected size ratio x(zend)/x(0) of a cold beam for tests/transport.sh (RK4 of x'' = -K(z) x)
# x'' = -K(z) x from z = 0 to zend for x(0) = 1, x'(0) = 0; K(z) given by kind
function K(z) {
  if (kind == "quad_x") return (z >= z0 && z < z1) ? k : 0;
  if (kind == "quad_y") return (z >= z0 && z < z1) ? -k : 0;
  if (kind == "sol") {  # Larmor frame: (B/(2 B rho))^2 with linear fringes
    f = 0; if (z >= z0 && z < z1) { f = 1; if (z - z0 < ramp) f = (z - z0) / ramp; if (z1 - z < ramp) f = (z1 - z) / ramp }
    return (f * k) * (f * k)
  }
}
BEGIN {
  x = 1; v = 0; n = 200000; h = zend / n
  for (i = 0; i < n; i++) {
    z = i * h
    a1 = -K(z) * x; x1 = x + 0.5 * h * v; v1 = v + 0.5 * h * a1
    a2 = -K(z + 0.5 * h) * x1; x2 = x + 0.5 * h * v1; v2 = v + 0.5 * h * a2
    a3 = -K(z + 0.5 * h) * x2; x3 = x + h * v2; v3 = v + h * a3
    a4 = -K(z + h) * x3
    x += h / 6 * (v + 2 * v1 + 2 * v2 + v3); v += h / 6 * (a1 + 2 * a2 + 2 * a3 + a4)
  }
  printf "%.5f\n", (x < 0 ? -x : x)
}
