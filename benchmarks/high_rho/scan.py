"""rho scan with fel1d.py: averaged vs unaveraged 1D steady-state planar FEL (gamma 500, 5 kA, lambda_u 5 mm;
rho set through the beam area). Writes scan.txt."""
import numpy as np
import fel1d as F

lu, g0, I = 5e-3, 500.0, 5000.0
A0 = 2 * np.pi * (10e-6) ** 2
rows = []
for K in (1.0, 2.0, 4.0):
    for rt in (0.002, 0.005, 0.01, 0.02, 0.04, 0.08, 0.15):
        A = A0 * (F.rho_1d(g0, I, A0, K, lu) / rt) ** 3
        rho = F.rho_1d(g0, I, A, K, lu)
        Lg1 = lu / (4 * np.pi * np.sqrt(3) * rho)
        Ps = rho * g0 * 0.511e6 * I
        r = []
        for m in ("averaged", "unaveraged"):
            z, P, b = F.run(g0, I, A, K, lu, 25 * Lg1, 1e-7 * Ps, m)
            r += [F.gain_length(z, P) / Lg1, P.max() / Ps, z[np.argmax(P)] / lu]
        rows.append([K, rho, Lg1 / lu] + r)
        print(*(f"{v:.4g}" for v in rows[-1]))
np.savetxt("scan.txt", rows, fmt="%.5g", header="K rho Lg1D/lambda_u | averaged: Lg/Lg1D Psat/(rho P_beam) "
           "z_sat/lambda_u | unaveraged: the same   (npart 256, lambda_u/96; converged to 1e-4 at rho 0.04)")
