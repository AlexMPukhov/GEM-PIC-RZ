# LWFA cases of the QuickPIC laser-module paper (Meng et al. 2026)

W. Meng et al., "Improvement of the laser module in the quasi-static PIC code QuickPIC", Plasma Phys. Control.
Fusion (2026), doi:10.1088/1361-6587/aebd1b, Sec. 4.1: two laser wakefield cases benchmarked against OSIRIS.
The same cases run in QUARZ (envelope code, `benchmarks/quickpic_laser` of AlexMPukhov/QUARZ, which also holds the
FBPIC reference script `fbpic_ref.py` and the comparison script). Here they run in GEM-PIC-RZ, which resolves the
laser period (Δζ = λ0/25) and keeps the slow-time derivative of the forward field.

| deck | case |
|---|---|
| `case1.deck` | a0 = 2.5, k_p w0 = 3.162, k0/k_p = 34 (0.8 µm, 1.5·10¹⁸ cm⁻³), ω_p τ_FWHM = 2.26 (FWHM of the field envelope), z = 0 … 500 c/ω_p (2.16 mm) |
| `case2.deck` | a0 = 4, k_p w0 = 4, k0/k_p = 20 (4.4·10¹⁸ cm⁻³), same pulse length, z = 0 … 200 c/ω_p (0.51 mm) |
| `run_gem.sbatch` | Slurm runner (16 OpenMP threads) |

Grid as in QuickPIC/QUARZ: Δr = 0.1 c/ω_p, window ξ = −5 … 12 c/ω_p around the initial pulse centre, slow step
Δz = 0.5 c/ω_p; plasma electrons cold, ions immobile, no smoothing.

Results (normalised rms error of the on-axis E_z, rms(ΔE_z)/max|E_z,ref|; FBPIC 0.27 with the same parameters as
the reference):

| case, z [c/ω_p] | GEM-PIC-RZ vs FBPIC | QUARZ vs FBPIC |
|---|---|---|
| 1, 50 | 0.44 % in front of the bubble closure (ξ < 5.5), 10.3 % with the closure | 0.73 %, 5.5 % |
| 1, 100 / 200 / 300 / 400 / 500 (front; whole box) | 0.49 / 0.51 / 0.72 / 0.86 / 0.93 %; 17.3 / 5.0 / 2.6 / 1.2 / 0.9 % | 0.89 / 1.04 / 1.23 / 1.39 / 1.44 %; 13.9 / 3.8 / 2.9 / 1.7 / 1.5 % |
| 2, 20 / 50 / 100 / 150 / 200 (ξ < 7) | 1.4 / 1.4 / 1.1 / 1.5 / 1.6 % | 1.0 / 1.1 / 1.3 / 1.5 / 1.6 % |

The laser envelope on the axis (max|E_x|/E_L for |ξ| < 3) agrees with FBPIC to 0.7–1.6 % in case 1 (2.60, 2.68, 2.34, 1.88, 1.53, 1.26 vs 2.58, 2.66, 2.31, 1.85, 1.50, 1.25 at z = 50 … 500) and to 1–3 % in case 2 (4.03, 4.36, 5.13, 5.26, 5.04 vs
4.02, 4.34, 5.07, 5.20, 4.89). At the bubble closure GEM-PIC-RZ has a sharper E_z spike and traps electrons
(compare README, "Bubble closure"; `plasma.smooth_length` regularises it). Cost on 16 threads of an AMD EPYC 7543:
case 1, 951 s (4.2 core-hours); case 2, 442 s (2.0 core-hours).
