"""Sensitivity of the FBPIC benchmarks to the source-smoothing length a (plasma.smooth_length, .smooth_length_m0) and
to the step Delta. Reads the GEM-PIC-RZ runs in lwfa_fbpic/, laser_evolution/, ionization_fbpic/, downramp_fbpic/
(see README.md, section "Sensitivity to the smoothing length") and prints one table per benchmark, the numbers of
Table "sensitivity" of the paper.
usage (in benchmarks/): python3 sensitivity.py  ->  stdout, sensitivity.txt"""
import glob
import os
import subprocess
import numpy as np

E0 = 4.0126e12  # m c omega0 / e at 0.8 um [V/m]
out = []


def say(s=""):
    print(s)
    out.append(s)


def axis_at(run, D):
    for fn in sorted(glob.glob(f"{run}/axis_*.txt")):
        z = float(open(fn).readline().split("z =")[1].split("mm")[0])
        if abs(z * 1000 - D) < 1e-6:
            return np.loadtxt(fn)
    return None


# --- LWFA wake (lwfa_fbpic): first wake period, wake behind the bubble, E_z at the bubble back, at 100 um
os.chdir("lwfa_fbpic")
f = np.loadtxt("fbpic_D100.txt")
zf, Ezf = f[:, 0], f[:, 1] / E0
say("LWFA wake at 100 um (a0 = 2, 0.002 n_c): rms(E_z - FBPIC) in the first period [% of peak], behind the bubble "
    "(zeta 1-12.5 um) [1e-3 m c w0/e], max E_z at the bubble back (FBPIC 0.042)")
for run, lab in [("gem_out", "a = 0"), ("gem_out_a1", "a = 1"), ("gem_out_a2", "a = 2"), ("gem_out_a4", "a = 4")]:
    row = f"  {lab:8s}"
    for sfx, dz in [("", 2.5), ("_dz1.25", 1.25)]:
        r = run + sfx if run != "gem_out" else ("gem_out" if not sfx else "gem_out_dz1.25")
        g = axis_at(r, 100)
        if g is None:
            row += f"   dz {dz}: --"
            continue
        z, Ez = g[:, 0], g[:, 2]
        m1 = (z > 14) & (z < 26)
        E1 = np.interp(z[m1], zf, Ezf)
        mb = (z > 1) & (z < 12.5)
        rb = np.sqrt(np.mean((Ez[mb] - np.interp(z[mb], zf, Ezf)) ** 2))
        pk = Ez[(z > 10) & (z < 14)].max()
        row += (f"   dz {dz}: {100 * np.sqrt(np.mean((Ez[m1] - E1) ** 2)) / abs(E1).max():4.1f} %"
                f" {1e3 * rb:5.2f} {pk:.3f}")
    say(row)
os.chdir("..")

# --- laser evolution (laser_evolution): W/W0, <k>/k0 and the wake at 300 um
os.chdir("laser_evolution")
H5 = "../../build/h5_extract"


def laser_table(d):
    cache = f"laser_{d}.txt"
    files = sorted(glob.glob(f"{d}/data_*.h5"))
    if os.path.exists(cache) and (not files or os.path.getmtime(cache) > max(os.path.getmtime(x) for x in files)):
        return np.loadtxt(cache)
    rows = []
    for fn in files:
        e = subprocess.run([H5, fn, "energy"], capture_output=True, text=True, check=True).stdout.split("\n")[1].split()
        sp = np.array([l.split() for l in subprocess.run([H5, fn, "spectrum", "-1e9", "1e9", "0.3", "1.6", "261"],
                                                          capture_output=True, text=True, check=True).stdout.split("\n")
                       if l and not l.startswith("#")], float)
        k, S = sp[:, 0], sp[:, 2]
        kbar = np.sum(k * S) / np.sum(S)
        step = int(fn.split("_")[-1].split(".")[0])
        ax = np.loadtxt(f"{d}/axis_{step:05d}.txt")
        rows.append((float(e[1]) * 1000, float(e[5]), np.abs(ax[:, 1]).max(), kbar,
                     np.sqrt(np.sum((k - kbar) ** 2 * S) / np.sum(S))))
    t = np.array(rows)
    np.savetxt(cache, t, header="z[um]  W(m>=1)[code]  max|Ex|_axis[m c w0/e]  <|k|>/k0  sigma_k/k0")
    return t


fb = np.loadtxt("fbpic_out/laser_fbpic.txt")
fa = np.loadtxt(sorted(glob.glob("fbpic_out/axis_*.txt"))[-1])
say()
say(f"Laser evolution at 300 um (a0 = 1.5, 0.01 n_c): W/W0 (FBPIC {fb[-1, 1] / fb[0, 1]:.3f}), <k>/k0 (FBPIC "
    f"{fb[-1, 3]:.4f}), rms(E_z - FBPIC) on the axis [% of the FBPIC peak]")
for d, lab in [("gem_out", "cold"), ("gem_out_a0.5m3", "a = 0.5 / 3"), ("gem_out_a1m2", "a = 1 / 2"),
               ("gem_out_a1m3", "a = 1 / 3"), ("gem_out_a1m4", "a = 1 / 4"), ("gem_out_a2m3", "a = 2 / 3"), ("gem_out_a3m3", "a = 3 / 3"), ("gem_out_a1m1", "a = 1 / 1"),
               ("gem_out_a1m3_dz1", "a = 1 / 3, dz 1")]:
    if not glob.glob(f"{d}/data_*.h5") and not os.path.exists(f"laser_{d}.txt"):
        say(f"  {lab:16s} --")
        continue
    t = laser_table(d)
    last = sorted(glob.glob(f"{d}/axis_*.txt"))[-1]
    g = np.loadtxt(last)
    m = (fa[:, 0] >= 0) & (fa[:, 0] <= 24.95)
    Ezi = np.interp(g[:, 0], fa[m, 0], fa[m, 1])
    rms = 100 * np.sqrt(np.mean((g[:, 2] - Ezi) ** 2)) / np.abs(fa[m, 1]).max()
    say(f"  {lab:16s} z {t[-1, 0]:5.0f}  W/W0 {t[-1, 1] / t[0, 1]:.3f}  <k>/k0 {t[-1, 3]:.4f}  wake rms {rms:6.1f} %")
os.chdir("..")


# --- injection (ionization_fbpic, downramp_fbpic)
def charge(run, species, z, cut):
    """charge [pC] and mean energy [MeV] above `cut` of the trapped electrons at the level s = z [um]"""
    for fn in sorted(glob.glob(f"{run}/beam_*.txt*")):
        zz = float(open(fn.replace("beam_", "axis_").replace(".gz", "")).readline().split("z =")[1].split("mm")[0])
        if abs(zz * 1000 - z) < 1e-6:
            a = np.loadtxt(fn, ndmin=2)
            a = a[a[:, 7] == species] if a.size else np.zeros((0, 8))
            E = (np.sqrt(1 + a[:, 3] ** 2 + a[:, 4] ** 2 + a[:, 5] ** 2) - 1) * 0.511
            q = np.abs(a[:, 6])
            m = E > cut
            return q[m].sum(), (np.sum(q[m] * E[m]) / q[m].sum() if q[m].sum() > 0 else 0)
    return None


os.chdir("ionization_fbpic")
fb = np.loadtxt("injection_fbpic.txt")
k = ~np.isnan(fb[:, 1])
Qf, Ef = np.interp(175, fb[k, 1], fb[k, 5]), np.interp(175, fb[k, 1], fb[k, 6])
say()
say(f"Ionization injection, bunch at z = 175 um: Q(E > 5 MeV) [pC], <E> [MeV] (FBPIC, interpolated: {Qf:.2f}, {Ef:.1f})")
for run, lab in [("gem_out_a0", "a = 0"), ("gem_out_a1", "a = 1"), ("gem_out", "a = 2"), ("gem_out_a4", "a = 4"),
                 ("gem_out_a1_dz1.25", "a = 1, dz 1.25"), ("gem_out_dz1.25", "a = 2, dz 1.25")]:
    r = charge(run, -30, 175, 5.0)
    say(f"  {lab:16s} " + ("--" if r is None else f"{r[0]:5.2f} ({100 * (r[0] / Qf - 1):+4.0f} %)  {r[1]:5.1f}"))
os.chdir("..")

os.chdir("downramp_fbpic")
fb = np.loadtxt("injection_fbpic.txt")
say()
say(f"Down-ramp injection at s = 150 um: Q(E > 2 MeV) [pC] (FBPIC at ct = 150 um: {fb[fb[:, 0] == 150][0, 2]:.1f})")
for runs, lab in [(("gem_out_a0", "gem_out_a0_dz1.25"), "a = 0"), (("gem_out_a1", None), "a = 1"),
                  (("gem_out", "gem_out_dz1.25"), "a = 2"), (("gem_out_a4", "gem_out_a4_dz1.25"), "a = 4")]:
    row = f"  {lab:8s}"
    for run, dz in zip(runs, (2.5, 1.25)):
        r = charge(run, -1, 150, 2.0) if run else None
        row += f"   dz {dz}: " + ("--" if r is None else f"{r[0]:5.1f}")
    say(row)
os.chdir("..")
open("sensitivity.txt", "w").write("\n".join(out) + "\n")
