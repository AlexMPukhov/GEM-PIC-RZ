"""Fig. bubble closure (benchmarks/lwfa_fbpic): (a) on-axis electron density at the back of the first bubble after
100 um, cold plasma and with source smoothing a = 2 k0^-1; (b) on-axis E_z behind the bubble, cold and smoothed, two
steps, and FBPIC; (c) peak on-axis density at the closure against the propagation distance for Delta = 5 um.
usage (from paper/): python3 scripts/fig_closure.py  ->  figs/fig_closure.pdf"""
import glob
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

B = "../benchmarks/lwfa_fbpic/"
E0 = 4.0126e12  # m c omega0 / e at 0.8 um [V/m]


def axis(run, D):
    for fn in sorted(glob.glob(B + run + "/axis_*.txt")):
        z = float(open(fn).readline().split("z =")[1].split("mm")[0])
        if abs(z * 1000 - D) < 1e-6:
            return np.loadtxt(fn)
    raise SystemExit(f"{run}: no axis file at {D} um")


def history(run):
    out = []
    for fn in sorted(glob.glob(B + run + "/axis_*.txt")):
        z = float(open(fn).readline().split("z =")[1].split("mm")[0]) * 1000
        a = np.loadtxt(fn)
        m = a[:, 0] < 14.0  # behind the laser-driven sheath, up to the back of the first bubble
        out.append((z, a[m, 3].max()))
    return np.array(out)


plt.rcParams.update({"font.size": 9, "font.family": "serif", "mathtext.fontset": "cm"})
fig, axs = plt.subplots(3, 1, figsize=(3.4, 6.2), gridspec_kw={"height_ratios": [1.0, 1.1, 0.9]})
runs = [("gem_out", "#bbbbbb", "-", 1.0, r"cold, $\Delta = 2.5\,\mu$m"),
        ("gem_out_dz1.25", "#f0a35e", "-", 1.0, r"cold, $\Delta = 1.25\,\mu$m"),
        ("gem_out_a2", "#1f4e9c", "-", 1.0, r"$a = 0.25\,\mu$m, $\Delta = 2.5\,\mu$m"),
        ("gem_out_a2_dz1.25", "#6fa8dc", ":", 1.4, r"$a = 0.25\,\mu$m, $\Delta = 1.25\,\mu$m")]
for run, c, ls, lw, lab in runs:
    g = axis(run, 100)
    m = (g[:, 0] > 11) & (g[:, 0] < 16)
    axs[0].semilogy(g[m, 0], np.maximum(g[m, 3], 1e-2), ls, color=c, lw=lw, label=lab)
    axs[1].plot(g[:, 0], g[:, 2], ls, color=c, lw=lw)
f = np.loadtxt(B + "fbpic_D100.txt")
m = (f[:, 0] >= 0) & (f[:, 0] <= 16)
axs[1].plot(f[m, 0], f[m, 1] / E0, "--", color="k", lw=0.8, label="FBPIC")
axs[0].set_xlim(11, 16)
axs[0].set_ylim(0.05, 30)
axs[0].set_ylabel(r"$n_e(r=0)/n_0$")
axs[0].legend(fontsize=5.5, loc="upper left", frameon=False)
axs[1].set_xlim(0, 16)
axs[1].set_xlabel(r"$\zeta$ [$\mu$m]")
axs[1].set_ylabel(r"$E_z(r=0)$ $[mc\omega_0/e]$")
axs[1].legend(fontsize=6, loc="upper left", frameon=False)
for run, c, mk, lab in [("gem_out_dz5", "#d1495b", "o", "cold"), ("gem_out_a2_dz5", "#1f4e9c", "s", r"$a = 0.25\,\mu$m")]:
    h = history(run)
    k = h[:, 0] >= 20
    axs[2].semilogy(h[k, 0], h[k, 1], "-" + mk, color=c, ms=2.5, lw=0.9, label=lab + r", $\Delta = 5\,\mu$m")
axs[2].set_xlabel(r"propagation distance [$\mu$m]")
axs[2].set_ylabel(r"max $n_e(r=0)/n_0$")
axs[2].set_xlim(20, 100)
axs[2].legend(fontsize=6, loc="upper left", frameon=False)
for a, t in zip(axs, ["(a)", "(b)", "(c)"]):
    a.text(0.97, 0.86, t, transform=a.transAxes, ha="right")
fig.tight_layout(h_pad=0.3)
fig.savefig("figs/fig_closure.pdf")
