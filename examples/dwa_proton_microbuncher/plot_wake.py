#!/usr/bin/env python3
"""Plots of a run of the DWA decks (text output of GEM-PIC-RZ: axis_*.txt, map_*.txt, bunches.txt, deck_used.txt).

    python3 plot_wake.py run                 # all output steps
    python3 plot_wake.py run --steps 50 250  # selected steps

Writes into the run directory:
    wake_axis.png    on-axis E_z(xi) for the steps; xi = zeta - zeta_driver (driver centre at 0, wake at xi < 0;
                     head of the window to the right); dashed: +-127 MV/m, the TM01 modal theory of the lining
    wake_map.png     E_z(xi, x) in the plane phi = 0 / pi at the last step (all modes summed)
    beam.png         driver: mean energy, rms energy spread, sigma_x and the normalized emittance vs z

Units: the text files hold fields in m c omega0 / e with omega0 = 2 pi c / lambda0 (grid.lambda0 from deck_used.txt).
Note: steps with z beyond the end of the dielectric (dielectric.z_end) show the field in the empty copper pipe, not
the wake of the structure (see README, "Structure exit").
"""
import argparse
import glob
import os
import re

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

MC2_EV = 0.51099895e6
WAKE_THEORY = 127.0  # MV/m, TM01 on axis behind 5 nC, sigma_z = 0.1 mm (README)


def deck_value(text, block, key):
    """first 'key = value' inside 'block { ... }' of deck_used.txt (value with unit as written), or None"""
    for m in re.finditer(block + r"\s*\{([^}]*)\}", text):
        k = re.search(r"\b" + key + r"\s*=\s*([-+0-9.eE]+)\s*([a-zA-Z]*)", m.group(1))
        if k:
            return float(k.group(1)), k.group(2)
    return None


def to_m(v):
    scale = {"m": 1.0, "mm": 1e-3, "um": 1e-6, "cm": 1e-2, "nm": 1e-9, "": 1.0}
    return v[0] * scale[v[1]]


def read_bunches(path, name):
    rows = []
    for line in open(path):
        if line.startswith("#"):
            continue
        f = line.split()
        if len(f) >= 13 and f[2] == name:
            rows.append([float(x) for x in (f[0], f[1], f[5], f[6], f[7], f[9], f[11])])
    return np.array(rows)  # step, z[mm], E[MeV], dE/E, sig_x[um], eps_nx[um], <zeta>[um]


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("run", help="output directory of the run (output.dir)")
    p.add_argument("--steps", type=int, nargs="*", help="output steps to plot (default: all)")
    p.add_argument("--bunch", default="driver", help="drive bunch name in bunches.txt (driver)")
    a = p.parse_args()

    deck = open(os.path.join(a.run, "deck_used.txt")).read()
    lam0 = to_m(deck_value(deck, "grid", "lambda0"))
    e0 = 2 * np.pi * MC2_EV / lam0 / 1e6  # MV/m per code unit
    zd = deck_value(deck, "dielectric", "z_end")
    z_diel_end = to_m(zd) * 1e3 if zd else None  # mm
    b = read_bunches(os.path.join(a.run, "bunches.txt"), a.bunch)
    zeta_of = {int(r[0]): r[6] / 1e3 for r in b}  # driver centre [mm] per step

    steps = sorted(int(re.search(r"(\d+)", os.path.basename(f)).group(1))
                   for f in glob.glob(os.path.join(a.run, "axis_*.txt")))
    if a.steps:
        steps = [s for s in steps if s in a.steps]
    steps = [s for s in steps if s > 0]
    if not steps:
        raise SystemExit("no axis_*.txt output steps found")

    # on-axis wake
    fig, ax = plt.subplots(figsize=(10, 4.2))
    cmap = plt.get_cmap("viridis")
    for j, s in enumerate(steps):
        f = os.path.join(a.run, "axis_%05d.txt" % s)
        z = float(re.search(r"z = ([-0-9.eE+]+) mm", open(f).readline()).group(1))
        d = np.loadtxt(f)
        xi = d[:, 0] / 1e3 - zeta_of.get(s, zeta_of[min(zeta_of)])
        after = z_diel_end is not None and z > z_diel_end + 1e-9
        lab = "z = %.0f mm%s" % (z, "  (%.0f mm after the structure end)" % (z - z_diel_end) if after else "")
        ax.plot(xi, d[:, 2] * e0, lw=1.0, ls="--" if after else "-",
                color="tab:red" if after else cmap(j / max(len(steps) - 1, 1)), label=lab)
    for v in (WAKE_THEORY, -WAKE_THEORY):
        ax.axhline(v, color="k", ls=":", lw=0.8)
    ax.text(0.01, 0.97, "dotted: TM01 modal theory, %.0f MV/m" % WAKE_THEORY, transform=ax.transAxes, va="top",
            fontsize=8)
    ax.set_xlabel(r"$\xi = \zeta - \zeta_{driver}$  [mm]   (driver centre at 0, wake at $\xi$ < 0)")
    ax.set_ylabel(r"$E_z$ on axis  [MV/m]")
    ax.legend(fontsize=7, loc="lower left", ncol=2)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(a.run, "wake_axis.png"), dpi=130)

    # map at the last step
    s = steps[-1]
    f = os.path.join(a.run, "map_%05d.txt" % s)
    if os.path.exists(f):
        z = float(re.search(r"z = ([-0-9.eE+]+) mm", open(f).readline()).group(1))
        d = np.loadtxt(f)
        zeta = np.unique(d[:, 0])
        x = d[: len(d) // len(zeta), 1]
        ez = d[:, 3].reshape(len(zeta), len(x)) * e0
        fig, ax = plt.subplots(figsize=(10, 4.2))
        vm = np.percentile(np.abs(ez), 99.5)
        im = ax.pcolormesh(zeta / 1e3 - zeta_of.get(s, 0.0), x / 1e3, ez.T, cmap="RdBu_r", vmin=-vm, vmax=vm,
                           shading="auto", rasterized=True)
        for r in (1.0, -1.0):
            ax.axhline(r, color="limegreen", lw=0.8)
        fig.colorbar(im, ax=ax, label=r"$E_z$  [MV/m]")
        ax.set_xlabel(r"$\xi$  [mm]")
        ax.set_ylabel(r"x  [mm]  (plane $\varphi$ = 0 / $\pi$)")
        ax.set_title("z = %.0f mm; green: inner radius of the quartz" % z, fontsize=10)
        fig.tight_layout()
        fig.savefig(os.path.join(a.run, "wake_map.png"), dpi=130)

    # driver
    fig, axs = plt.subplots(1, 3, figsize=(12, 3.4))
    axs[0].plot(b[:, 1] / 1e3, b[:, 2])
    axs[0].set_ylabel("mean energy [MeV]")
    axs[1].plot(b[:, 1] / 1e3, b[:, 3] * 100)
    axs[1].set_ylabel("rms energy spread [%]")
    axs[2].plot(b[:, 1] / 1e3, b[:, 4], label=r"$\sigma_x$ [µm]")
    axs[2].plot(b[:, 1] / 1e3, b[:, 5] * 10, label=r"$10\,\epsilon_{nx}$ [µm]")
    axs[2].legend(fontsize=8)
    for q in axs:
        q.set_xlabel("z [m]")
        q.grid(alpha=0.3)
        if z_diel_end is not None:
            q.axvline(z_diel_end / 1e3, color="tab:red", ls="--", lw=0.8)
    fig.tight_layout()
    fig.savefig(os.path.join(a.run, "beam.png"), dpi=130)
    print("wrote", ", ".join(os.path.join(a.run, n) for n in ("wake_axis.png", "wake_map.png", "beam.png")))


if __name__ == "__main__":
    main()
