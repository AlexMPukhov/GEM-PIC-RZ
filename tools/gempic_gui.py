#!/usr/bin/env python3
"""GEM-PIC-RZ GUI: input-deck editor, run control and live plots in the web browser.

    python3 tools/gempic_gui.py [run.deck] [--port 8766] [--exe PATH] [--no-browser]

Opens http://127.0.0.1:8766/?token=... in the default browser. Only the Python standard library and
numpy are needed, plus h5py for the HDF5 frames (without it, or for a build without HDF5, the GUI
reads the text maps map_NNNNN.txt). The page has no external dependencies, so the GUI also works
offline. On a remote machine, start it with --no-browser and forward the port:
    ssh -L 8766:127.0.0.1:8766 host     then open the printed URL locally.

Live plots: the GUI runs the deck with command-line overrides (the deck file is not changed):
output.every (frames every N steps, if given) and output.format with hdf5 added, so that every output
step writes data_NNNNN.h5 with all modes. It shows the frames as they appear: 2D maps in the plane
y = 0 (x = +r at phi = 0, x = -r at phi = pi), cuts along zeta at a given x and transverse cuts at a
given zeta, of E_x, E_y, E_z, B_x, B_y, B_z (code units m c omega0 / e, m omega0 / e) and of the
densities n_e/n0 and Z n_i/n0, and the bunch evolution from bunches.txt. Finished runs can be browsed
the same way. The deck is checked with `gempic_rz deck --check` before every run.

The server listens on 127.0.0.1 only and every request needs the random token of the URL.
"""
import argparse
import base64
import json
import os
import re
import secrets
import shutil
import signal
import subprocess
import sys
import threading
import time
import webbrowser
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from urllib.parse import parse_qs, urlparse

import numpy as np

try:
    import h5py
except ImportError:   # text maps only
    h5py = None

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

STEP_RE = re.compile(r"^step\s+(\d+)\s+z\s*=\s*([-+0-9.eE]+)\s*mm")
NSTEPS_RE = re.compile(r"\((\d+) steps\)")
OUTPUT_RE = re.compile(r"^output: dir '([^']*)', every (\d+) steps")
FORMAT_RE = re.compile(r"^output format: (.*?)(?: \(|$)")
H5_RE = re.compile(r"^data_(\d+)\.h5$")
MAP_RE = re.compile(r"^map_(\d+)\.txt$")
SIGNED = {"ex", "ey", "ez", "bx", "by", "bz"}          # diverging colour map; densities sequential
FIELDS = ["ex", "ey", "ez", "bx", "by", "bz", "ne", "zni"]


# ----------------------------------------------------------------------------- deck check
def check_deck(exe, deck_path, overrides=(), cwd=None):
    """run `gempic_rz deck --check`: output dir (relative to cwd), output.every, number of steps,
    output formats and the messages (errors)"""
    p = subprocess.run([exe, deck_path, "--check"] + list(overrides), cwd=cwd or os.path.dirname(deck_path),
                       capture_output=True, text=True, timeout=120)
    info = {"ok": p.returncode == 0, "outdir": None, "every": None, "steps": 0, "formats": ["text"],
            "messages": (p.stdout + p.stderr).splitlines()}
    for line in info["messages"]:
        m = OUTPUT_RE.match(line)
        if m:
            info["outdir"], info["every"] = m.group(1), int(m.group(2))
        m = NSTEPS_RE.search(line)
        if m and "z_end" in line:
            info["steps"] = int(m.group(1))
        m = FORMAT_RE.match(line)
        if m:
            info["formats"] = [f.strip() for f in m.group(1).split("+")]
    return info


# ----------------------------------------------------------------------------- frames
def _plane(F, phi):
    """F(r, phi) = F_0 + 2 Re sum_m F_m e^{i m phi} for modes F[m, zeta, r]"""
    out = np.real(F[0]).copy()
    for m in range(1, F.shape[0]):
        out += 2.0 * np.real(F[m] * np.exp(1j * m * phi))
    return out


def read_h5(fn):
    """fields in the plane y = 0 from data_NNNNN.h5: for each field the two half planes
    (phi = 0: x = +r, phi = pi: x = -r) on its own radial nodes, and the zeta of its slices"""
    with h5py.File(fn, "r") as f:
        z_mm = float(f.attrs["z_m"]) * 1e3
        um = float(f.attrs["k0inv_m"]) * 1e6
        rf = f["grid/r_faces"][()] * um
        rc = f["grid/r_centres"][()] * um
        zs = f["grid/zeta_slices"][()][:-1] * um     # the last slice is the empty ghost ahead of the window
        zh = f["grid/zeta_halfcells"][()] * um
        raw = {rec: f["fields/" + rec][()][:, :-1, :] for rec in ("Er", "Ephi", "Ez", "Br", "Bphi", "Bz")}
        rc0 = np.concatenate([[0.0], rc])   # centre fields get an axis node (regularity, as the code's gather)
        c0, c1 = rc[0], rc[1]
        ext = {}   # m = 1 values of E_phi, B_r on the axis: linear extrapolation of the centres
        for rec in ("Ephi", "Br"):
            F = raw[rec]
            ext[rec] = (c1 * F[1, :, 0] - c0 * F[1, :, 1]) / (c1 - c0) if F.shape[0] > 1 else None
        for rec in ("Ephi", "Ez", "Br"):
            F = raw[rec]
            ax = np.zeros(F.shape[:2], dtype=F.dtype)
            if rec == "Ez":
                ax[0] = F[0, :, 0]                      # m = 0 finite on the axis, m >= 1 vanish
            elif F.shape[0] > 1:
                ax[1] = ext[rec]                        # m = 1 finite, m = 0 and m >= 2 vanish
            raw[rec] = np.concatenate([ax[:, :, None], F], axis=2)
        # faces: face 0 is not used by the solver; its value follows from regularity
        # (E_r = -i E_phi, B_phi = i B_r on the axis for m = 1; 0 for the other modes; B_z as stored)
        for rec, src, fac in (("Er", "Ephi", -1j), ("Bphi", "Br", 1j)):
            F = raw[rec].copy()
            F[:, :, 0] = 0.0
            if F.shape[0] > 1:
                F[1, :, 0] = fac * ext[src]
            raw[rec] = F
        comps = {}
        # vector components in the plane: x = r cos(phi), so E_x = E_r and E_y = E_phi at phi = 0,
        # E_x = -E_r and E_y = -E_phi at phi = pi
        for name, rec, nodes in [("ex", "Er", rf), ("ey", "Ephi", rc0), ("ez", "Ez", rc0),
                                 ("bx", "Br", rc0), ("by", "Bphi", rf), ("bz", "Bz", rf)]:
            F = raw[rec]
            p0, p1 = _plane(F, 0.0), _plane(F, np.pi)
            if name in ("ex", "ey", "bx", "by"):
                p1 = -p1
            comps[name] = (nodes, zs, p0, p1)
        for name, rec in [("ne", "ne_over_n0"), ("zni", "Zni_over_n0")]:
            if "densities/" + rec in f:
                d = f["densities/" + rec][()]          # m = 0: [half-cell][centre]
                comps[name] = (rc, zh, d, d)
    return {"z_mm": z_mm, "comps": comps}


def read_map(fn):
    """text map map_NNNNN.txt (builds without HDF5): zeta x Ex Ez ne Zni in the plane phi = 0 / pi"""
    a = np.loadtxt(fn, ndmin=2)
    head = open(fn).readline()
    m = re.search(r"z = ([-+0-9.eE]+) mm", head)
    zeta = np.unique(a[:, 0])
    nx = len(a) // len(zeta)
    a = a[: nx * len(zeta)].reshape(len(zeta), nx, a.shape[1])
    x = a[0, :, 1]
    pos, neg = x >= 0, x < 0
    rp, rn = x[pos], -x[neg][::-1]
    comps = {}
    for name, col in [("ex", 2), ("ez", 3), ("ne", 4), ("zni", 5)]:
        p0 = a[:, pos, col]
        p1 = a[:, neg, col][:, ::-1]
        # store both halves on the nodes of the positive half (the map is symmetric in |x|)
        if len(rn) == len(rp) and np.allclose(rn, rp):
            comps[name] = (rp, zeta, p0, p1)
        else:
            comps[name] = (rp, zeta, p0, np.array([np.interp(rp, rn, row) for row in p1]))
    return {"z_mm": float(m.group(1)) if m else 0.0, "comps": comps}


class FrameCache:
    def __init__(self, n=4):
        self.n, self.d, self.lock = n, {}, threading.Lock()

    def get(self, fn):
        mt = os.path.getmtime(fn)
        with self.lock:
            e = self.d.get(fn)
            if e and e[0] == mt:
                return e[1]
        try:
            rec = read_h5(fn) if fn.endswith(".h5") else read_map(fn)
        except Exception:   # being written (text maps) or unreadable
            return None
        with self.lock:
            self.d[fn] = (mt, rec)
            while len(self.d) > self.n:
                self.d.pop(next(iter(self.d)))
        return rec


CACHE = FrameCache()


def half_planes(rec, name, xs):
    """field `name` on the transverse positions xs (x >= 0: phi = 0, x < 0: phi = pi), [zeta][x]"""
    nodes, zeta, p0, p1 = rec["comps"][name]
    xs = np.asarray(xs, float)
    out = np.empty((len(zeta), len(xs)))
    for half, sel in ((p0, xs >= 0), (p1, xs < 0)):
        r = np.abs(xs[sel])
        j = np.clip(np.searchsorted(nodes, r) - 1, 0, len(nodes) - 2)
        w = np.clip((r - nodes[j]) / (nodes[j + 1] - nodes[j]), 0, 1)
        out[:, sel] = half[:, j] * (1 - w) + half[:, j + 1] * w
    return zeta, out


# ----------------------------------------------------------------------------- run control
class Runner:
    def __init__(self, exe):
        self.exe = exe
        self.proc = None
        self.lock = threading.Lock()
        self.log = []
        self.step = -1
        self.z = 0.0
        self.rc = None
        self.workdir = os.getcwd()
        self.outdir = None
        self.steps = 0
        self.started = None
        self.ended = None

    def running(self):
        return self.proc is not None and self.proc.poll() is None

    def overrides(self, every, extra, formats):
        ov = [e for e in extra if "=" in e]
        if every > 0:
            ov.append(f"output.every={every}")
        if h5py is not None and "hdf5" not in formats:
            ov.append("output.format=" + ", ".join(formats + ["hdf5"]))
        return ov

    def start(self, deck_path, text, threads, ranks, every, extra):
        if self.running():
            raise RuntimeError("a run is already in progress")
        deck_path = os.path.abspath(deck_path)
        with open(deck_path, "w") as f:
            f.write(text)
        wd = os.path.dirname(deck_path)
        base = [e for e in extra if "=" in e]
        info = check_deck(self.exe, deck_path, base, wd)
        if not info["ok"]:
            raise RuntimeError("the deck has errors:\n" + "\n".join(m for m in info["messages"] if m.strip())[-3000:])
        ov = self.overrides(every, extra, info["formats"])
        if threads > 1 and not any(o.startswith("numerics.backend") for o in ov):
            ov.append("numerics.backend=openmp")
        info = check_deck(self.exe, deck_path, ov, wd)   # with the overrides (e.g. a build without HDF5)
        if not info["ok"]:
            ov = [o for o in ov if not o.startswith("output.format")]
            info = check_deck(self.exe, deck_path, ov, wd)
            if not info["ok"]:
                raise RuntimeError("\n".join(info["messages"])[-3000:])
        cmd = [self.exe, deck_path] + ov
        if ranks > 1:
            mpi = shutil.which("mpirun") or shutil.which("mpiexec")
            if mpi is None:
                raise RuntimeError("mpirun not found")
            cmd = [mpi, "-np", str(ranks)] + (["--allow-run-as-root"] if os.geteuid() == 0 else []) + cmd
        env = dict(os.environ)
        env["OMP_NUM_THREADS"] = str(max(1, threads))
        env.setdefault("OMP_PROC_BIND", "false")
        outdir = os.path.join(wd, info["outdir"] or "diag")
        # remove old frames of the same output directory (a new run starts at step 0)
        if os.path.isdir(outdir):
            for fn in os.listdir(outdir):
                if H5_RE.match(fn) or MAP_RE.match(fn) or fn.endswith(".h5.tmp"):
                    os.remove(os.path.join(outdir, fn))
        with self.lock:
            self.log = ["$ " + " ".join(cmd), ""]
            self.step, self.z, self.rc = -1, 0.0, None
            self.workdir, self.outdir, self.steps = wd, outdir, info["steps"]
            self.started, self.ended = time.time(), None
        self.proc = subprocess.Popen(cmd, cwd=wd, env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                     text=True, bufsize=1, start_new_session=True)
        threading.Thread(target=self._reader, args=(self.proc,), daemon=True).start()

    def _reader(self, proc):
        for line in proc.stdout:
            line = line.rstrip("\n")
            m = STEP_RE.match(line)
            with self.lock:
                self.log.append(line)
                if len(self.log) > 4000:
                    del self.log[:1000]
                if m:
                    self.step, self.z = int(m.group(1)), float(m.group(2))
        proc.wait()
        with self.lock:
            self.rc = proc.returncode
            self.ended = time.time()
            self.log.append(f"[exit code {proc.returncode}]")

    def fresh(self, fn):
        """the file exists and, during a run, was written by this run"""
        if not os.path.exists(fn):
            return False
        return not (self.running() and self.started and os.path.getmtime(fn) < self.started)

    def stop(self):
        if self.running():
            try:
                os.killpg(self.proc.pid, signal.SIGTERM)
            except ProcessLookupError:
                pass

    def frames(self, outdir=None):
        """(steps, kind): HDF5 data files if there are any, else text maps"""
        d = outdir or self.outdir
        if d is None or not os.path.isdir(d):
            return [], None
        files = os.listdir(d)
        h5 = sorted(int(m.group(1)) for m in (H5_RE.match(f) for f in files) if m)
        if h5 and h5py is not None:
            return h5, "h5"
        maps = sorted(int(m.group(1)) for m in (MAP_RE.match(f) for f in files) if m)
        if self.running() and maps:
            maps = maps[:-1]   # the newest text map may still be written
        return maps, "map"

    def frame_file(self, step, outdir=None):
        d = outdir or self.outdir
        steps, kind = self.frames(d)
        if kind == "h5":
            return os.path.join(d, f"data_{step:05d}.h5")
        return os.path.join(d, f"map_{step:05d}.txt")

    def status(self, nlog=300):
        with self.lock:
            return {"running": self.running(), "rc": self.rc, "step": self.step, "z": self.z, "steps": self.steps,
                    "outdir": self.outdir, "workdir": self.workdir, "log": self.log[-nlog:],
                    "elapsed": ((self.ended or time.time()) - self.started) if self.started else 0.0}


def read_bunches(fn):
    """bunches.txt -> {name: {column: list}} (one line per step and bunch)"""
    cols = ["step", "z_mm", "name", "N", "Q_pC", "E_MeV", "dE_E", "sig_x_um", "sig_y_um", "eps_nx_um",
            "eps_ny_um", "zeta_um", "sig_z_um"]
    out = {}
    with open(fn) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            p = line.split()
            if len(p) < len(cols):
                continue
            d = out.setdefault(p[2], {c: [] for c in cols if c != "name"})
            for c, v in zip(cols, p):
                if c != "name":
                    d[c].append(float(v))
    return out


# ----------------------------------------------------------------------------- HTTP
class Handler(BaseHTTPRequestHandler):
    server_version = "GEM-PIC-RZ-GUI"

    def log_message(self, *a):
        pass

    def _send(self, code, body, ctype="application/json"):
        if isinstance(body, (dict, list)):
            body = json.dumps(body).encode()
        elif isinstance(body, str):
            body = body.encode()
        self.send_response(code)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(body)

    def _auth(self, q):
        tok = self.headers.get("X-Token") or q.get("token", [""])[0]
        return secrets.compare_digest(tok, self.server.token)

    def do_GET(self):
        u = urlparse(self.path)
        q = parse_qs(u.query)
        if u.path == "/":
            if not self._auth(q):
                return self._send(403, "forbidden: open the URL with ?token=... printed by gempic_gui.py", "text/plain")
            with open(os.path.join(HERE, "gempic_gui.html"), "rb") as f:
                return self._send(200, f.read(), "text/html; charset=utf-8")
        if not self._auth(q):
            return self._send(403, {"error": "forbidden"})
        try:
            return self._send(200, self.api_get(u.path, q))
        except Exception as e:  # report to the page
            return self._send(400, {"error": str(e)})

    def do_POST(self):
        u = urlparse(self.path)
        if not self._auth(parse_qs(u.query)):
            return self._send(403, {"error": "forbidden"})
        n = int(self.headers.get("Content-Length", "0"))
        data = json.loads(self.rfile.read(n) or b"{}")
        try:
            return self._send(200, self.api_post(u.path, data))
        except Exception as e:
            return self._send(400, {"error": str(e)})

    # ------------------------------------------------------------------ GET API
    def api_get(self, path, q):
        R = self.server.runner
        a = lambda k, d=None: q.get(k, [d])[0]   # noqa: E731
        if path == "/api/init":
            p = self.server.deck
            text = open(p).read() if p and os.path.exists(p) else ""
            return {"deck": p or os.path.join(os.getcwd(), "run.deck"), "text": text, "examples": self.examples(),
                    "exe": R.exe, "cpus": os.cpu_count() or 1, "mpirun": bool(shutil.which("mpirun")),
                    "h5py": h5py is not None}
        if path == "/api/deck":
            p = os.path.abspath(a("path"))
            return {"path": p, "text": open(p).read()}
        if path == "/api/status":
            return R.status()
        od = a("outdir")
        if path == "/api/frames":
            steps, kind = R.frames(od)
            comps, info = [], {}
            if steps:
                rec = CACHE.get(R.frame_file(steps[-1], od))
                if rec is not None:
                    comps = [c for c in FIELDS if c in rec["comps"]]
                    nodes, zeta, _, _ = rec["comps"][comps[0]]
                    info = {"zeta": [float(zeta[0]), float(zeta[-1])], "rmax": float(nodes[-1])}
            base = od or R.outdir
            beams = bool(base) and R.fresh(os.path.join(base, "bunches.txt"))
            return {"dir": base, "frames": steps, "kind": kind, "comps": comps, **info, "beams": beams}
        if path in ("/api/map", "/api/lines", "/api/rline"):
            rec = CACHE.get(R.frame_file(int(a("step")), od))
            if rec is None:
                raise RuntimeError("frame not available")
            comp = a("comp", "ez")
            if path == "/api/map":
                nodes = rec["comps"][comp][0]
                rmax = min(float(a("rmax", 0) or nodes[-1]), nodes[-1])
                nx = int(a("nr", 240))
                upper = a("upper", "0") == "1"
                xs = np.linspace(0 if upper else -rmax, rmax, nx)
                zeta, g = half_planes(rec, comp, xs)
                nz = int(a("nz", 900))
                if len(zeta) > nz:   # thin out in zeta by block maxima of |f| (keeps narrow features)
                    k = int(np.ceil(len(zeta) / nz))
                    n = len(zeta) // k
                    b = g[: n * k].reshape(n, k, nx)
                    idx = np.abs(b).argmax(axis=1)
                    g = np.take_along_axis(b, idx[:, None, :], axis=1)[:, 0, :]
                    zs = zeta[: n * k].reshape(n, k).mean(axis=1)
                else:
                    zs = zeta
                fin = g[np.isfinite(g)]
                vmin, vmax = (float(fin.min()), float(fin.max())) if fin.size else (0.0, 1.0)
                # robust range: 99.5 % of the values (a closure spike would otherwise hide everything)
                if fin.size and comp in SIGNED:
                    a99 = float(np.percentile(np.abs(fin), 99.5))
                    rlo, rhi = -a99, a99
                elif fin.size:
                    rlo, rhi = float(np.percentile(fin, 0.5)), float(np.percentile(fin, 99.5))
                else:
                    rlo, rhi = vmin, vmax
                return {"z_mm": rec["z_mm"], "comp": comp, "signed": comp in SIGNED, "zeta0": float(zs[0]),
                        "zeta1": float(zs[-1]), "nz": int(len(zs)), "x0": float(xs[0]), "x1": float(xs[-1]), "nx": nx,
                        "vmin": vmin, "vmax": vmax, "rlo": rlo, "rhi": rhi,
                        "data": base64.b64encode(np.ascontiguousarray(g, dtype=np.float32).tobytes()).decode()}
            comps = [c for c in a("comps", "ez").split(",") if c in rec["comps"]]
            if path == "/api/lines":   # along zeta at the transverse position x
                x0 = float(a("x", 0))
                out, zeta = {}, None
                for c in comps:
                    zeta, v = half_planes(rec, c, [x0])
                    out[c] = v[:, 0].tolist()
                return {"z_mm": rec["z_mm"], "x": zeta.tolist() if zeta is not None else [], "series": out, "xpos": x0}
            zc = float(a("zeta", 0))   # transverse cut at zeta
            nodes = rec["comps"][comps[0]][0] if comps else np.array([1.0])
            rmax = min(float(a("rmax", 0) or nodes[-1]), nodes[-1])
            xs = np.linspace(-rmax, rmax, 401)
            out, zk = {}, zc
            for c in comps:
                zeta, v = half_planes(rec, c, xs)
                k = int(np.argmin(np.abs(zeta - zc)))
                zk = float(zeta[k])
                out[c] = v[k].tolist()
            return {"z_mm": rec["z_mm"], "x": xs.tolist(), "zeta": zk, "series": out}
        if path == "/api/beams":
            base = od or R.outdir
            fn = os.path.join(base or ".", "bunches.txt")
            if not R.fresh(fn):
                return {}
            try:
                return read_bunches(fn)
            except Exception:
                return {}
        raise RuntimeError("unknown request " + path)

    def examples(self):
        out = []
        for sub in ("examples", "tests"):
            d = os.path.join(ROOT, sub)
            if os.path.isdir(d):
                out += sorted(os.path.join(d, f) for f in os.listdir(d) if f.endswith(".deck"))
        for sub in sorted(os.listdir(os.path.join(ROOT, "benchmarks"))) if os.path.isdir(os.path.join(ROOT, "benchmarks")) else []:
            p = os.path.join(ROOT, "benchmarks", sub, "gem.deck")
            if os.path.exists(p):
                out.append(p)
        return out

    # ------------------------------------------------------------------ POST API
    def api_post(self, path, d):
        R = self.server.runner
        if path == "/api/save":
            p = os.path.abspath(d["path"])
            with open(p, "w") as f:
                f.write(d["text"])
            return {"path": p}
        if path == "/api/check":
            p = os.path.abspath(d["path"])
            with open(p, "w") as f:
                f.write(d["text"])
            extra = [s.strip() for s in d.get("extra", "").split("\n") if "=" in s]
            info = check_deck(R.exe, p, extra)
            return {"ok": info["ok"], "messages": info["messages"]}
        if path == "/api/run":
            R.start(d["path"], d["text"], int(d.get("threads", 1)), int(d.get("ranks", 1)), int(d.get("every", 0) or 0),
                    [s.strip() for s in d.get("extra", "").split("\n") if s.strip()])
            return {"ok": True}
        if path == "/api/stop":
            R.stop()
            return {"ok": True}
        if path == "/api/open":   # browse the output of a finished run
            p = os.path.abspath(d["path"])
            if d.get("text"):
                with open(p, "w") as f:
                    f.write(d["text"])
            info = check_deck(R.exe, p)
            if info["outdir"] is None:
                raise RuntimeError("cannot find the output directory of this deck:\n" + "\n".join(info["messages"])[-2000:])
            R.outdir, R.workdir = os.path.join(os.path.dirname(p), info["outdir"]), os.path.dirname(p)
            return {"outdir": R.outdir}
        raise RuntimeError("unknown request " + path)


def find_exe(arg):
    for c in [arg, os.path.join(ROOT, "build", "gempic_rz"), shutil.which("gempic_rz")]:
        if c and os.path.isfile(c) and os.access(c, os.X_OK):
            return os.path.abspath(c)
    return arg or "gempic_rz"


def main():
    ap = argparse.ArgumentParser(description="GEM-PIC-RZ GUI (deck editor, run control, live plots in the browser)")
    ap.add_argument("deck", nargs="?", help="input deck to open")
    ap.add_argument("--port", type=int, default=8766)
    ap.add_argument("--exe", help="gempic_rz executable (default: build/gempic_rz of this repository, then PATH)")
    ap.add_argument("--no-browser", action="store_true", help="do not open a browser (e.g. on a remote machine)")
    ap.add_argument("--token", help=argparse.SUPPRESS)
    args = ap.parse_args()
    srv = ThreadingHTTPServer(("127.0.0.1", args.port), Handler)
    srv.daemon_threads = True
    srv.token = args.token or secrets.token_urlsafe(16)
    srv.runner = Runner(find_exe(args.exe))
    srv.deck = os.path.abspath(args.deck) if args.deck else None
    if srv.deck and os.path.exists(srv.deck):
        try:
            info = check_deck(srv.runner.exe, srv.deck)
            if info["outdir"]:
                srv.runner.outdir = os.path.join(os.path.dirname(srv.deck), info["outdir"])
                srv.runner.workdir = os.path.dirname(srv.deck)
        except (OSError, subprocess.SubprocessError):
            pass
    url = f"http://127.0.0.1:{args.port}/?token={srv.token}"
    print(f"GEM-PIC-RZ GUI: {url}\n  executable: {srv.runner.exe}\n  Ctrl-C to quit", flush=True)
    if not args.no_browser:
        threading.Timer(0.5, lambda: webbrowser.open(url)).start()
    try:
        srv.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        srv.runner.stop()


if __name__ == "__main__":
    main()
