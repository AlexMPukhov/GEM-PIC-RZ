"""GUI test (tools/gempic_gui.py), headless through its HTTP API:
  1. a run of tests/regression.deck (2 steps, frames every step) started from the GUI finishes, with
     three HDF5 frames and the bunch table;
  2. the fields that the GUI builds from data_NNNNN.h5 (modes, staggered faces/centres, axis regularity)
     equal the code's own gather: on the axis (axis_NNNNN.txt: E_x, E_z, n_e) and in the plane y = 0 at
     both half planes (map_NNNNN.txt: E_x, E_z), to the print precision of the text files;
  3. the GUI's text-map reader (builds without HDF5) gives the same values;
  4. map, transverse cut and bunch requests answer; a deck with an error is refused with the message;
  5. with mpirun: a 2-rank run gives the same frames.
usage: python3 tests/gui_check.py <gempic_rz>     (exit 77: skipped without numpy/h5py)"""
import json
import os
import shutil
import socket
import subprocess
import sys
import tempfile
import time
import urllib.error
import urllib.parse
import urllib.request

try:
    import numpy as np
    import h5py  # noqa: F401
except ImportError:
    print("SKIP: needs numpy and h5py")
    sys.exit(77)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "tools"))
import gempic_gui  # noqa: E402

exe = os.path.abspath(sys.argv[1])
fails = []


def check(cond, msg):
    print(("ok    " if cond else "FAIL  ") + msg)
    if not cond:
        fails.append(msg)


def free_port():
    s = socket.socket()
    s.bind(("127.0.0.1", 0))
    p = s.getsockname()[1]
    s.close()
    return p


work = tempfile.mkdtemp(prefix="gempic_gui_")
deck = os.path.join(work, "run.deck")
shutil.copy(os.path.join(ROOT, "tests", "regression.deck"), deck)
port = free_port()
srv = subprocess.Popen([sys.executable, os.path.join(ROOT, "tools", "gempic_gui.py"), deck, "--port", str(port),
                        "--no-browser", "--token", "T", "--exe", exe], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
B = f"http://127.0.0.1:{port}"


def get(path, **q):
    u = B + path + "?" + urllib.parse.urlencode(dict(q, token="T"))
    try:
        return json.load(urllib.request.urlopen(u, timeout=60))
    except urllib.error.HTTPError as e:
        return {"HTTP": e.code, **json.load(e)}


def post(path, d):
    r = urllib.request.Request(B + path + "?token=T", data=json.dumps(d).encode(),
                               headers={"Content-Type": "application/json"})
    try:
        return json.load(urllib.request.urlopen(r, timeout=120))
    except urllib.error.HTTPError as e:
        return {"HTTP": e.code, **json.load(e)}


def run(extra, ranks=1, outdir=None):
    text = open(deck).read()
    r = post("/api/run", {"path": deck, "text": text, "threads": 1, "ranks": ranks, "every": 1, "extra": extra})
    if "error" in r:
        return r
    for _ in range(300):
        s = get("/api/status")
        if not s["running"] and s["rc"] is not None:
            return s
        time.sleep(0.5)
    return {"error": "timeout"}


try:
    for _ in range(50):
        try:
            urllib.request.urlopen(B + "/api/init?token=T", timeout=2)
            break
        except OSError:
            time.sleep(0.2)
    check(get("/api/init").get("deck") == deck, "server answers, deck opened")
    check(urllib.request.urlopen(B + "/?token=T").status == 200, "page served")
    try:
        urllib.request.urlopen(B + "/api/init?token=wrong")
        check(False, "requests without the token are refused")
    except urllib.error.HTTPError as e:
        check(e.code == 403, "requests without the token are refused")

    s = run("grid.z_end=0.01mm\noutput.dir=out")
    check(s.get("rc") == 0 and s.get("step") == 2, f"run from the GUI finished (rc {s.get('rc')}, step {s.get('step')})")
    out = os.path.join(work, "out")
    f = get("/api/frames")
    check(f.get("frames") == [0, 1, 2] and f.get("kind") == "h5", f"frames {f.get('frames')} ({f.get('kind')})")
    check(f.get("comps") == ["ex", "ey", "ez", "bx", "by", "bz", "ne", "zni"], "field list")

    # 2. against the code's gather: axis and map files of the same run
    worst = {}
    for step in (1, 2):
        ax = np.loadtxt(os.path.join(out, f"axis_{step:05d}.txt"))
        line = get("/api/lines", step=step, comps="ex,ez,ne", x=0)
        for c, col in (("ex", 1), ("ez", 2), ("ne", 3)):
            v = np.array(line["series"][c])
            ref = ax[:, col]
            worst[f"axis {c}"] = max(worst.get(f"axis {c}", 0), np.abs(v - ref).max() / max(np.abs(ref).max(), 1e-30))
        mp = np.loadtxt(os.path.join(out, f"map_{step:05d}.txt"))
        zeta = np.unique(mp[:, 0])
        rec = gempic_gui.read_h5(os.path.join(out, f"data_{step:05d}.h5"))
        nx = len(mp) // len(zeta)
        m3 = mp.reshape(len(zeta), nx, mp.shape[1])
        xs = m3[0, :, 1]
        stride = int(round((zeta[1] - zeta[0]) / (rec["comps"]["ez"][1][1] - rec["comps"]["ez"][1][0])))
        for c, col in (("ex", 2), ("ez", 3)):
            _, v = gempic_gui.half_planes(rec, c, xs)
            v = v[::stride][: len(zeta)]
            ref = m3[:, :, col]
            worst[f"map {c}"] = max(worst.get(f"map {c}", 0), np.abs(v - ref).max() / max(np.abs(ref).max(), 1e-30))
        # 3. text-map reader
        tm = gempic_gui.read_map(os.path.join(out, f"map_{step:05d}.txt"))
        for c, col in (("ex", 2), ("ez", 3)):
            _, v = gempic_gui.half_planes(tm, c, xs)
            worst[f"text map {c}"] = max(worst.get(f"text map {c}", 0),
                                         np.abs(v - m3[:, :, col]).max() / max(np.abs(m3[:, :, col]).max(), 1e-30))
    for k, v in worst.items():
        check(v < 2e-5, f"GUI field vs code ({k}): max rel. difference {v:.1e}")

    # 4. other requests
    m = get("/api/map", step=2, comp="ez", nr=100, nz=300)
    check("data" in m and m["nz"] <= 300 and m["nx"] == 100 and m["signed"], "map request (thinned, signed colour map)")
    r = get("/api/rline", step=2, comps="ez,ne", zeta=20)
    check(len(r.get("x", [])) == 401 and abs(r["zeta"] - 20) < 0.1, "transverse cut")
    b = get("/api/beams")
    check("witness" in b and len(b["witness"]["z_mm"]) == 3, "bunch table (witness, 3 steps)")
    bad = post("/api/run", {"path": deck, "text": open(deck).read() + "\nnot_a_block { x = 1 }\n", "threads": 1,
                            "ranks": 1, "every": 1, "extra": ""})
    check(bad.get("HTTP") == 400 and "error" in bad, "a deck with an error is refused: " + bad.get("error", "")[:60].replace("\n", " "))
    with open(deck, "w") as fh:   # the GUI saved the broken deck: restore it
        fh.write(open(os.path.join(ROOT, "tests", "regression.deck")).read())
    c = post("/api/check", {"path": deck, "text": open(deck).read(), "extra": ""})
    check(c.get("ok") is True, "deck check")

    # 5. two MPI ranks
    ldd = subprocess.run(["ldd", exe], capture_output=True, text=True).stdout if shutil.which("ldd") else ""
    if shutil.which("mpirun") and "libmpi" in ldd:
        ref = gempic_gui.read_h5(os.path.join(out, "data_00002.h5"))
        s = run("grid.z_end=0.01mm\noutput.dir=out_p2", ranks=2)
        f2 = get("/api/frames")
        ok = s.get("rc") == 0 and f2.get("frames") == [0, 1, 2]
        d = 0.0
        if ok:
            p2 = gempic_gui.read_h5(os.path.join(work, "out_p2", "data_00002.h5"))
            for c in ("ex", "ez", "ne"):
                a0, a1 = ref["comps"][c][2], p2["comps"][c][2]
                d = max(d, np.abs(a0 - a1).max() / max(np.abs(a0).max(), 1e-30))
        check(ok and d < 1e-9, f"2 MPI ranks: frames and fields as serial (max rel. difference {d:.1e})")
finally:
    srv.terminate()
    srv.wait()
    shutil.rmtree(work, ignore_errors=True)

print("gui test: " + ("PASS" if not fails else f"FAIL ({len(fails)})"))
sys.exit(1 if fails else 0)
