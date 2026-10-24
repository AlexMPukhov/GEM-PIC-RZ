#!/bin/sh
# Lasers from and to LASY / openPMD envelope files (laser.file, output.lasy), files from tests/lasy_files.py:
#   1. Gaussian pulse at focus from a thetaMode file = the built-in Gaussian (axis E_x to 2e-3 of the peak);
#      the same on a cartesian (t, y, x) grid (5e-3); circular polarization (1, i)/sqrt(2) = built-in circular
#   2. the file at focus put at file_z = 0.3 mm: the field propagated back to z = 0 has the peak of the
#      built-in Gaussian with z_focus = 0.3 mm (0.5 %)
#   3. export at z = 0.3 mm (output.lasy) and re-import at its plane: the initial field equals the exported
#      level (2e-3 of the peak); with LASY installed, lasy.profiles.FromOpenPMDProfile reads the file
#   4. inject = boundary: the luminal Bessel pulse of the flying_focus block, given as a file and fed through the
#      absorber layer, reproduces the built-in flying focus after 0.3 mm (peak 1 %, rms of E_x 3e-3)
# usage: sh tests/lasy.sh <gempic_rz>     (exit 77: skipped, needs python3 with numpy and h5py)
set -e
abspath() { (cd "$(dirname "$1")" && echo "$(pwd)/$(basename "$1")"); }
EXE=$(abspath "$1"); shift
SRC=$(cd "$(dirname "$0")/.." && pwd)
python3 -c "import numpy, h5py" 2>/dev/null || { echo "lasy test skipped: needs python3 with numpy and h5py"; exit 77; }
WORK=$(mktemp -d)
cd "$WORK"
python3 "$SRC/tests/lasy_files.py"
cat > g.deck <<'EOF'
grid { lambda0 = 0.8 um; m_max = 1; radial_grid = uniform; dr = 0.4 um; r_max = 40 um
       dzeta = 0.0333 um; nzeta = 1800; dz = 10 um; z_end = 0.3 mm }
laser { a0 = 1; waist = 8 um; length_1e = 5.99585 um; zeta_center = 30 um }
output { every = 0; maps = false; particles = false }
numerics { iterations = 1 }
EOF
sed -e 's|^laser.*|laser { file = gauss_rt.h5; a0 = 1; zeta_center = 30 um }|' g.deck > f.deck
# compare the E_x columns of two axis files: max |difference| / max |E_x| of the first <= tol
cmpaxis() {  # <label> <file a> <file b> <tol>
  paste "$2" "$3" | awk -v l="$1" -v t="$4" '!/^#/ { d = $2 - $7; if (d < 0) d = -d; if (d > m) m = d;
    a = $2 < 0 ? -$2 : $2; if (a > p) p = a }
    END { r = m / p; printf "  %s: max |dE_x| / peak = %.2e (tolerance %g)", l, r, t;
          if (r <= t) print "  ok"; else { print "  FAILED"; exit 1 } }'
}
peak() { awk '!/^#/ { v = $2 < 0 ? -$2 : $2; if (v > m) m = v } END { printf "%.5f", m }' "$1"; }
check() {  # <label> <value> <expected> <relative tolerance>
  if awk -v v="$2" -v e="$3" -v t="$4" 'BEGIN { d = (v - e) / e; if (d < 0) d = -d; exit !(d <= t) }'; then
    echo "  $1 = $2 (expected $3, relative tolerance $4)  ok"
  else
    echo "  $1 = $2, expected $3 (relative tolerance $4)  FAILED"; exit 1
  fi
}
R="grid.z_end=1um"
# ---- 1. at the file's plane
"$EXE" g.deck $R output.dir=g0 > g0.log
"$EXE" f.deck $R output.dir=f0 > f0.log
cmpaxis "thetaMode file = built-in Gaussian" g0/axis_00000.txt f0/axis_00000.txt 2e-3
"$EXE" f.deck $R laser.file=gauss_xyt.h5 output.dir=x0 > x0.log
cmpaxis "cartesian file = built-in Gaussian" g0/axis_00000.txt x0/axis_00000.txt 5e-3
"$EXE" g.deck $R laser.polarization=circular output.dir=gc > gc.log
"$EXE" f.deck $R laser.file=gauss_circ.h5 laser.a0=1.41421356 output.dir=fc > fc.log
cmpaxis "circular polarization" gc/axis_00000.txt fc/axis_00000.txt 2e-3
# ---- 2. propagated from the file's plane (exact) vs the paraxial Gaussian before its focus
"$EXE" g.deck $R laser.z_focus=0.3mm output.dir=g3 > g3.log
"$EXE" f.deck $R laser.file_z=0.3mm output.dir=f3 > f3.log
grep -q "propagated from its plane" f3.log
check "peak E_x on the axis 0.3 mm before focus" "$(peak f3/axis_00000.txt)" "$(peak g3/axis_00000.txt)" 0.005
# ---- 3. export and re-import
"$EXE" f.deck laser.file_z=0.3mm output.lasy=true output.dir=fp > fp.log
X=fp/lasy/laser_00030.h5
test -f "$X"
"$EXE" f.deck $R laser.file=$X laser.a0=-1 output.dir=re > re.log 2>&1 && exit 1  # a0 must be positive
sed -e "s|^laser.*|laser { file = $X; zeta_center = 30 um; t_center = 0 fs }|" g.deck > re.deck
"$EXE" re.deck $R output.dir=re > re.log
cmpaxis "re-imported LASY export = exported level" fp/axis_00030.txt re/axis_00000.txt 2e-3
if python3 -c "import lasy" 2>/dev/null; then
  python3 -c "
from lasy.profiles.from_openpmd_profile import FromOpenPMDProfile
p = FromOpenPMDProfile('$X', 'laserEnvelope')
assert abs(p.lambda0 - 0.8e-6) < 1e-12, p.lambda0
print('  LASY reads the export: lambda0 = %.4g m, polarization %s  ok' % (p.lambda0, p.pol))"
else
  echo "  (LASY not installed: reading the export with LASY skipped)"
fi
# ---- 4. inject = boundary: Bessel pulse from a file against the built-in flying focus
V="$SRC/benchmarks/flying_focus/vacuum.deck"
sed -e 's|^flying_focus.*|laser { file = ffm.h5; zeta_center = 90 um; inject = boundary }|' "$V" > ffl.deck
O="grid.z_end=0.3mm grid.absorber_strength=1 output.every=0 output.maps=false"
"$EXE" "$V" $O output.dir=ffref "$@" > ffref.log
"$EXE" ffl.deck $O output.dir=fflasy "$@" > fflasy.log
check "flying focus from a file: peak E_x on the axis" "$(peak fflasy/axis_00030.txt)" "$(peak ffref/axis_00030.txt)" 0.01
paste ffref/axis_00030.txt fflasy/axis_00030.txt | awk '!/^#/ { d = $2 - $7; s += d * d; n++ }
  END { r = sqrt(s / n); printf "  flying focus from a file: rms of the E_x difference = %.2e (tolerance 3e-3)", r;
        if (r <= 3e-3) print "  ok"; else { print "  FAILED"; exit 1 } }'
cd - > /dev/null
rm -rf "$WORK"
echo "LASY test passed"
