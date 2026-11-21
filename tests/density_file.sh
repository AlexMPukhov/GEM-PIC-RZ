#!/bin/sh
# Density profiles from files (plasma.density_file, ion.density_file) against the same profiles given in the deck
# (tests/regression.deck: laser, plasma up-ramp 0 -> 1 over 400 k0^-1, H+ ions, witness bunch; 4 steps):
#   1. n(z) text table, a factor on `density`, with the points of the deck's z_points / z_values:
#      identical to the inline table (serial, openmp; cuda: to rounding)
#   2. the same ramp as absolute density in cm^-3 over um (other units: the numbers differ in the last digit), as a
#      CSV file with a header line, and shifted with file_z: equal to rounding
#   3. n(r, z) text table (ramp in z x a hollow radial profile, linear between the grid points) against the same
#      profile given as density_function(z, r); restart from a checkpoint = uninterrupted run
#   4. openPMD meshes (tests/density_files.py; needs python3 with numpy and h5py, else skipped): the profile of 3. on
#      a regular thetaMode mesh, on a rectilinear mesh (gridCoordinates) as absolute density in m^-3, and a smooth
#      gas jet sampled on 81 x 41 points against its density_function (axis fields to 1e-3 of the peak)
#   5. an ion species with its own profile from a file (pre-ionized H, free_electrons = false) = the inline table
#   6. start-up checks: summary (ranges, peak), warnings (file cut at a finite density, gradient shorter than
#      4 steps), rejected files (negative, nan, z not increasing, incomplete grid, no header, wrong unit) and
#      rejected combinations (z_points or radial_profile with the file, file outside the run)
# usage: sh tests/density_file.sh <gempic_rz> <compare_columns> [backend]
set -e
abspath() { (cd "$(dirname "$1")" && echo "$(pwd)/$(basename "$1")"); }
EXE=$(abspath "$1"); CMP=$(abspath "$2"); BACKEND=${3:-serial}
SRC=$(cd "$(dirname "$0")/.." && pwd)
WORK=$(mktemp -d)
cd "$WORK"
R="$SRC/tests/regression.deck"
O="numerics.backend=$BACKEND output.every=1 output.map_stride=8 grid.z_end=160k0^-1"
# serial / openmp: identical; cuda (atomic deposition in a nondeterministic order): equal to rounding
if [ "$BACKEND" = cuda ]; then T0="1e-6 1e-9"; else T0="0 0"; fi
T1="1e-6 1e-9"
FILES="axis_00002.txt axis_00004.txt map_00004.txt bunches.txt"
same() { for f in $FILES; do "$CMP" $1/$f $2/$f $3 $4; done; }   # <dir a> <dir b> <tolerances>
# replace the inline table of the plasma block by a file
withfile() { sed -e "s|^  z_points .*|  density_file  = $1|" -e '/^  z_values /d' "$R"; }

# ---- 1. n(z), factor on density, the points of the deck
cat > ramp.txt <<'EOF'
# up-ramp of tests/regression.deck
# z[k0^-1]  n
0    0
400  1
EOF
withfile ramp.txt > f1.deck
grep -q density_file f1.deck && ! grep -q z_values f1.deck
"$EXE" "$R" $O output.dir=ref > ref.log
"$EXE" f1.deck $O output.dir=f1 > f1.log
same ref f1 $T0
grep -q "density file ramp.txt (text; # z\[k0^-1\]  n): n(z) on 2 points, factor on 'density'" f1.log
echo "  n(z) text table (factor) = inline z_points / z_values: tolerances $T0  ok"

# ---- 2. absolute density, other units, CSV, file_z
# 0.002 n_c at 0.8 um = 0.002 * 1.1148e21 / 0.64 cm^-3; 400 k0^-1 = 400 * 0.8 / (2 pi) um
awk 'BEGIN { n = 0.002 * 1.1148e21 / 0.64; z = 400 * 0.8 / (2 * 3.14159265358979323846)
  printf "z[um],n_e[cm^-3]\n0,0\n%.17g,%.17g\n%.17g,%.17g\n", z, n, 3 * z, n > "ramp_abs.csv"
  printf "# z [mm] ; n [m^-3]\n%.17g ; 0\n%.17g ; %.17g\n%.17g ; %.17g\n", 1, 1 + z * 1e-3, n * 1e6, 1 + 3e-3 * z, n * 1e6 > "ramp_shift.txt" }'
withfile ramp_abs.csv > f2.deck
"$EXE" f2.deck $O output.dir=f2 > f2.log
same ref f2 $T1
grep -q "absolute density \[cm^-3\]" f2.log
grep -q "peak 3.484e+18 cm^-3 (0.002 n_c)" f2.log
withfile ramp_shift.txt > f2s.deck
"$EXE" f2s.deck $O plasma.file_z=-1mm output.dir=f2s > f2s.log
same ref f2s $T1
echo "  n(z) in cm^-3 over um (CSV) and in m^-3 over mm shifted with file_z = inline table: tolerances $T1  ok"

# ---- 3. n(r, z) against density_function
# f(z) = min(z / 400, 1), g(r) = 0.6 + 0.8 min(r / 75, 1): linear between the grid points of the file
awk 'BEGIN { print "# z[k0^-1] r[k0^-1] n"
  nz = split("0 400 1000", z, " "); nr = split("0 75 150 400", r, " ")
  for (i = 1; i <= nz; i++) for (j = 1; j <= nr; j++) {
    f = z[i] < 400 ? z[i] / 400 : 1; g = 0.6 + 0.8 * (r[j] < 75 ? r[j] / 75 : 1)
    printf "%g %g %.17g\n", z[i], r[j], f * g } }' > jet2d.txt
withfile jet2d.txt > f3.deck
sed -e 's|^  z_points .*|  density_function = "min(z / (400 k0^-1), 1) * (0.6 + 0.8 * min(r / (75 k0^-1), 1))"|' \
    -e '/^  z_values /d' "$R" > fn3.deck
"$EXE" fn3.deck $O output.dir=fn3 > fn3.log
"$EXE" f3.deck $O output.dir=f3 > f3.log
same fn3 f3 $T1
grep -q "n(r, z) on 4 x 3 points" f3.log
grep -q "radial profile: from the density file" f3.log
"$EXE" f3.deck $O output.dir=f3r grid.z_end=80k0^-1 output.checkpoint_every=2 > f3r1.log
"$EXE" f3.deck $O output.dir=f3r --restart latest > f3r2.log
for f in axis_00004.txt map_00004.txt bunches.txt; do "$CMP" f3/$f f3r/$f $T0; done
echo "  n(r, z) text table = density_function(z, r): tolerances $T1; restart = uninterrupted run ($T0)  ok"

# ---- 4. openPMD meshes
if python3 -c "import numpy, h5py" 2>/dev/null; then
  python3 "$SRC/tests/density_files.py"
  withfile jet2d_theta.h5 > f4.deck
  "$EXE" f4.deck $O output.dir=f4a > f4a.log
  same fn3 f4a $T1
  grep -q "openPMD; mesh 'density', iteration 0" f4a.log
  "$EXE" f4.deck $O plasma.density_file=jet2d_rect.h5 plasma.file_mesh=n_gas output.dir=f4b > f4b.log
  same fn3 f4b $T1
  grep -q "absolute density \[m^-3\]" f4b.log
  # smooth gas jet: linear interpolation of 81 x 41 samples against the function
  JET='exp(-((z - 150 k0^-1) / (200 k0^-1))^2 - (r / (120 k0^-1))^2)'
  sed -e "s|^  z_points .*|  density_function = \"$JET\"|" -e '/^  z_values /d' "$R" > fnj.deck
  "$EXE" fnj.deck $O output.dir=fnj > fnj.log
  "$EXE" f4.deck $O plasma.density_file=gauss_jet.h5 output.dir=f4j > f4j.log
  paste fnj/axis_00004.txt f4j/axis_00004.txt | awk '!/^#/ { h = NF / 2
      for (c = 2; c <= h; c++) { d = $c - $(c + h); if (d < 0) d = -d; if (d > m[c]) m[c] = d
        a = $c < 0 ? -$c : $c; if (a > p[c]) p[c] = a } }
    END { bad = 0
      for (c = 2; c <= 4; c++) { r = m[c] / p[c]; if (r > w) w = r }
      printf "  smooth jet, 81 x 41 samples: largest deviation of E_x, E_z, n_e on the axis = %.2e of the peak (tolerance 1e-3)", w
      if (w <= 1e-3) print "  ok"; else { print "  FAILED"; exit 1 } }'
  echo "  openPMD: regular thetaMode mesh and rectilinear mesh in m^-3 = density_function: tolerances $T1  ok"
else
  echo "  openPMD part skipped: needs python3 with numpy and h5py"
fi

# ---- 5. ion block with its own profile from a file
sed -e 's/^  density       = 0.002 n_c/  density       = 0.002 n_c\n  free_electrons = false/' \
    -e 's/^  element = H/  element = H\n  profile = own\n  density = 0.002 n_c\n  XXPROF\n  r_max = 150 k0^-1/' "$R" > pre.deck
sed 's/^  XXPROF/  z_points = 0, 400 k0^-1\n  z_values = 0, 1/' pre.deck > i0.deck
sed 's/^  XXPROF/  density_file = ramp.txt/' pre.deck > i1.deck
"$EXE" i0.deck $O output.dir=i0 > i0.log
"$EXE" i1.deck $O output.dir=i1 > i1.log
same i0 i1 $T0
# absolute atom density: the block's density is then not given
sed -e 's/^  XXPROF/  density_file = ramp_abs.csv/' -e '/^  density = 0.002 n_c/d' pre.deck > i2.deck
"$EXE" i2.deck $O output.dir=i2 > i2.log
same i0 i2 $T1
echo "  ion block, profile = own, from the files = inline table: tolerances $T0 (factor), $T1 (cm^-3)  ok"

# ---- 6. start-up checks
expect_fail() {  # <label> <text expected in the message> <deck> [overrides]
  l=$1; t=$2; shift 2
  if "$EXE" "$@" --check > chk.log 2>&1; then echo "  $l: accepted, expected an error  FAILED"; exit 1; fi
  grep -q -- "$t" chk.log || { echo "  $l: error message without '$t':"; cat chk.log; exit 1; }
}
expect_warn() {  # <label> <text> <deck> [overrides]
  l=$1; t=$2; shift 2
  "$EXE" "$@" --check > chk.log 2>&1 || { echo "  $l: rejected:"; cat chk.log; exit 1; }
  grep "WARNING" chk.log | grep -q -- "$t" || { echo "  $l: no warning with '$t':"; cat chk.log; exit 1; }
}
"$EXE" f1.deck --check > chk.log 2>&1
grep -q "z = 0 ... 0.0509296 mm (run: 0 ... 0.0254648 mm)" chk.log
grep -q "peak 3.484e+18 cm^-3 (0.002 n_c) at z = 0.0509296 mm" chk.log
! grep -q "WARNING.*density file" chk.log
printf '# z[k0^-1] n\n0 0\n180 1\n' > short.txt
withfile short.txt > c.deck
expect_warn "file ends inside the run" "ends at 0.0229183 mm with 100 % of its peak" c.deck
"$EXE" c.deck plasma.file_outside=hold --check > chk.log 2>&1
grep -q "edge values are continued" chk.log && ! grep -q "WARNING.*density file" chk.log
printf '# z[k0^-1] n\n0 0\n100 0\n120 1\n400 1\n' > steep.txt
withfile steep.txt > c.deck
expect_warn "steep gradient" "less than 4 steps" c.deck
printf '# z[k0^-1] n\n0 0\n100 -1\n' > bad.txt; withfile bad.txt > c.deck
expect_fail "negative density" "line 3: negative density" c.deck
printf '# z[k0^-1] n\n0 0\n100 nan\n' > bad.txt
expect_fail "nan" "not finite" c.deck
printf '# z[k0^-1] n\n0 0\n100 1\n100 1\n' > bad.txt
expect_fail "z not increasing" "strictly increasing" c.deck
printf '# z[k0^-1] r[k0^-1] n\n0 0 1\n0 10 1\n100 0 1\n' > bad.txt
expect_fail "incomplete grid" "do not fill a rectilinear grid" c.deck
printf '0 0\n100 1\n' > bad.txt
expect_fail "no header" "no header" c.deck
printf '# z[fs] n\n0 0\n100 1\n' > bad.txt
expect_fail "wrong unit" "is a time" c.deck
printf '# z n\n0 0\n100 1\n' > bad.txt
expect_fail "no length unit" "needs a unit" c.deck
expect_fail "missing file" "cannot open" c.deck plasma.density_file=nothing_here.txt
expect_fail "z_points with the file" "given by density_file" f1.deck plasma.z_points=0,1mm plasma.z_values=0,1
expect_fail "radial_profile with n(r, z)" "includes the radial profile" f3.deck plasma.radial_profile=uniform
expect_fail "file outside the run" "outside the run" f1.deck plasma.file_z=1mm
expect_fail "file_z without a file" "only used with 'density_file'" "$R" plasma.file_z=1mm
expect_fail "absolute file and ion density" "remove 'density'" i2.deck ion.density=0.002n_c
echo "  summary, 2 warnings, 12 rejected files and combinations  ok"

# ---- 7. the example (examples/gas_jet): the files are found beside the deck; no warning about the profile
G="$SRC/examples/gas_jet/gas_jet.deck"
"$EXE" "$G" --check > ex1.log 2>&1
grep -q "n(z) on 41 points, absolute density \[cm^-3\]" ex1.log
grep -q "peak 3e+18 cm^-3" ex1.log
"$EXE" "$G" plasma.density_file=jet_rz.txt --check > ex2.log 2>&1
grep -q "n(r, z) on 6 x 41 points, absolute density \[cm^-3\]" ex2.log
grep -q "peak 4.406e+18 cm^-3 .* r = 25 um" ex2.log
! grep -q "WARNING.*density file" ex1.log && ! grep -q "WARNING.*density file" ex2.log
"$EXE" "$G" numerics.backend=$BACKEND grid.z_end=15.27um output.every=0 output.dir=ex > ex3.log
awk '!/^#/ { for (c = 1; c <= NF; c++) if ($c != $c + 0 || $c == "nan" || $c == "inf") bad = 1; n++ }
  END { exit !(n > 0 && !bad) }' ex/axis_00003.txt
echo "  examples/gas_jet: n(z) and n(r, z) sample files accepted without warnings, 3 steps run  ok"
cd - > /dev/null
[ -n "$KEEP" ] && echo "$WORK" || rm -rf "$WORK"
echo "density_file test passed ($BACKEND)"
