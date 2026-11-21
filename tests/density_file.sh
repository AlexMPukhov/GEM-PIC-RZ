#!/bin/sh
# Density profiles from files (plasma.density_file, ion.density_file; format: docs/density_file.md) against the same
# profiles given in the deck (tests/regression.deck: laser, plasma up-ramp 0 -> 1 over 400 k0^-1, H+ ions, witness
# bunch; 4 steps):
#   1. geometry z, columns layout, a factor on `density`, with the points of the deck's z_points / z_values:
#      identical to the inline table (serial, openmp; cuda: to rounding)
#   2. the same ramp as absolute density in cm^-3 over um (other units: the numbers differ in the last digit), with
#      commas, in the grid layout (uniform axis) and shifted with file_z: equal to rounding
#   3. geometry rz (ramp in z x a hollow radial profile, linear between the grid points), columns layout and grid
#      layout (list axes, both orders), against the same profile given as density_function(z, r); restart from a
#      checkpoint = uninterrupted run
#   4. geometry z with a transverse function in the header (Gaussian in r; one that depends on x) against the
#      inline table times the same density_function
#   5. openPMD meshes (tests/density_files.py; needs python3 with numpy and h5py, else skipped): the profile of 3. on
#      a regular thetaMode mesh, on a rectilinear mesh (gridCoordinates) as absolute density in m^-3, a mesh n(z)
#      with the attribute transverseFunction, and a smooth gas jet sampled on 321 x 161 points against its
#      density_function (axis fields to 1e-3 of the peak)
#   6. an ion species with its own profile from a file (pre-ionized H, free_electrons = false) = the inline table
#   7. start-up checks: summary (ranges, peak), warnings (file cut at a finite density, gradient shorter than
#      4 steps), rejected files (no geometry, 3D and Cartesian geometries, negative, nan, z not increasing,
#      incomplete grid, wrong count, wrong unit, unknown key) and rejected combinations (z_points or radial_profile
#      with the file, file outside the run)
#   8. the example examples/gas_jet
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
# geometry = z
# columns  = z[k0^-1]  n
0    0
400  1
EOF
withfile ramp.txt > f1.deck
grep -q density_file f1.deck && ! grep -q z_values f1.deck
"$EXE" "$R" $O output.dir=ref > ref.log
"$EXE" f1.deck $O output.dir=f1 > f1.log
same ref f1 $T0
grep -q "density file ramp.txt (text; geometry z, columns): n(z) on 2 points, factor on 'density'" f1.log
echo "  n(z), columns layout (factor) = inline z_points / z_values: tolerances $T0  ok"

# ---- 2. absolute density, other units, commas, grid layout, file_z
# 0.002 n_c at 0.8 um = 0.002 * 1.1148e21 / 0.64 cm^-3; 400 k0^-1 = 400 * 0.8 / (2 pi) um
awk 'BEGIN { n = 0.002 * 1.1148e21 / 0.64; z = 400 * 0.8 / (2 * 3.14159265358979323846)
  printf "# geometry = z\n# columns = z[um], n_e[cm^-3]\n0,0\n%.17g,%.17g\n%.17g,%.17g\n", z, n, 3 * z, n > "ramp_abs.csv"
  printf "# density_profile = 1\n# geometry = z\n# density_unit = m^-3\n# z = uniform 1 %.17g 4 mm\n0 %.17g\n%.17g %.17g\n", z * 1e-3, n * 1e6, n * 1e6, n * 1e6 > "ramp_grid.txt" }'
withfile ramp_abs.csv > f2.deck
"$EXE" f2.deck $O output.dir=f2 > f2.log
same ref f2 $T1
grep -q "absolute density \[cm^-3\]" f2.log
grep -q "peak 3.484e+18 cm^-3 (0.002 n_c)" f2.log
withfile ramp_grid.txt > f2s.deck
"$EXE" f2s.deck $O plasma.file_z=-1mm output.dir=f2s > f2s.log
same ref f2s $T1
grep -q "geometry z, grid layout): n(z) on 4 points, absolute density \[m^-3\]" f2s.log
echo "  n(z) in cm^-3 over um (columns) and in m^-3 over mm (grid layout, file_z) = inline table: tolerances $T1  ok"

# ---- 3. n(r, z) against density_function
# f(z) = min(z / 400, 1), g(r) = 0.6 + 0.8 min(r / 75, 1): linear between the grid points of the file
awk 'function f(z) { return z < 400 ? z / 400 : 1 } function g(r) { return 0.6 + 0.8 * (r < 75 ? r / 75 : 1) }
BEGIN { nz = split("0 400 1000", z, " "); nr = split("0 75 150 400", r, " ")
  a = "jet2d.txt"; b = "jet2d_zr.txt"; c = "jet2d_rz.txt"
  print "# geometry = rz\n# columns = z[k0^-1] r[k0^-1] n" > a
  h = "# geometry = rz\n# density_unit = 1\n# z = list 0 400 1000 k0^-1\n# r = list 0 75 150 400 k0^-1"
  print h "\n# order = z r" > b; print h "\n# order = r z" > c
  for (i = 1; i <= nz; i++) { for (j = 1; j <= nr; j++) {
      printf "%g %g %.17g\n", z[i], r[j], f(z[i]) * g(r[j]) > a; printf "%.17g ", f(z[i]) * g(r[j]) > b }
    printf "\n" > b }
  for (j = 1; j <= nr; j++) { for (i = 1; i <= nz; i++) printf "%.17g\n", f(z[i]) * g(r[j]) > c } }'
withfile jet2d.txt > f3.deck
sed -e 's|^  z_points .*|  density_function = "min(z / (400 k0^-1), 1) * (0.6 + 0.8 * min(r / (75 k0^-1), 1))"|' \
    -e '/^  z_values /d' "$R" > fn3.deck
"$EXE" fn3.deck $O output.dir=fn3 > fn3.log
"$EXE" f3.deck $O output.dir=f3 > f3.log
same fn3 f3 $T1
grep -q "n(r, z) on 4 x 3 points" f3.log
grep -q "radial profile: from the density file" f3.log
for l in zr rz; do
  "$EXE" f3.deck $O plasma.density_file=jet2d_$l.txt output.dir=f3$l > f3$l.log
  same f3 f3$l $T0
done
"$EXE" f3.deck $O output.dir=f3r grid.z_end=80k0^-1 output.checkpoint_every=2 > f3r1.log
"$EXE" f3.deck $O output.dir=f3r --restart latest > f3r2.log
for f in axis_00004.txt map_00004.txt bunches.txt; do "$CMP" f3/$f f3r/$f $T0; done
echo "  n(r, z), columns layout = density_function(z, r): tolerances $T1; grid layout in both orders = columns"
echo "  layout and restart = uninterrupted run ($T0)  ok"

# ---- 4. n(z) with a transverse function in the header
for c in gauss side; do
  if [ $c = gauss ]; then T='exp(-(r / (120 k0^-1))^2)'; else T='(1 + 0.3 * x / (150 k0^-1)) * exp(-(r / (120 k0^-1))^2)'; fi
  { printf '# geometry = z\n# columns = z[k0^-1] n\n# transverse = "%s"\n0 0\n400 1\n' "$T"; } > t_$c.txt
  withfile t_$c.txt > ft.deck
  sed -e "s|^  r_max  |  density_function = \"$T\"\n  r_max  |" "$R" > fnt.deck
  grep -q density_function fnt.deck
  "$EXE" fnt.deck $O output.dir=fnt_$c > fnt_$c.log
  "$EXE" ft.deck $O output.dir=ft_$c > ft_$c.log
  same fnt_$c ft_$c $T1
  grep -q "times the transverse function of the file" ft_$c.log
done
echo "  n(z) x transverse function of the header (Gaussian in r; depending on x) = inline table x density_function:"
echo "  tolerances $T1  ok"

# ---- 5. openPMD meshes
if python3 -c "import numpy, h5py" 2>/dev/null; then
  python3 "$SRC/tests/density_files.py"
  withfile jet2d_theta.h5 > f4.deck
  "$EXE" f4.deck $O output.dir=f4a > f4a.log
  same fn3 f4a $T1
  grep -q "openPMD; mesh 'density', iteration 0, geometry rz" f4a.log
  "$EXE" f4.deck $O plasma.density_file=jet2d_rect.h5 plasma.file_mesh=n_gas output.dir=f4b > f4b.log
  same fn3 f4b $T1
  grep -q "absolute density \[m^-3\]" f4b.log
  "$EXE" f4.deck $O plasma.density_file=ramp_gauss.h5 output.dir=f4t > f4t.log
  same fnt_gauss f4t $T1
  # smooth gas jet: linear interpolation of 321 x 161 samples against the function
  JET='exp(-((z - 150 k0^-1) / (200 k0^-1))^2 - (r / (120 k0^-1))^2)'
  sed -e "s|^  z_points .*|  density_function = \"$JET\"|" -e '/^  z_values /d' "$R" > fnj.deck
  "$EXE" fnj.deck $O output.dir=fnj > fnj.log
  "$EXE" f4.deck $O plasma.density_file=gauss_jet.h5 output.dir=f4j > f4j.log
  paste fnj/axis_00004.txt f4j/axis_00004.txt | awk '!/^#/ { h = NF / 2
      for (c = 2; c <= h; c++) { d = $c - $(c + h); if (d < 0) d = -d; if (d > m[c]) m[c] = d
        a = $c < 0 ? -$c : $c; if (a > p[c]) p[c] = a } }
    END { for (c = 2; c <= 4; c++) { r = m[c] / p[c]; if (r > w) w = r }
      printf "  smooth jet, 321 x 161 samples: deviation on the axis / peak: E_x %.2e, E_z %.2e, n_e %.2e (tolerance 1e-3)", m[2] / p[2], m[3] / p[3], m[4] / p[4]
      if (w <= 1e-3) print "  ok"; else { print "  FAILED"; exit 1 } }'
  echo "  openPMD: regular thetaMode mesh, rectilinear mesh in m^-3 = density_function; n(z) with transverseFunction ="
  echo "  the text file: tolerances $T1  ok"
else
  echo "  openPMD part skipped: needs python3 with numpy and h5py"
fi

# ---- 6. ion block with its own profile from a file
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

# ---- 7. start-up checks
nfail=0
expect_fail() {  # <label> <text expected in the message> <deck> [overrides]
  l=$1; t=$2; shift 2
  if "$EXE" "$@" --check > chk.log 2>&1; then echo "  $l: accepted, expected an error  FAILED"; exit 1; fi
  grep -q -- "$t" chk.log || { echo "  $l: error message without '$t':"; cat chk.log; exit 1; }
  nfail=$((nfail + 1))
}
expect_warn() {  # <label> <text> <deck> [overrides]
  l=$1; t=$2; shift 2
  "$EXE" "$@" --check > chk.log 2>&1 || { echo "  $l: rejected:"; cat chk.log; exit 1; }
  grep "WARNING" chk.log | grep -q -- "$t" || { echo "  $l: no warning with '$t':"; cat chk.log; exit 1; }
}
H='# geometry = z\n# columns = z[k0^-1] n\n'
"$EXE" f1.deck --check > chk.log 2>&1
grep -q "z = 0 ... 0.0509296 mm (run: 0 ... 0.0254648 mm)" chk.log
grep -q "peak 3.484e+18 cm^-3 (0.002 n_c) at z = 0.0509296 mm" chk.log
! grep -q "WARNING.*density file" chk.log
printf "${H}0 0\n180 1\n" > short.txt
withfile short.txt > c.deck
expect_warn "file ends inside the run" "ends at 0.0229183 mm with 100 % of its peak" c.deck
"$EXE" c.deck plasma.file_outside=hold --check > chk.log 2>&1
grep -q "edge values are continued" chk.log && ! grep -q "WARNING.*density file" chk.log
printf "${H}0 0\n100 0\n120 1\n400 1\n" > steep.txt
withfile steep.txt > c.deck
expect_warn "steep gradient" "less than 4 steps" c.deck
withfile bad.txt > c.deck
printf '# z[k0^-1] n\n0 0\n100 1\n' > bad.txt
expect_fail "no geometry" "no line '# geometry = ...'" c.deck
printf '# geometry = xyz\n# columns = x[um] y[um] z[um] n\n0 0 0 1\n' > bad.txt
expect_fail "3D file" "3D profile n(x, y, z) is for the 3D codes" c.deck
printf '# geometry = xz\n# columns = x[um] z[um] n\n0 0 1\n' > bad.txt
expect_fail "Cartesian 2D file" "Cartesian 2D profile" c.deck
printf '# geometry = z\n0 0\n100 1\n' > bad.txt
expect_fail "no layout" "neither '# columns = ...' nor the axes" c.deck
printf '# geometry = rz\n# columns = z[k0^-1] n\n0 0\n100 1\n' > bad.txt
expect_fail "columns do not fit the geometry" "does not fit geometry rz" c.deck
printf "${H}# presure = 1 bar\n0 0\n100 1\n" > bad.txt
expect_fail "unknown key" "unknown header key 'presure'" c.deck
printf "${H}0 0\n100 -1\n" > bad.txt
expect_fail "negative density" "line 4: negative density" c.deck
printf "${H}0 0\n100 nan\n" > bad.txt
expect_fail "nan" "not finite" c.deck
printf "${H}0 0\n100 1\n100 1\n" > bad.txt
expect_fail "z not increasing" "strictly increasing" c.deck
printf '# geometry = rz\n# columns = z[k0^-1] r[k0^-1] n\n0 0 1\n0 10 1\n100 0 1\n' > bad.txt
expect_fail "incomplete grid" "do not fill a rectilinear grid" c.deck
printf '# geometry = z\n# density_unit = 1\n# z = uniform 0 100 3 k0^-1\n0 1\n' > bad.txt
expect_fail "wrong number of values" "2 values, but the axes of the header give 3" c.deck
printf '# geometry = z\n# columns = z[fs] n\n0 0\n100 1\n' > bad.txt
expect_fail "wrong unit" "is a time" c.deck
printf '# geometry = z\n# columns = z n\n0 0\n100 1\n' > bad.txt
expect_fail "no length unit" "needs a unit" c.deck
printf '# geometry = rz\n# columns = z[k0^-1] r[k0^-1] n\n# transverse = "1"\n0 0 1\n0 10 1\n100 0 1\n100 10 1\n' > bad.txt
expect_fail "transverse function with rz" "only for geometry z" c.deck
expect_fail "missing file" "cannot open" c.deck plasma.density_file=nothing_here.txt
expect_fail "z_points with the file" "given by density_file" f1.deck plasma.z_points=0,1mm plasma.z_values=0,1
expect_fail "radial_profile with n(r, z)" "includes the radial profile" f3.deck plasma.radial_profile=uniform
expect_fail "radial_profile with a transverse function" "carries a transverse function" ft.deck plasma.radial_profile=uniform
expect_fail "file outside the run" "outside the run" f1.deck plasma.file_z=1mm
expect_fail "file_z without a file" "only used with 'density_file'" "$R" plasma.file_z=1mm
expect_fail "absolute file and ion density" "remove 'density'" i2.deck ion.density=0.002n_c
echo "  summary, 2 warnings, $nfail rejected files and combinations  ok"

# ---- 8. the example (examples/gas_jet): the files are found beside the deck; no warning about the profile
G="$SRC/examples/gas_jet/gas_jet.deck"
"$EXE" "$G" --check > ex1.log 2>&1
grep -q "n(z) on 41 points, absolute density \[cm^-3\]" ex1.log
grep -q "peak 3e+18 cm^-3" ex1.log
"$EXE" "$G" plasma.density_file=jet_rz.txt --check > ex2.log 2>&1
grep -q "n(r, z) on 6 x 41 points, absolute density \[cm^-3\]" ex2.log
grep -q "peak 4.406e+18 cm^-3 .* r = 25 um" ex2.log
"$EXE" "$G" plasma.density_file=jet_gauss.txt --check > ex3.log 2>&1
grep -q "times the transverse function of the file" ex3.log
! grep -q "WARNING.*density file" ex1.log && ! grep -q "WARNING.*density file" ex2.log
"$EXE" "$G" numerics.backend=$BACKEND grid.z_end=15.27um output.every=0 output.dir=ex > ex4.log
awk '!/^#/ { for (c = 1; c <= NF; c++) if ($c != $c + 0 || $c == "nan" || $c == "inf") bad = 1; n++ }
  END { exit !(n > 0 && !bad) }' ex/axis_00003.txt
echo "  examples/gas_jet: the three sample files accepted, no warnings for the profile, 3 steps run  ok"
cd - > /dev/null
[ -n "$KEEP" ] && echo "$WORK" || rm -rf "$WORK"
echo "density_file test passed ($BACKEND)"
