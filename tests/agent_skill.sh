#!/bin/sh
# Agent skills (skills/*/SKILL.md, docs/agent_skill.md) against the code, so that they cannot drift from it:
#   1. format of the open Agent Skills specification (agentskills.io): frontmatter name = directory name
#      (lower case, digits, hyphens, <= 64), description 1 ... 1024 characters, SKILL.md <= 500 lines; every file
#      of the skill that a SKILL.md names (references/, scripts/, examples/) exists
#   2. the copies of docs/input_deck.md and of the example decks equal their sources (tools/sync_skills.sh --check)
#   3. every deck key the hand-written skill files name as block.key is a key the parser reads (string literal in
#      src/config.cpp) and is documented for that block in docs/input_deck.md
#   4. every ```deck snippet of the skill files passes  gempic_rz --check
#   5. every deck shipped in a skill passes --check and runs 3 steps with the real executable (finite numbers in
#      bunches.txt and in the axis files); the LASY deck with the file written by the skill's own script
#   6. the helper scripts: plasma_numbers.py against lambda_p of the code's summary
# usage: sh tests/agent_skill.sh <gempic_rz> [hdf5: 1 | 0]
#
# Reuse in another code repository (QUARZ, ...): copy this file and tools/sync_skills.sh and change only the
# block "code-specific settings" below (and the list of copied files in tools/sync_skills.sh).
abspath() { (cd "$(dirname "$1")" && echo "$(pwd)/$(basename "$1")"); }
EXE=$(abspath "$1")
HDF5=${2:-1}
SRC=$(cd "$(dirname "$0")/.." && pwd)
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT
cd "$WORK"
fail=0
bad() { echo "  FAILED: $*"; fail=1; }
# ---- code-specific settings -----------------------------------------------------------------------------------
# blocks of the input deck: a word of the skill files written as block.key is taken as a deck key
BLOCKS="grid laser flying_focus plasma ion bunch undulator quadrupole plasma_lens solenoid chicane wall dielectric bound_electrons output numerics"
PARSER="$SRC/src/config.cpp"          # the source that reads the keys: "key" must occur there as a string literal
REFERENCE="$SRC/docs/input_deck.md"   # the input reference: sections "## ... `block` ..." must name `key`
MINSKILLS=2; MINKEYS=30; MINDECKS=4   # fewer found = the extraction is broken
check_deck() { "$EXE" "$@" --check; }                       # <deck> [overrides]: parse and validate only
# <deck> <output dir> [overrides]: a run of 3 steps (z_end = 3 dz, with dz and its unit taken from the deck)
short_run() {
  deck=$1; out=$2; shift 2
  zend=$(awk '/^[ \t]*dz[ \t]*=/ { printf "%g%s", 3 * $3, $4; exit }' "$deck")
  [ -n "$zend" ] || { echo "no uniform grid.dz in the deck"; return 1; }
  OMP_PROC_BIND=false "$EXE" "$deck" grid.z_end="$zend" output.dir="$out" output.every=0 "$@"
}
# ---- end of the code-specific settings ------------------------------------------------------------------------
PY=""
python3 -c "import sys" 2>/dev/null && PY=python3

echo "1. format (agentskills.io specification)"
nskill=0
for dir in "$SRC"/skills/*/; do
  dir=${dir%/}; name=$(basename "$dir"); f="$dir/SKILL.md"; nskill=$((nskill + 1))
  [ -f "$f" ] || { bad "$name: no SKILL.md"; continue; }
  [ "$(head -1 "$f")" = "---" ] || bad "$name: SKILL.md must start with the YAML frontmatter (---)"
  n=$(awk '/^---$/ { c++; next } c == 1 && /^name:/ { sub(/^name:[ \t]*/, ""); print; exit }' "$f")
  d=$(awk '/^---$/ { c++; next } c == 1 && /^description:/ { sub(/^description:[ \t]*/, ""); print; exit }' "$f")
  [ "$n" = "$name" ] || bad "$name: frontmatter name '$n' differs from the directory name"
  echo "$n" | grep -Eq '^[a-z0-9]+(-[a-z0-9]+)*$' || bad "$name: name must be lower case letters, digits, single hyphens"
  [ ${#n} -le 64 ] || bad "$name: name longer than 64 characters"
  [ ${#d} -ge 1 ] && [ ${#d} -le 1024 ] || bad "$name: description has ${#d} characters (1 ... 1024)"
  l=$(wc -l < "$f"); [ "$l" -le 500 ] || bad "$name: SKILL.md has $l lines (<= 500)"
  for ref in $(grep -oE '(references|scripts|examples)/[A-Za-z0-9_.-]+\.[a-z0-9]+' "$f" | sort -u); do
    [ -f "$dir/$ref" ] || bad "$name: SKILL.md names $ref, which is not in the skill"
  done
  echo "  $name: name ok, description ${#d} characters, $l lines"
done
[ "$nskill" -ge $MINSKILLS ] || bad "expected at least $MINSKILLS skills in skills/, found $nskill"

echo "2. copies of repository files"
sh "$SRC/tools/sync_skills.sh" --check || fail=1

echo "3. deck keys named in the skill files"
# hand-written files only (references/input_deck.md is the copy checked in 2.)
MD=$(ls "$SRC"/skills/*/SKILL.md "$SRC"/skills/*/references/*.md | grep -v '/input_deck.md$')
alt=$(echo $BLOCKS | tr ' ' '|')
cat $MD | grep -oE "(^|[^A-Za-z0-9_.])($alt)(\[[0-9]+\])?\.[a-z_][a-z_0-9]*" |
  sed -E 's/^[^a-z]*//; s/\[[0-9]+\]//' | grep -vE '\.(deck|txt|h5|md|sh|py|cpp|hpp)$' | sort -u > keys.txt
nkey=0
while IFS=. read -r block key; do
  nkey=$((nkey + 1))
  grep -q "\"$key\"" "$PARSER" || bad "$block.$key: the parser (${PARSER#$SRC/}) reads no key \"$key\""
  # the section of the block in docs/input_deck.md: from its "## " heading to the next one
  awk -v b="\`$block\`" -v k="\`$key\`" '/^## / { on = index($0, b) > 0 } on && index($0, k) { found = 1 }
    END { exit !found }' "$REFERENCE" || bad "$block.$key: not documented for the block $block in ${REFERENCE#$SRC/}"
done < keys.txt
echo "  $nkey keys checked against ${PARSER#$SRC/} and ${REFERENCE#$SRC/}"
[ "$nkey" -ge $MINKEYS ] || bad "only $nkey keys found in the skill files: the extraction is broken"

echo "4. deck snippets of the skill files"
nsnip=0
for f in $MD; do
  awk -v base="snip_$(basename "$(dirname "$f")")_$(basename "$f" .md)" '
    /^```deck[ \t]*$/ { n++; out = base "_" n ".deck"; next }
    /^```/ { out = ""; next }
    out != "" { print > out }' "$f"
done
for s in snip_*.deck; do
  [ -f "$s" ] || continue
  nsnip=$((nsnip + 1))
  if check_deck "$s" > "$s.log" 2>&1; then echo "  $s: ok"
  else bad "$s does not pass --check:"; tail -3 "$s.log"; fi
done
[ "$nsnip" -ge 1 ] || bad "no \`\`\`deck snippet found in the skill files"

echo "5. decks shipped in the skills: --check and 3 steps"
ndeck=0
for deck in "$SRC"/skills/*/examples/*.deck; do
  b=$(basename "$deck" .deck); ndeck=$((ndeck + 1))
  scripts="$(dirname "$(dirname "$deck")")/scripts"
  if grep -q 'lasy_gaussian.h5' "$deck"; then          # needs the envelope file of the skill's script
    if [ "$HDF5" != 1 ] || ! python3 -c "import numpy, h5py" 2>/dev/null; then
      echo "  $b: skipped (needs HDF5 and python3 with numpy and h5py)"; continue
    fi
    python3 "$scripts/make_lasy_gaussian.py" > "$b.make.log" 2>&1 || { bad "$b: make_lasy_gaussian.py"; continue; }
  fi
  extra=""
  [ "$HDF5" = 1 ] || extra="output.format=text"
  check_deck "$deck" $extra > "$b.check.log" 2>&1 || { bad "$b: --check"; tail -3 "$b.check.log"; continue; }
  if ! short_run "$deck" "out_$b" $extra > "$b.run.log" 2>&1; then
    bad "$b: the run stopped"; tail -5 "$b.run.log"; continue
  fi
  steps=$(grep -c '^step' "$b.run.log")
  grep -q '^done' "$b.run.log" || bad "$b: no 'done' line"
  [ -s "out_$b/bunches.txt" ] || bad "$b: no bunches.txt"
  [ "$HDF5" = 1 ] && grep -q 'hdf5' "$deck" || ls "out_$b"/axis_*.txt > /dev/null 2>&1 || bad "$b: no axis files"
  if cat "out_$b"/bunches.txt "out_$b"/axis_*.txt 2>/dev/null | grep -v '^#' | grep -qiE 'nan|inf'; then
    bad "$b: NaN or Inf in the output"
  fi
  echo "  $b: ok ($(grep '^done' "$b.run.log" | cut -c1-40))"
done
[ "$ndeck" -ge $MINDECKS ] || bad "expected at least $MINDECKS example decks in the skills, found $ndeck"

echo "6. helper scripts"
if [ -n "$PY" ]; then
  S="$SRC/skills/gempic-rz-simulation/scripts"
  for s in "$S"/*.py; do
    python3 -c "import ast, sys; ast.parse(open(sys.argv[1]).read())" "$s" || bad "$(basename "$s"): syntax"
  done
  # lambda_p of plasma_numbers.py against the summary of the code (3.48e18 cm^-3: 17.9 um)
  lp=$(python3 "$S/plasma_numbers.py" --density 3.48e18 --a0 2 | awk '/^plasma wavelength/ { print $4 }')
  ref=$(awk '/lambda_p =/ { for (i = 1; i < NF; i++) if ($i == "lambda_p") print $(i + 2); exit }' lwfa_witness.check.log)
  if awk -v a="$lp" -v b="$ref" 'BEGIN { d = (a - b) / b; if (d < 0) d = -d; exit !(b > 0 && d < 2e-3) }'; then
    echo "  plasma_numbers.py: lambda_p = $lp um (code: $ref um)  ok"
  else bad "plasma_numbers.py: lambda_p = $lp um, the code's summary has $ref um"; fi
else
  echo "  skipped (no python3)"
fi

if [ $fail = 0 ]; then echo "agent skill test passed"; else echo "agent skill test FAILED"; fi
exit $fail
