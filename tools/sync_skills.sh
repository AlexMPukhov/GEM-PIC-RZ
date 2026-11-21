#!/bin/sh
# Agent skills (skills/, see docs/agent_skill.md): the parts that are copies of repository files.
#   sh tools/sync_skills.sh           regenerate the copies (run it after editing docs/input_deck.md or the examples)
#   sh tools/sync_skills.sh --check   exit 1 if a copy differs from its source (ctest agent_skill)
# A skill is installed without the rest of the repository, so it carries its own copy of the deck reference and of
# the example decks; the check keeps them identical to the sources.
SRC=$(cd "$(dirname "$0")/.." && pwd)
SK="$SRC/skills/gempic-rz-simulation"
URL=https://github.com/AlexMPukhov/GEM-PIC-RZ/blob/main
DECKS="lwfa_witness pwfa_driver_witness fel_sase undulator_emission"
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

# the deck reference: a header line, links to other repository files made absolute
{
  echo "<!-- Generated from docs/input_deck.md by tools/sync_skills.sh. Do not edit this copy. -->"
  sed -e "s|](openpmd_rectilinear.md)|]($URL/docs/openpmd_rectilinear.md)|g" \
      -e "s|](density_file.md)|]($URL/docs/density_file.md)|g" \
      -e "s|](\.\./examples/|]($URL/examples/|g" "$SRC/docs/input_deck.md"
} > "$TMP/input_deck.md"

rc=0
put() {  # <generated file> <destination>
  if [ "$MODE" = "--check" ]; then
    if cmp -s "$1" "$2"; then echo "  in sync: ${2#$SRC/}"
    else echo "  OUT OF SYNC: ${2#$SRC/}  (run: sh tools/sync_skills.sh)"; rc=1; fi
  else
    mkdir -p "$(dirname "$2")" && cp "$1" "$2" && echo "  written: ${2#$SRC/}"
  fi
}
MODE=$1
put "$TMP/input_deck.md" "$SK/references/input_deck.md"
for d in $DECKS; do put "$SRC/examples/$d.deck" "$SK/examples/$d.deck"; done
exit $rc
