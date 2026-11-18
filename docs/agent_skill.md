# Agent skills for GEM-PIC-RZ

An *agent skill* is a folder with a `SKILL.md` file (and optional references, scripts and examples) that an AI
coding agent loads when a task calls for it. The repository ships two, so that an agent can build GEM-PIC-RZ and
write, check, run and debug input decks from a request in plain words.

## Format and where it works

The skills follow the open [Agent Skills specification](https://agentskills.io/specification):

- a directory named like the skill, with `SKILL.md`: YAML frontmatter with `name` (lower case, digits, hyphens,
  ≤ 64 characters, equal to the directory name) and `description` (≤ 1024 characters: what the skill does and when
  to use it), then the instructions in Markdown;
- optional `references/` (read on demand), `scripts/` (run by the agent), and other files (`examples/` here);
- *progressive disclosure*: the agent always sees only name and description; it loads `SKILL.md` (kept below 500
  lines) when the task matches, and the other files only when `SKILL.md` sends it there.

The format is not tied to one vendor or model. It is read natively by Claude Code, OpenCode, Codex, Cursor,
Gemini CLI, Goose, Cline, Qwen Code, OpenHands, Crush, Roo Code, Kilo Code, Continue and others, so it works with
open-source agents and open models. An agent without skill support can still use it: tell it to read
`skills/gempic-rz-simulation/SKILL.md` and follow it.

## Install

```
npx skills add AlexMPukhov/GEM-PIC-RZ                       # asks which skills and which agents
npx skills add AlexMPukhov/GEM-PIC-RZ -a opencode -g        # one agent, for all projects of the user
npx skills add AlexMPukhov/GEM-PIC-RZ --skill gempic-rz-simulation
```

`skills` is the open-source installer of [skills.sh](https://skills.sh) (needs Node.js). It looks for
`skills/<name>/SKILL.md` in the repository, which is where they are. Without Node.js copy the folders:

| agent | project directory | user directory |
|---|---|---|
| Claude Code | `.claude/skills/` | `~/.claude/skills/` |
| OpenCode | `.agents/skills/` | `~/.config/opencode/skills/` |
| Codex, Cline | `.agents/skills/` | `~/.agents/skills/` |
| Cursor | `.agents/skills/` | `~/.cursor/skills/` |
| Gemini CLI | `.agents/skills/` | `~/.gemini/skills/` |
| Goose | `.goose/skills/` | `~/.config/goose/skills/` |

```
cp -r GEM-PIC-RZ/skills/gempic-rz-* ~/.claude/skills/        # e.g. Claude Code, for all projects
```

Then ask for a simulation: "LWFA, a0 = 2, 800 nm, n = 1e18 cm^-3, 2 mm", "a seeded FEL at 13 nm with a 1 GeV, 3 kA
beam", "why does my witness lose charge?".

## The two skills

```
skills/
  gempic-rz-setup/
    SKILL.md                    clone, build, ctest, backends (serial, OpenMP, MPI, CUDA), output, cluster jobs
  gempic-rz-simulation/
    SKILL.md                    deck in one page, 8 rules, "verify before delivering"
    references/input_deck.md    all blocks, keys, units, defaults         (copy of docs/input_deck.md)
    references/physics_rules.md how to choose grid, box and steps; what the code can and cannot model
    references/diagnostics.md   output files, columns, units, plotting, what to look at
    references/debugging.md     errors, warnings, symptoms and fixes
    examples/*.deck             lwfa_witness, pwfa_driver_witness, fel_sase, undulator_emission  (copies of examples/)
                                lwfa_plain (plain LWFA, all numbers from plasma_numbers.py)
                                lasy_laser (laser from a LASY file, focus 0.3 mm inside the run)
    scripts/plasma_numbers.py   λ_p, c/ω_p, P/P_c, matched spot; starting values of grid, box, ramp, focus
    scripts/plot_run.py         table per output step (numpy); quick-look figure (matplotlib)
    scripts/make_lasy_gaussian.py  a Gaussian pulse as a LASY envelope file (numpy, h5py)
```

What the simulation skill makes an agent do: copy the closest example, change the physics, choose the numerical
parameters from the rules (carrier resolved, window, `dz` as an accuracy parameter, a radial grid that does not
reflect radiation, the box, the FEL resolution, the limits of the model, e.g. density down-ramp injection belongs to
a full PIC code), and then verify: `gempic_rz deck --check`, every warning fixed or explained, a cost estimate, a
smoke run of a few steps, the full run, numbers reported.

## Keeping the skills in step with the code

A skill is installed without the rest of the repository, so it carries copies. Two mechanisms keep them right:

- `sh tools/sync_skills.sh` regenerates `references/input_deck.md` from `docs/input_deck.md` and the four example
  decks from `examples/`. Run it after editing those files.
- ctest `agent_skill` (`tests/agent_skill.sh`, 25 s):
  1. the format of the specification (names, description length, `SKILL.md` ≤ 500 lines, every file a `SKILL.md`
     names exists);
  2. the copies equal their sources (`tools/sync_skills.sh --check`);
  3. every deck key that the hand-written skill files name as `block.key` (62 keys) is read by the parser (string
     literal in `src/config.cpp`) and documented for that block in `docs/input_deck.md`;
  4. every ```` ```deck ```` snippet of the skill files passes `gempic_rz --check`;
  5. every deck of the skills passes `--check` and runs 3 steps with the real executable, with finite numbers in
     `bunches.txt` and the axis files (the LASY deck with the file written by the skill's script; skipped without
     HDF5 or h5py);
  6. `plasma_numbers.py` gives the plasma wavelength of the code's summary.

  With a wrong key (`grid.dx_fake`), an edited copy or a snippet with a renamed key the test fails and names the
  place.

Rules for editing the skills: write deck keys as `block.key` (so that the test sees them); put a complete small
deck into a ```` ```deck ```` block if it should be checked; keep `SKILL.md` short and move details to
`references/`; the repository is public, so only published material goes in.

## Trial with a fresh agent

Two fresh agent sessions (2026-11-18; a large and a small Claude model, run as sub-agents) got the two skill
folders, an existing build and the one-line request "LWFA, a0 = 2, 800 nm, n = 1e18 cm^-3, 2 mm", with the
instruction to use nothing but the skill files. They could not ask questions.

| | large model | small model |
|---|---|---|
| `--check` attempts until the deck was clean | 2 (two warnings at the first) | 1 |
| errors of the code | none | none |
| smoke run, `\|Ex\|max` at step 0 | 1.987 (= a0/√(1 + z_focus²/z_R²), focus 100 µm ahead) | 2.002 |
| full run, 209 steps, 3125 slices × 85–90 radial cells | 394 s on 8 threads | 428 s on 6 threads |
| laser amplitude at 0.5 mm → 2 mm | 2.22 → 1.44 | 2.10 → 1.38 |
| accelerating field of the first bucket at 0.5 mm | 0.88 E0 (85 GV/m) | 0.84 E0 (81 GV/m) |
| trapped charge | 0 | 0 |

Both decks were valid and ran to the end at the first full attempt; both agents chose the matched spot (15 µm), a
100 µm up-ramp, a uniform radial grid out to the absorber, `plasma.smooth_length` and `numerics.trap_gamma`
without being told, and neither added a witness. What the first version of the skill got wrong, by their reports,
and what was changed:

| finding | change |
|---|---|
| the helper script proposed a box of 2.5 w0 without regard to the run length, while the rule in the text asked for 2.5 × the largest spot; over 2 mm = 2.3 z_R the pulse was not guided (spot 15 → 22 µm) and the laser field at the absorber rose to 0.16 (11 % of the peak) at the end | `plasma_numbers.py --length` prints both boxes and says when the larger one is needed; `plot_run.py` prints the field at the absorber relative to the peak for every output step; rule 5, the recipe and the debugging table say that a matched spot is no guarantee of guiding |
| default duration cτ = λ_p/2 (56 fs at 10¹⁸ cm⁻³) was long, and "τ" was not defined | default c τ_FWHM = λ_p/3 (intensity FWHM), one definition of `length_1e` everywhere |
| the script's window (not a multiple of dζ) and laser position (rounded) gave two warnings when copied | the script prints `grid.nzeta` and a position with a margin |
| the cost estimate from cell counts was off by a factor 5 | the cost comes from the time per step of the smoke run |
| the quick-look script needed matplotlib | it prints a table with numpy alone; the figure is optional |
| no example without witness; focus position, ramp length, `trap_gamma` value left to guess | `examples/lwfa_plain.deck`; the script prints ramp, focus, `trap_gamma`, `smooth_length` |
| paths `scripts/…`, `examples/…` taken as relative to the working directory | stated: relative to the skill directory |
| wake estimate did not cover a0 ≈ 2; density spikes of 10–40 on the axis not mentioned; no convergence step in the checklist | estimate (a0²/2)/√(1 + a0²/2) E0; note on the spikes; step 6 "Convergence" |

Limits of this trial: the agents ran inside the session that wrote the skill (as sub-agents with the instruction
to ignore everything else), not in a separately installed agent, and the `skills` installer itself was not run
here (no Node.js on the machines used); the layout follows its documented discovery rule (`skills/<name>/SKILL.md`).

## The same for QUARZ and other codes

The pattern carries over as it is: `skills/<code>-setup` and `skills/<code>-simulation`, a copy of the input
reference kept in sync by a script, example inputs that are copies of validated examples, and a test that runs
every input of the skill with the real executable and checks every named parameter against the parser.
