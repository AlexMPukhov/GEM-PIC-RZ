#!/bin/sh
# Energy spread dE/E in bunches.txt against the particle data of the same run (tests/energy_spread.py): the
# 4 kA, gamma = 1000 beam of benchmarks/harmonics/fel.deck in the planar undulator, 6 steps, dE/E 0 ... 1e-5.
# dE/E is the relative rms of gamma, here ~1e-5: from sum gamma^2 - (sum gamma)^2 / Q the variance cancels to
# ~1e-11 of the sums (it was 27 % low at step 2 and 0 at step 1); the code combines centred moments.
# Requires the printed value to 2e-4 (4 digits) and 0 for the monoenergetic start.
# usage: sh tests/energy_spread.sh <gempic_rz>     (exit 77: skipped, needs python3 with numpy and h5py)
set -e
abspath() { (cd "$(dirname "$1")" && echo "$(pwd)/$(basename "$1")"); }
EXE=$(abspath "$1"); shift
SRC=$(cd "$(dirname "$0")/.." && pwd)
python3 -c "import numpy, h5py" 2>/dev/null || { echo "energy-spread test skipped: needs python3 with numpy and h5py"; exit 77; }
WORK=$(mktemp -d)
cd "$WORK"
"$EXE" "$SRC/benchmarks/harmonics/fel.deck" numerics.backend=serial grid.z_end=3.75mm grid.dzeta=1.328nm \
  bunch.macroparticles=24000 output.format=hdf5 output.particles=true output.every=1 output.dir=out "$@" > run.log
python3 "$SRC/tests/energy_spread.py" out b
cd - > /dev/null
rm -rf "$WORK"
echo "energy-spread test passed"
