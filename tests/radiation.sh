#!/bin/sh
# Radiation reaction of the beam list (numerics.radiation_reaction; tests/radiation.deck, tests/radiation.py):
# 100 GeV electrons in a helical wiggler, classical / stochastic (zero-inflated Gamma) / Gaussian recoil against an
# RK4 of the Landau-Lifshitz energy and the Sokolov-Ternov variance (on the given backend); the deterministic
# classical run against serial (1e-9). The stochastic runs of two backends agree only statistically: the random
# numbers hash the particle state, which the engine and the serial code reproduce only to rounding.
# Then (tests/radiation_extra.py) a planar undulator with K = 1 (variance of the radiation integral over all
# harmonics, 1.29 x the local constant field) and a thin chicane (one draw of Lambda = 2.1 photons: exp(-Lambda)).
# With an mpiexec as 4th argument all runs are made on 3 ranks (MPI pipeline).
# usage: sh tests/radiation.sh <gempic_rz> <compare_columns> <backend> [mpiexec]   (exit 77: skipped, needs numpy and h5py)
set -e
abspath() { (cd "$(dirname "$1")" && echo "$(pwd)/$(basename "$1")"); }
EXE=$(abspath "$1"); CMP=$(abspath "$2"); B=${3:-serial}; MPIEXEC=${4:-}
RUN="$EXE"
if [ -n "$MPIEXEC" ]; then
  [ "$(id -u)" = 0 ] && export OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
  F=""; if $MPIEXEC --version 2>&1 | grep -qiE "open mpi|openrte|open-mpi"; then F="--oversubscribe"; fi
  RUN="$MPIEXEC $F -n 3 $EXE"
fi
SRC=$(cd "$(dirname "$0")/.." && pwd)
python3 -c "import numpy, h5py" 2>/dev/null || { echo "radiation test skipped: needs python3 with numpy and h5py"; exit 77; }
WORK=$(mktemp -d)
cd "$WORK"
for m in classical stochastic gaussian; do
  $RUN "$SRC/tests/radiation.deck" numerics.backend=$B numerics.radiation_reaction=$m output.dir=out_$m > run_$m.log
done
python3 "$SRC/tests/radiation.py" out_classical out_stochastic out_gaussian 2e5 20.6 0.01 0.07 0.8e-6
# planar undulator, K = 1: the variance of the radiation integral (1.29 x the local constant field)
U="undulator.type=planar undulator.period=1mm undulator.K=1 undulator.periods=100 undulator.natural_focusing=false"
U="$U grid.dz=62.5um grid.z_end=10cm numerics.beam_dz_max=62.5um bunch.macroparticles=40000"
for m in classical stochastic; do
  $RUN "$SRC/tests/radiation.deck" $U numerics.backend=$B numerics.radiation_reaction=$m output.dir=und_$m > und_$m.log
done
python3 "$SRC/tests/radiation_extra.py" undulator und_classical und_stochastic 2e5 1 1e-3 100 1
# thin chicane: Lambda = 2.1 photons in one draw
for m in classical stochastic gaussian; do
  $RUN "$SRC/tests/radiation_chicane.deck" numerics.backend=$B numerics.radiation_reaction=$m output.dir=chic_$m > chic_$m.log
done
python3 "$SRC/tests/radiation_extra.py" chicane chic_classical chic_stochastic chic_gaussian 2e4 2.5e-3 0.2 4
if [ "$B" != serial ] || [ -n "$MPIEXEC" ]; then
  "$EXE" "$SRC/tests/radiation.deck" numerics.backend=serial numerics.radiation_reaction=classical \
    output.dir=out_serial > run_serial.log
  echo " ${MPIEXEC:+3 ranks, }$B vs serial (classical; expected: equal to 1e-9):"
  "$CMP" out_serial/bunches.txt out_classical/bunches.txt 1e-9 1e-12
fi
cd - > /dev/null
rm -rf "$WORK"
echo "radiation test passed"
