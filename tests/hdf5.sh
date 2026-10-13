#!/bin/sh
# HDF5 output: run a short version of tests/regression.deck with
# output.format = both and check that the HDF5 files, read back with
# tools/h5_extract, reproduce the text output (axis, x-z map, beam list);
# then the same in single precision with compression.
# usage: sh tests/hdf5.sh <gempic_rz> <h5_extract> <compare_columns>
set -e
abspath() { (cd "$(dirname "$1")" && echo "$(pwd)/$(basename "$1")"); }
EXE=$(abspath "$1"); H5X=$(abspath "$2"); CMP=$(abspath "$3")
SRC=$(cd "$(dirname "$0")/.." && pwd)
WORK=$(mktemp -d)
cd "$WORK"
RUN="numerics.backend=serial grid.z_end=40k0^-1 output.every=1 output.format=both"
"$EXE" "$SRC/tests/regression.deck" $RUN output.dir=d > d.log
"$EXE" "$SRC/tests/regression.deck" $RUN output.dir=s output.hdf5_precision=single output.hdf5_compression=6 > s.log
for dir in d s; do
  if [ $dir = d ]; then TOL="1e-6 1e-9"; else TOL="1e-5 1e-6"; fi
  "$H5X" $dir/data_00001.h5 axis > axis.txt
  "$H5X" $dir/data_00001.h5 map 4 > map.txt
  "$H5X" $dir/data_00001.h5 beam > beam.txt
  "$CMP" $dir/axis_00001.txt axis.txt $TOL
  "$CMP" $dir/map_00001.txt map.txt $TOL
  "$CMP" $dir/beam_00001.txt beam.txt $TOL
done
"$H5X" d/data_00001.h5 info
ls -l d/data_00001.h5 s/data_00001.h5 | awk '{print "  " $5 " bytes  " $9}'
cd - > /dev/null
rm -rf "$WORK"
echo "hdf5 test passed"
