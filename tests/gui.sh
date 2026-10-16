#!/bin/sh
# GUI (tools/gempic_gui.py) through its HTTP API, see tests/gui_check.py. Exit code 77 (skipped)
# without python3 / numpy / h5py.
# usage: sh tests/gui.sh <gempic_rz>
command -v python3 > /dev/null 2>&1 || { echo "SKIP: python3 not found"; exit 77; }
exec python3 "$(dirname "$0")/gui_check.py" "$1"
