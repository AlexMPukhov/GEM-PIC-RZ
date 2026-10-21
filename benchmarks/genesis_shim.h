// Force-included when building GENESIS 4 against serial HDF5 (see docs/fel_handoff.md): MPI-IO calls become no-ops.
// Adjust the mpi.h path to your OpenMPI installation.
#pragma once
#if defined(__cplusplus)
#ifndef OMPI_SKIP_MPICXX
#define OMPI_SKIP_MPICXX 1
#endif
#include "/usr/lib/x86_64-linux-gnu/openmpi/include/mpi.h"
#define H5Pset_fapl_mpio(a, b, c) (0)
#define H5Pset_dxpl_mpio(a, b) (0)
#endif
