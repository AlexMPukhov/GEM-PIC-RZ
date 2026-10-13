// ============================================================================
//  pipeline.hpp -- MPI parallelization: domain decomposition along zeta,
//  pipelined in the propagation steps (pipeline.cpp)
// ============================================================================
#pragma once
#include <string>

namespace gem {

class Simulation;

struct RunControl {
  std::string deckPath;  // for the restart hint at the wall-time limit
  bool restarted = false;
};

// runs the whole simulation on all ranks of MPI_COMM_WORLD; returns the exit code
int runPipeline(Simulation& sim, const RunControl& rc);

}  // namespace gem
