// ============================================================================
//  engine.hpp -- parallel time step of GEM-PIC-RZ (GPU or CPU threads)
//
//  The Engine keeps the complete state of a Simulation (fields, currents,
//  particles) in device memory and performs Simulation::step() with kernels
//  (engine.cpp). The algorithm is the one of simulation.cpp; the parallel
//  structure per step is
//
//    beam push            one thread per beam particle
//    zeta march, per slice n and predictor-corrector pass:
//      prepare            one team: predicted currents, predictor fields,
//                         bound of |E| for the ionization
//      push + deposit     one thread per streaming particle (atomic deposition,
//                         freed electrons appended with an atomic counter)
//      field solve        one team per azimuthal mode (operators in parallel
//                         over r, tridiagonal solves by one thread)
//      trapped / freed    one thread per streaming particle
//
//  Nothing is copied between host and device inside the march; after every
//  step download() provides the host copies used by the diagnostics.
//  Compiled by nvcc the engine runs on an NVIDIA GPU (backend "cuda"),
//  compiled by a C++ compiler on CPU threads with OpenMP (backend "openmp").
// ============================================================================
#pragma once
#include <memory>
#include <string>

namespace gem {

class Simulation;
struct EngineState;

class Engine {
 public:
  explicit Engine(Simulation& sim);  // uploads the state of sim (initial or from a checkpoint)
  ~Engine();
  void step();      // one step s^0 -> s^1 (updates sim.s0, sim.step_no)
  void download();  // device -> sim: fields, densities, beam list, counters
  void downloadFull();  // in addition all currents and the pending particles (checkpoints)

  static bool gpuBuild();          // compiled for CUDA
  static std::string describe();   // device / thread information

 private:
  Simulation& sim;
  std::unique_ptr<EngineState> st;
};

}  // namespace gem
