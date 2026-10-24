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
#include <complex>
#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "particles.hpp"

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

  // ---- MPI pipeline (pipeline.cpp): the zeta march of the block of a rank on the device. The host keeps the
  // beam list, the exchanges and the output; the state is moved between sim and the device by the calls below
  // (levels: 0 = sim.fold / c0b / c0s, 1 = sim.fnew / c1b / c1s; slice and half-cell ranges [lo, hi) of the
  // stored range). The level index of the device follows swapLevels().
  enum Curr { C0B, C0S, C1B, C1S };
  void putFields(int level, int lo, int hi);
  void getFields(int level, int lo, int hi);
  void putCurrents(Curr which, int lo, int hi);
  void getCurrents(Curr which, int lo, int hi);
  void swapLevels();
  struct BlockMarch {
    int a = 0, b = 0;  // own slices: the march goes from b-1 down to a
    double s1 = 0.0;   // level of the step
    // prediction of j^1 at slice a (filtered beam current of half-cell a-1 + c0s), when a rank is behind
    const std::vector<std::complex<double>>* predBelowR = nullptr;
    const std::vector<std::complex<double>>* predBelowP = nullptr;
    // the front rank's prediction at slice b (the R correction of slice b), when a rank is in front
    const std::vector<std::complex<double>>* predFrontR = nullptr;
    const std::vector<std::complex<double>>* predFrontP = nullptr;
    bool absorbFirst = false;  // absorb slice a after the march (rank 0: slice 0)
    int iF = 0, iC = 0;        // pre-absorber values of faces >= iF, centres >= iC are kept (beam_field = mean)
    int callAfter = -1;        // after this many slices: call after() (the R message to the front rank)
    std::function<void()> after;
  };
  // the march; streaming particles: sim.stream on entry (statuses as in Simulation::march), the survivors in
  // sim.stream on return; particles trapped as beam particles: trapped() (in the order of the march; their
  // beam current is deposited into c1b), slower ones in sim.pendingNext; sim.predRprev / predPprev: the
  // prediction made at slice a; sim.freedCharge
  void marchBlock(const BlockMarch& bm);
  std::vector<Particle> trapped(size_t from = 0);  // particles trapped so far in this march (index >= from)
  // values of the absorber layer before the radial absorber, per absorbed slice (layout of PipelineRun::preAbs)
  bool preAbsorber(int n, std::vector<std::complex<double>>& out);

  // internal (public only because nvcc does not allow device lambdas in private member functions)
  struct MarchSpec;
  void march(const MarchSpec& ms);  // the zeta march over slices ms.nHi ... ms.nLo (engine.cpp)
  void absorbSlices(int n0, int n1, double s1, bool keep);  // radial absorber of fnew on [n0, n1)

 private:
  Simulation& sim;
  std::unique_ptr<EngineState> st;
};

}  // namespace gem
