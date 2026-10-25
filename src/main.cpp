// ============================================================================
//  main.cpp -- command line driver of GEM-PIC-RZ
//
//  gempic_rz <deck> [block.key=value ...] [--check]
//
//    <deck>              input deck (see docs/input_deck.md)
//    block.key=value     override a deck value, e.g.  grid.z_end=2mm
//                        laser[1].a0=1.5   "bunch.name=witness 2"
//    --check             read and check the deck, print the summary, exit
//    --template          print an annotated example deck and exit
//    --version           print the version and exit
// ============================================================================
#include "version.hpp"
#include <chrono>
#if defined(_OPENMP)
#include <omp.h>
#endif
#include <memory>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#include "config.hpp"
#include "deck.hpp"
#include "comm.hpp"
#include "pipeline.hpp"
#include "simulation.hpp"

namespace {

const char* TEMPLATE = R"(# GEM-PIC-RZ input deck (template)       lengths need units: nm um mm cm m, or k0^-1
grid {
  lambda0          = 0.8 um        # reference wavelength (defines code units)
  m_max            = 1             # azimuthal modes 0..m_max (1: linearly polarized laser)
  radial_grid      = table         # table | stretched | uniform
  r_points         = 0,  1.5, 4,   24 um   # reference radii ...
  dr_points        = 40, 40,  400, 400 nm  # ... and the radial step there (linear in between)
  r_max            = 24 um         # outer radius (or: nr = 110)
  dzeta            = 32 nm         # longitudinal step
  window_length    = 38.4 um       # (or: nzeta = 1200)
  dz               = 5 um          # step along the propagation distance
  z_end            = 0.5 mm
  r_absorber       = 20 um         # radial absorber outside this radius
}

laser {                            # repeat the block for more pulses
  a0           = 2
  waist        = 8 um              # 1/e field radius at focus
  tau_fwhm     = 12.5 fs           # intensity FWHM (or: length_1e = 3.2 um)
  zeta_center  = 28 um             # position in the window (0 = rear end)
  z_focus      = 0 mm
  polarization = linear            # linear | circular
}

plasma {
  density        = 3.5e18 cm^-3
  z_points       = 0, 50, 10000 um # piecewise-linear longitudinal profile
  z_values       = 0, 1, 1
  radial_profile = uniform         # uniform | parabolic | table
  r_max          = 19 um
}

bunch {                            # repeat the block for more bunches
  name          = witness
  charge        = 5 pC
  energy        = 510 MeV
  energy_spread = 0.1 %
  sigma_x       = 0.38 um
  emittance_x   = 0.3 um           # normalized rms emittance (mm mrad)
  sigma_z       = 0.5 um
  zeta_center   = 16.5 um
  macroparticles = 20000
}

output {
  dir    = diag
  every  = 25
  format = text                    # text | hdf5 | both  (hdf5: data_NNNNN.h5 with all modes)
  checkpoint_every = 0             # checkpoint_NNNNN.h5 for  gempic_rz deck --restart latest  (0: none)
}
)";

void usage() {
  std::fprintf(stderr,
               "usage: gempic_rz <deck> [block.key=value ...] [--check] [--restart <checkpoint.h5> | latest]\n"
               "       gempic_rz --template        (print an example deck)\n"
               "       gempic_rz --version\n");
}

}  // namespace

int run(int argc, char** argv) {
  std::string deckPath;
  std::vector<std::string> overrides;
  bool checkOnly = false;
  std::string restart;  // checkpoint file or "latest"
  for (int a = 1; a < argc; ++a) {
    const std::string s(argv[a]);
    if (s == "--check") checkOnly = true;
    else if (s == "--restart" && a + 1 < argc) restart = argv[++a];
    else if (s.rfind("--restart=", 0) == 0) restart = s.substr(10);
    else if (s == "--template") { std::fputs(TEMPLATE, stdout); return 0; }
    else if (s == "--version") { std::printf("GEM-PIC-RZ %s\n", GEMPIC_VERSION); return 0; }
    else if (s == "-h" || s == "--help") { usage(); return 0; }
    else if (s.find('=') != std::string::npos) overrides.push_back(s);
    else if (deckPath.empty()) deckPath = s;
    else { std::fprintf(stderr, "unexpected argument '%s'\n", s.c_str()); usage(); return 1; }
  }
  if (deckPath.empty()) { usage(); return 1; }

  gem::Config cfg;
  gem::Deck deck;
  try {
    deck = gem::Deck::parseFile(deckPath);
    for (const auto& o : overrides) deck.applyOverride(o);
    std::vector<std::string> warnings;
    cfg = gem::Config::fromDeck(deck, &warnings);
    std::printf("GEM-PIC-RZ %s   deck: %s\n%s", GEMPIC_VERSION, deckPath.c_str(), cfg.summary().c_str());
    for (const auto& w : warnings) std::printf("WARNING: %s\n", w.c_str());
  } catch (const gem::DeckError& e) {
    std::fprintf(stderr, "error in input deck: %s\n", e.what());
    return 1;
  }
  if (restart == "latest") {
    const std::string p = gem::Simulation::latestCheckpoint(cfg.output.dir);
    if (p.empty()) {
      std::fprintf(stderr, "error: --restart latest: no checkpoint_NNNNN.h5 in '%s'\n", cfg.output.dir.c_str());
      return 1;
    }
    restart = p;
  }
  if (!restart.empty()) std::printf("restart from %s\n", restart.c_str());
  if (checkOnly) return 0;
  const int nranks = gem::comm::size();
  const char* forcePipe = std::getenv("GEMPIC_PIPELINE");  // pipelined step also on one rank (testing)
  const bool pipeline = nranks > 1 || (forcePipe && forcePipe[0] == '1');

  // keep a record of the run
  if (gem::comm::rank() == 0) {
    std::filesystem::create_directories(cfg.output.dir);
    const std::string sfx = restart.empty() ? "" : "_restart";  // keep the record of the first run
    std::ofstream(cfg.output.dir + "/deck_used" + sfx + ".txt")
        << "# deck " << deckPath << " after command-line overrides" << (restart.empty() ? "" : ", restart from ")
        << restart << "\n\n" << deck.dump();
    std::ofstream(cfg.output.dir + "/summary" + sfx + ".txt") << cfg.summary();
  }

  std::unique_ptr<gem::Simulation> simp;
  try {
    if (pipeline) {  // every rank builds only its own block of the window
      const gem::StorageBlock blk = gem::pipelineBlock(cfg, nranks, gem::comm::rank());
      simp = std::make_unique<gem::Simulation>(cfg, restart, false, &blk);
    } else {
      simp = std::make_unique<gem::Simulation>(cfg, restart, true);
    }
  } catch (const std::exception& e) {
    std::fprintf(stderr, "error: %s\n", e.what());
    return 1;
  }
  gem::Simulation& sim = *simp;
  sim.setDeckText(deck.dump());
  if (pipeline) {
    int nth = 1;
#if defined(_OPENMP)
    if (cfg.numerics.backend == gem::NumericsConfig::Backend::OpenMP) nth = omp_get_max_threads();
#endif
    const bool gpu = cfg.numerics.backend == gem::NumericsConfig::Backend::Cuda;
    std::printf("backend: %s on %d MPI rank%s%s\n", gpu ? "cuda" : nth > 1 ? "openmp" : "serial", nranks,
                nranks > 1 ? "s" : "",
                nth > 1 ? (", " + std::to_string(nth) + " threads per rank").c_str() : "");
    gem::RunControl rc;
    rc.deckPath = deckPath;
    rc.restarted = !restart.empty();
    return gem::runPipeline(sim, rc);
  }
  std::printf("backend: %s\n", sim.backendInfo().c_str());
  const double zend = cfg.grid.z_end;
  const int every = cfg.output.every;
  auto t0 = std::chrono::steady_clock::now();
  const int firstStep = sim.stepNumber();
  if (restart.empty()) sim.diagnostics(true);  // (after a restart: written by the interrupted run)
  else std::printf("continuing at step %d, z = %.6g mm\n", firstStep, sim.z() * cfg.units.k0inv_um() * 1e-3);
  const auto& O = cfg.output;
  auto hours = [&] { return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count() / 3600.0; };
  while (sim.z() < zend - cfg.grid.endTol()) {
    try {
      sim.step();
    } catch (const std::exception& e) {
      std::fprintf(stderr, "error in step %d: %s\n", sim.stepNumber() + 1, e.what());
      return 2;
    }
    if (!sim.finite()) {
      std::fprintf(stderr, "non-finite field at step %d\n", sim.stepNumber());
      return 2;
    }
    const bool last = sim.z() >= zend - cfg.grid.endTol();
    sim.diagnostics((every > 0 && sim.stepNumber() % every == 0) || last);
    // checkpoints; stop before the wall-time limit (1.5 mean steps of margin)
    const int done = sim.stepNumber() - firstStep;
    const bool wall = O.walltime_hours > 0 && !last && hours() + 1.5 * hours() / done > O.walltime_hours;
    if ((O.checkpoint_every > 0 && (sim.stepNumber() % O.checkpoint_every == 0 || last)) || wall) {
      try {
        std::printf("checkpoint: %s\n", sim.writeCheckpoint().c_str());
      } catch (const std::exception& e) {
        std::fprintf(stderr, "error writing the checkpoint: %s\n", e.what());
        return 2;
      }
    }
    if (wall) {
      std::printf("wall-time limit (%.3g h) reached at step %d; continue with:  gempic_rz %s --restart latest\n",
                  O.walltime_hours, sim.stepNumber(), deckPath.c_str());
      return 3;
    }
  }
  const double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  std::printf("done: %d steps in %.1f s, particles returned from the beam list to streaming: %ld\n", sim.stepNumber(),
              sec, sim.restreamedParticles());
  return 0;
}

int main(int argc, char** argv) {
  gem::comm::init(&argc, &argv);
  // only rank 0 writes to the console (errors go to stderr on every rank)
  if (gem::comm::rank() != 0 && !std::freopen("/dev/null", "w", stdout)) return 1;
  const int rc = run(argc, argv);
  gem::comm::finalize();
  return rc;
}
