// ============================================================================
//  output_checkpoint.cpp -- checkpoints and restarts (HDF5)
// ============================================================================
#include "output_hdf5_internal.hpp"

#if defined(GEMPIC_HAVE_HDF5)

namespace gem {
using namespace h5detail;

// ============================================================================
//  checkpoints: the complete state, in double precision, so that a restarted
//  run continues exactly as the uninterrupted one
//
//  /                         attributes: kind = "checkpoint", step, z, seed, counters,
//                            grid parameters (checked at the restart)
//  /grid/r_faces, r_centres  radial grid (checked at the restart)
//  /grid/zeta_slices         longitudinal grid (checked at the restart)
//  /state/fields/Rr ... Ez   raw characteristics R, L and E_z, B_z of level s^0
//  /state/beam_currents, /state/stream_currents   j_r, j_phi, j_z, rho, rho_i at s^0
//  /state/beam               beam list (every member of the particle record)
//  /state/pending            particles that continue streaming in the next step
//  /state/rng                state of the random generator of the plasma loading
//  /input_deck               deck of the run that wrote the checkpoint
// ============================================================================
std::string Simulation::writeCheckpoint() {
  namespace fs = std::filesystem;
  if (engine) engine->downloadFull();  // currents and pending particles of the device
  const auto& O = cfg.output;
  fs::create_directories(O.dir);
  char name[1024];
  std::snprintf(name, sizeof name, "%s/checkpoint_%05d.h5", O.dir.c_str(), step_no);
  const std::string path = name, tmp = path + ".tmp";
  H5Eset_auto2(H5E_DEFAULT, nullptr, nullptr);
  {
    Writer W(tmp, /*single=*/false, /*compression=*/0);
    const hid_t root = W.file();
    Writer::attr(root, "code", std::string("GEM-PIC-RZ"));
    Writer::attr(root, "kind", std::string("checkpoint"));
    Writer::attr(root, "format_version", 1LL);
    Writer::attr(root, "step", (long long)step_no);
    Writer::attr(root, "z", s0);
    Writer::attr(root, "seed", (long long)cfg.numerics.seed);
    Writer::attr(root, "n_restreamed", (long long)n_restreamed);
    Writer::attr(root, "n_trapped_total", (long long)n_trapped_total);
    Writer::attr(root, "lambda0_um", cfg.units.lambda0_um);
    Writer::attr(root, "modes", (long long)(M + 1));
    Writer::attr(root, "Nr", (long long)g.N);
    Writer::attr(root, "Nzeta", (long long)Nz);
    Writer::attr(root, "dzeta", cfg.grid.dzeta);
    Writer::attr(root, "dz", cfg.grid.stepAt(s0));  // informative: the z grid may change at a restart
    std::string ions;
    for (size_t s = 0; s < cfg.ions.size(); ++s) ions += (s ? "," : "") + cfg.ions[s].name;
    Writer::attr(root, "ion_names", ions);
    Writer::attr(root, "bunches", (long long)cfg.bunches.size());
    Writer::attr(root, "backend", backendInfo());
    W.group("/grid");
    W.real("/grid/r_faces", {hsize_t(g.N + 1)}, g.rho.data(), 1.0, "1/k0", "faces");
    W.real("/grid/r_centres", {hsize_t(g.N)}, g.c.data(), 1.0, "1/k0", "centres");
    W.real("/grid/zeta_slices", {hsize_t(Nz + 1)}, zg.z.data(), 1.0, "1/k0", "slices");
    W.group("/state");
    W.group("/state/fields");
    const char* desc = "[m][slice][r], level s^0";
    W.complex("/state/fields/Rr", {fold.Rr.size()}, fold.Rr.data(), 1.0, "code units", desc);
    W.complex("/state/fields/Lr", {fold.Lr.size()}, fold.Lr.data(), 1.0, "code units", desc);
    W.complex("/state/fields/Bz", {fold.Bz.size()}, fold.Bz.data(), 1.0, "code units", desc);
    W.complex("/state/fields/Rp", {fold.Rp.size()}, fold.Rp.data(), 1.0, "code units", desc);
    W.complex("/state/fields/Lp", {fold.Lp.size()}, fold.Lp.data(), 1.0, "code units", desc);
    W.complex("/state/fields/Ez", {fold.Ez.size()}, fold.Ez.data(), 1.0, "code units", desc);
    writeCurrents(W, "/state/beam_currents", c0b);
    writeCurrents(W, "/state/stream_currents", c0s);
    writeParticles(W, "/state/beam", beam);
    writeParticles(W, "/state/pending", pendingNext);
    std::ostringstream os;
    os << rng;
    W.text("/state/rng", os.str());
    W.text("/input_deck", deckText);
  }
  fs::rename(tmp, path);  // a crash while writing never leaves a broken checkpoint_NNNNN.h5
  // keep only the newest checkpoint_keep files
  if (O.checkpoint_keep > 0) {
    std::vector<std::pair<int, fs::path>> files;
    for (const auto& e : fs::directory_iterator(O.dir)) {
      const int st = checkpointStep(e.path().filename().string());
      if (st >= 0) files.push_back({st, e.path()});
    }
    std::sort(files.begin(), files.end());
    for (size_t k = 0; k + O.checkpoint_keep < files.size(); ++k) fs::remove(files[k].second);
  }
  return path;
}

void Simulation::readCheckpoint(const std::string& path) {
  H5Eset_auto2(H5E_DEFAULT, nullptr, nullptr);
  Reader R(path);
  auto fail = [&](const std::string& what) {
    throw std::runtime_error("restart from " + path + ": " + what);
  };
  if (!R.has("/state") || R.attrS("kind") != "checkpoint") fail("not a checkpoint file");
  // the grid must be the one of the checkpoint
  auto differs = [](double a, double b) { return std::fabs(a - b) > 1e-12 * std::max(std::fabs(a), std::fabs(b)); };
  if (R.attrI("modes") != M + 1) fail("m_max differs from the checkpoint");
  if (R.attrI("Nr") != g.N) fail("number of radial cells differs from the checkpoint");
  if (R.attrI("Nzeta") != Nz) fail("nzeta differs from the checkpoint");
  if (R.has("/grid/zeta_slices")) {
    const std::vector<double> zs = R.reals("/grid/zeta_slices", zg.z.size());
    for (size_t n = 0; n < zs.size(); ++n)
      if (differs(zs[n], zg.z[n])) fail("the zeta grid differs from the checkpoint");
  } else if (!zg.uniform || differs(R.attrD("dzeta"), cfg.grid.dzeta)) {
    fail("dzeta differs from the checkpoint");
  }
  // (dz and the z grid may change: the scheme takes any step from s^0)
  if (differs(R.attrD("lambda0_um"), cfg.units.lambda0_um)) fail("lambda0 differs from the checkpoint");
  const std::vector<double> rf = R.reals("/grid/r_faces", g.rho.size());
  for (size_t i = 0; i < rf.size(); ++i)
    if (differs(rf[i], g.rho[i])) fail("the radial grid differs from the checkpoint");
  std::string ions;
  for (size_t s = 0; s < cfg.ions.size(); ++s) ions += (s ? "," : "") + cfg.ions[s].name;
  if (R.attrS("ion_names") != ions) fail("the ion blocks differ from the checkpoint ('" + R.attrS("ion_names") + "')");
  if (R.attrI("bunches") != (long long)cfg.bunches.size())
    std::printf("WARNING: the deck has %zu bunch blocks, the checkpoint %lld (bunch names in the output may be wrong)\n",
                cfg.bunches.size(), R.attrI("bunches"));
  if (R.attrI("seed") != (long long)cfg.numerics.seed)
    std::printf("WARNING: numerics.seed differs from the checkpoint: the ionization random numbers differ from the "
                "uninterrupted run\n");
  // state
  R.complexes("/state/fields/Rr", fold.Rr);
  R.complexes("/state/fields/Lr", fold.Lr);
  R.complexes("/state/fields/Bz", fold.Bz);
  R.complexes("/state/fields/Rp", fold.Rp);
  R.complexes("/state/fields/Lp", fold.Lp);
  R.complexes("/state/fields/Ez", fold.Ez);
  readCurrents(R, "/state/beam_currents", c0b);
  readCurrents(R, "/state/stream_currents", c0s);
  beam = readParticles(R, "/state/beam");
  pendingNext = readParticles(R, "/state/pending");
  std::istringstream is(R.text("/state/rng"));
  is >> rng;
  if (!is) fail("cannot restore the random generator");
  step_no = int(R.attrI("step"));
  s0 = R.attrD("z");
  n_restreamed = long(R.attrI("n_restreamed"));
  n_trapped_total = long(R.attrI("n_trapped_total"));
}

std::string Simulation::latestCheckpoint(const std::string& dir) {
  namespace fs = std::filesystem;
  int best = -1;
  std::string path;
  if (fs::is_directory(dir))
    for (const auto& e : fs::directory_iterator(dir)) {
      const int st = checkpointStep(e.path().filename().string());
      if (st > best) { best = st; path = e.path().string(); }
    }
  return path;
}

}  // namespace gem

#endif
