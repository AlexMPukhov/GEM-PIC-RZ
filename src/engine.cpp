// ============================================================================
//  engine.cpp -- the parallel engine: setup, upload and download (kernels in engine_kernels.hpp,
//  the time step in engine_step.cpp)
// ============================================================================
#include "engine_state.hpp"

namespace gem {

bool Engine::gpuBuild() {
#if defined(GEM_BACKEND_CUDA)
  return true;
#else
  return false;
#endif
}

std::string Engine::describe() { return par::describe(); }

Engine::~Engine() { st->prof.report(); }

Engine::Engine(Simulation& s) : sim(s), st((par::initDevice(), new EngineState)) {
  EngineState& E = *st;
  if (const char* e = std::getenv("GEMPIC_PROFILE")) E.prof.on = (e[0] == '1');
  const RadialGrid& g = sim.g;
  const int M = sim.M, Nz = sim.Nz, N = g.N;
  E.nf = size_t(M + 1) * (Nz + 1) * (N + 1);
  E.nc = size_t(M + 1) * (Nz + 1) * N;
  // grid
  E.rho.assign(g.rho.data(), g.rho.size());
  E.c.assign(g.c.data(), g.c.size());
  E.Af.assign(g.Af.data(), g.Af.size());
  E.Ac.assign(g.Ac.data(), g.Ac.size());
  E.Vf.assign(g.Vf.data(), g.Vf.size());
  E.Vc.assign(g.Vc.data(), g.Vc.size());
  E.kC.assign(g.kC.data(), g.kC.size());
  E.kF.assign(g.kF.data(), g.kF.size());
  // ions and ADK tables
  const auto& ionsCfg = sim.cfg.ions;
  std::vector<dev::IonV> iv;
  std::vector<dev::LevelV> lv;
  std::vector<size_t> levOff;
  for (size_t s = 0; s < ionsCfg.size(); ++s) {
    const auto& I = ionsCfg[s];
    iv.push_back({I.mobile ? 1 : 0, I.ionization ? 1 : 0, I.Zmax, I.mass});
    levOff.push_back(lv.size());
    const AdkRates& R = sim.adk[s];
    for (const auto& L : R.lev) lv.push_back({L.logpre, L.power, L.expo, L.Epeak, L.Emin});
    E.zStride = std::max(E.zStride, int(sim.exitStates[s].size()));
    if (I.ionization) E.maxSpan = std::max(E.maxSpan, I.Zmax);
  }
  E.levels.assign(lv.data(), lv.size());
  std::vector<dev::AdkV> av;
  for (size_t s = 0; s < ionsCfg.size(); ++s) {
    const AdkRates& R = sim.adk[s];
    av.push_back({R.Z0, R.Zmax, R.Emin_all, E.levels.data() + levOff[s]});
  }
  E.ions.assign(iv.data(), iv.size());
  E.adk.assign(av.data(), av.size());
  // tridiagonal matrices (factorized on the host, per matrix class and mode)
  E.uploadMatrices(sim.solver, M);
  {
    const ZetaGrid& zg = sim.zg;
    E.zz.assign(zg.z.data(), zg.z.size());
    E.zhs.assign(zg.hs.data(), zg.hs.size());
    E.zzc.assign(zg.zc.data(), zg.zc.size());
    E.zwlo.assign(zg.wlo.data(), zg.wlo.size());
    E.zwhi.assign(zg.whi.data(), zg.whi.size());
    std::vector<int> cl(Nz);
    for (int n = 0; n < Nz; ++n) cl[n] = sim.solver.classOf(n);
    E.cls.assign(cl.data(), cl.size());
  }
  // parameters
  dev::Params& P = E.P;
  P.M = M; P.Nz = Nz; P.N = N; P.Delta = sim.Delta;
  P.trapExtra = (sim.cfg.numerics.trap_advance - 1.0) * sim.Delta;
  P.theta = sim.cfg.numerics.theta;
  P.zv = ZetaView{Nz, sim.zg.uniform ? 1 : 0, sim.zg.h, E.zz.data(), E.zhs.data(), E.zzc.data(), E.zwlo.data(),
                  E.zwhi.data()};
  P.cls = E.cls.data();
  P.vay = sim.cfg.numerics.beam_vay ? 1 : 0;
  P.beamDzMax = sim.cfg.numerics.beam_dz_max;
  if (!sim.undV.empty()) E.und.assign(sim.undV.data(), sim.undV.size());
  if (!sim.undK.empty()) E.undK.assign(sim.undK.data(), sim.undK.size());
  if (!sim.magV.empty()) E.mag.assign(sim.magV.data(), sim.magV.size());
  P.ext = ExternalV{int(sim.undV.size()), sim.undV.empty() ? nullptr : E.und.data(),
                    sim.undK.empty() ? nullptr : E.undK.data(), int(sim.magV.size()),
                    sim.magV.empty() ? nullptr : E.mag.data()};
  P.g = {N, E.rho.data(), E.c.data(), E.Af.data(), E.Ac.data(), E.Vf.data(), E.Vc.data(), E.kC.data(), E.kF.data()};
  P.nIons = int(ionsCfg.size());
  P.ions = E.ions.data();
  P.adk = E.adk.data();
  P.vz_min = sim.cfg.numerics.vz_min;
  P.vz_beam = sim.cfg.numerics.vz_beam;
  P.seed = sim.cfg.numerics.seed;
  // fields: fold from the host (lasers, bunch fields), fnew = 0
  for (auto& F : E.F) F.init(E.nf, E.nc);
  auto up = [](par::Array<cplx>& d, const std::vector<cd>& h) {
    par::upload(d.data(), reinterpret_cast<const cplx*>(h.data()), h.size());
  };
  dev::FieldStore& Fo = E.F[E.iold];
  up(Fo.Rr, sim.fold.Rr); up(Fo.Lr, sim.fold.Lr); up(Fo.Bz, sim.fold.Bz);
  up(Fo.Rp, sim.fold.Rp); up(Fo.Lp, sim.fold.Lp); up(Fo.Ez, sim.fold.Ez);
  if (sim.cfg.flying_focus.present) {
    E.ffGf.reserve(sim.ffGf.size()); up(E.ffGf, sim.ffGf);
    E.ffGc.reserve(sim.ffGc.size()); up(E.ffGc, sim.ffGc);
  }
  const size_t cf = size_t(M + 1) * Nz * (N + 1), cc = size_t(M + 1) * Nz * N;
  E.nd = size_t(Nz) * N;
  for (int k = 0; k < 2; ++k) { E.Cb[k].init(cf, cc, E.nd); E.Cs[k].init(cf, cc, E.nd); }
  up(E.Cb[E.iold].jr, sim.c0b.jr); up(E.Cb[E.iold].jp, sim.c0b.jp); up(E.Cb[E.iold].jz, sim.c0b.jz);
  par::upload(E.Cb[E.iold].rho.data(), sim.c0b.rho.data(), E.nd);
  par::upload(E.Cb[E.iold].rhoi.data(), sim.c0b.rhoi.data(), E.nd);
  // streaming currents of s^0 (zero at the start, from the checkpoint after a restart)
  up(E.Cs[E.iold].jr, sim.c0s.jr); up(E.Cs[E.iold].jp, sim.c0s.jp); up(E.Cs[E.iold].jz, sim.c0s.jz);
  par::upload(E.Cs[E.iold].rho.data(), sim.c0s.rho.data(), E.nd);
  par::upload(E.Cs[E.iold].rhoi.data(), sim.c0s.rhoi.data(), E.nd);
  for (int k = 0; k < 2; ++k) { E.pred[k].clear(size_t(M + 1) * (N + 1)); E.predP[k].clear(size_t(M + 1) * N); }
  E.solveScratch.clear(size_t(M + 1) * dev::SOLVE_SCRATCH * (N + 1));
  if (sim.sfilter.on()) {  // factors of the source filter (constant for the run)
    std::vector<cplx> band;
    for (const Banded5& B : sim.sfilter.vectorFactors())
      for (const cd& v : B.A) band.push_back(cplx(v.real(), v.imag()));
    E.filtBand.assign(band.data(), band.size());
    std::vector<double> sv;
    for (int k = 0; k < 3; ++k)
      for (int m = 0; m <= M; ++m) {
        const auto& v = sim.sfilter.scalarFactors(k)[m];
        sv.insert(sv.end(), v.begin(), v.end());
      }
    E.filtS.assign(sv.data(), sv.size());
    const size_t ns = size_t(M + 1) * N;
    E.filt = dev::FiltV{1, N, E.filtBand.data(), E.filtS.data(), E.filtS.data() + ns, E.filtS.data() + 2 * ns};
    E.filtScratch.clear(size_t(M + 1) * (4 * size_t(N) + 1));
  }
#if !defined(GEM_BACKEND_CUDA) && !defined(GEM_EMULATE_GPU)
#if defined(_OPENMP)
  E.local = std::vector<dev::CurrStore>(size_t(omp_get_max_threads()));
#else
  E.local = std::vector<dev::CurrStore>(1);
#endif
  for (auto& L : E.local) L.init(size_t(M + 1) * (N + 1), size_t(M + 1) * N, size_t(N));
#endif
  // particles
  E.nBeam = int(sim.beam.size());
  E.beam.reserve(std::max<size_t>(1, sim.beam.size()));
  par::upload(E.beam.data(), sim.beam.data(), sim.beam.size());
  E.cnt.clear(dev::C_COUNT);
  E.lcnt.clear(dev::L_COUNT);
  E.bound.clear(1);
  E.freed.clear(std::max<size_t>(1, ionsCfg.size()));
  E.exitSt.clear(std::max<size_t>(1, ionsCfg.size() * E.zStride));
  E.hcnt[dev::C_BEAM] = E.nBeam;
  E.pushCounters();
  // particles that continue streaming in the next step (after a restart)
  E.nPprev = int(sim.pendingNext.size());
  E.pprev.reserve(std::max<size_t>(1, sim.pendingNext.size()));
  par::upload(E.pprev.data(), sim.pendingNext.data(), sim.pendingNext.size());
  const unsigned long long lc[dev::L_COUNT] = {(unsigned long long)sim.n_restreamed,
                                               (unsigned long long)sim.n_trapped_total};
  par::upload(E.lcnt.data(), lc, dev::L_COUNT);
  par::fence();
}

// ----------------------------------------------------------------------------
void Engine::downloadFull() {
  download();
  EngineState& E = *st;
  auto dn = [](std::vector<cd>& h, par::Array<cplx>& d) {
    par::download(reinterpret_cast<cplx*>(h.data()), d.data(), h.size());
  };
  dn(sim.c0b.jr, E.Cb[E.iold].jr); dn(sim.c0b.jp, E.Cb[E.iold].jp); dn(sim.c0b.jz, E.Cb[E.iold].jz);
  dn(sim.c0s.jr, E.Cs[E.iold].jr); dn(sim.c0s.jp, E.Cs[E.iold].jp); dn(sim.c0s.jz, E.Cs[E.iold].jz);
  sim.pendingNext.resize(size_t(E.nPprev));
  par::download(sim.pendingNext.data(), E.pprev.data(), sim.pendingNext.size());
}

// ----------------------------------------------------------------------------
void Engine::download() {
  EngineState& E = *st;
  E.prof.start();
  par::fence();
  dev::FieldStore& Fo = E.F[E.iold];
  auto dn = [](std::vector<cd>& h, par::Array<cplx>& d) {
    par::download(reinterpret_cast<cplx*>(h.data()), d.data(), h.size());
  };
  dn(sim.fold.Rr, Fo.Rr); dn(sim.fold.Lr, Fo.Lr); dn(sim.fold.Bz, Fo.Bz);
  dn(sim.fold.Rp, Fo.Rp); dn(sim.fold.Lp, Fo.Lp); dn(sim.fold.Ez, Fo.Ez);
  par::download(sim.c0b.rho.data(), E.Cb[E.iold].rho.data(), E.nd);
  par::download(sim.c0b.rhoi.data(), E.Cb[E.iold].rhoi.data(), E.nd);
  par::download(sim.c0s.rho.data(), E.Cs[E.iold].rho.data(), E.nd);
  par::download(sim.c0s.rhoi.data(), E.Cs[E.iold].rhoi.data(), E.nd);
  sim.beam.resize(size_t(E.nBeam));
  par::download(sim.beam.data(), E.beam.data(), sim.beam.size());
  std::vector<double> ex(E.exitSt.capacity());
  par::download(ex.data(), E.exitSt.data(), ex.size());
  for (size_t s = 0; s < sim.exitStates.size(); ++s)
    for (size_t z = 0; z < sim.exitStates[s].size(); ++z) sim.exitStates[s][z] = ex[s * E.zStride + z];
  if (!sim.freedCharge.empty()) par::download(sim.freedCharge.data(), E.freed.data(), sim.freedCharge.size());
  unsigned long long lc[dev::L_COUNT];
  par::download(lc, E.lcnt.data(), dev::L_COUNT);
  sim.n_restreamed = long(lc[dev::L_RESTREAM]);
  sim.n_trapped_total = long(lc[dev::L_TRAPPED]);
  E.prof.stop(Profile::DOWNLOAD);
}

}  // namespace gem
