// ============================================================================
//  config.cpp -- reading the deck blocks into typed configuration structures
// ============================================================================
#include "config_internal.hpp"

namespace gem {
using namespace cfgdetail;

// ============================================================================
//  helper methods of the config structures
// ============================================================================
double GridConfig::mapping(double x) const {
  if (radial == Radial::Table) return rtable.coord(x);
  return dr_min * x + 0.5 * (dr_max - dr_min) * (x + wcells * (logcosh((x - xi0) / wcells) - logcosh(-xi0 / wcells)));
}

double GridConfig::zetaNode(int n) const { return zetaTable ? ztable.coord(n) : n * dzeta; }

double GridConfig::zetaStep(double zeta) const { return zetaTable ? ztable.step(zeta) : dzeta; }

double GridConfig::stepAt(double s) const {
  if (!zTable) return dz;
  return stable.coord(stable.index(s) + 1.0) - s;
}

int GridConfig::stepsFrom(double s) const {
  int k = 0;
  while (s < z_end - endTol()) { s = s + stepAt(s); ++k; }  // the arithmetic of the run loop
  return k;
}

double UndulatorConfig::ku() const { return 2.0 * PI / period; }

double UndulatorConfig::resonance(double gamma) const {
  const double kk = (type == Type::Planar) ? 0.5 * K * K : K * K;
  return period * (1.0 + kk) / (2.0 * gamma * gamma);
}

double PlasmaConfig::longitudinal(double z) const {
  if (z_points.empty()) return 1.0;
  return interpTable(z_points, z_values, z);
}

double PlasmaConfig::radialFactor(double r) const {
  switch (radial) {
    case Radial::Uniform: return 1.0;
    case Radial::Parabolic: return 1.0 + channel_depth * r * r / (r_channel * r_channel);
    case Radial::Table: return interpTable(r_points, r_values, r);
  }
  return 1.0;
}

double PlasmaConfig::density(double r, double z) const {
  if (!present || r >= r_max) return 0.0;
  return n0 * longitudinal(z) * radialFactor(r);
}

double PlasmaConfig::density(double x, double y, double z) const {
  const double r = std::sqrt(x * x + y * y);
  const double n = density(r, z);
  if (n <= 0.0 || fn.empty()) return n;
  const double v[4] = {x, y, z, r};
  return std::max(0.0, n * fn(v));
}

// ============================================================================
//  box size for bunches with v < c
//
//  Radial: the field of a bunch of length L reaches the radius ~ gamma L (rest
//  frame: rod of length gamma L). A wall at r_max < gamma L shields the
//  longitudinal space charge, g = 1 + 2 ln(min(r_max, gamma L)/a) (a beam
//  pipe); r_max >= gamma L gives the free-space field. In an undulator the
//  radiation needs r_max >= 5 w with w^2 = a^2 + (lambda_r Z / (4 pi a))^2 at
//  the end of the undulator (Z: distance from its entrance).
//  (Checked with a gamma = 100 bunch: r_max = gamma L shields the energy
//  spread from the space charge by 1.4 %, 2 gamma L by 0.3 %.)
//  Window: the longitudinal field reaches the distance d beyond the ends of
//  the bunch. In free space E_z(d) = Q/(4 pi eps0 gamma^2 d^2), against
//  E_z(edge) ~ Q/(4 pi eps0 a max(a, gamma L)) (rod or disk in the rest
//  frame); inside the wall it decays as exp(-2.405 gamma d / r_max). 1 %:
//      d = min( 10 sqrt(a max(a, gamma L)) / gamma, ln(100) r_max / (2.405 gamma) ).
//  Inside a plasma the field is screened over c/omega_p: no advice then.
//  The bunch slips back by int (1 + K_eff^2)/(2 gamma^2) dz (K_eff^2 = K^2/2
//  planar, K^2 helical), so the gap behind it is checked at the end of the run.
// ============================================================================
std::vector<Config::BoxAdvice> Config::boxAdvice() const {
  std::vector<BoxAdvice> out;
  if (plasma.present) return out;  // fields screened over c / omega_p
  for (const auto& B : bunches) {
    if (!(B.gamma < 1e5)) continue;
    BoxAdvice A;
    A.name = B.name;
    const double half = (B.shape == BunchConfig::Shape::Gaussian) ? 3.0 * B.sigma_z : 0.5 * B.length;
    const double L = 2.0 * half, a = std::max(B.sigma_x, B.sigma_y), g2 = 2.0 * B.gamma * B.gamma;
    A.rmax_need = B.gamma * L;
    A.rmax_why = "gamma L";
    for (const auto& U : undulators) {
      const double Z = std::max(0.0, std::min(grid.z_end, U.z_start + U.length) - U.z_start);
      const double lr = U.resonance(B.gamma), d = lr * Z / (4.0 * PI * a);
      const double w5 = 5.0 * std::sqrt(a * a + d * d);
      if (w5 > A.rmax_need) { A.rmax_need = w5; A.rmax_why = "5 x the size of the undulator radiation"; }
    }
    A.rmax_now = grid.r_max;
    A.gap_need = std::min(10.0 * std::sqrt(a * std::max(a, B.gamma * L)) / B.gamma,
                          std::log(100.0) * grid.r_max / (2.405 * B.gamma));
    // slippage over the run
    double slip = grid.z_end / g2;
    for (const auto& U : undulators) {
      const double Z = std::max(0.0, std::min(grid.z_end, U.z_start + U.length) - U.z_start);
      double k2 = 0.0;  // mean K^2 over the part inside the run (taper: K(z))
      for (int j = 0; j < 200; ++j) {
        const double Kz = U.Kat(U.z_start + (j + 0.5) / 200.0 * Z);
        k2 += Kz * Kz / 200.0;
      }
      slip += Z * ((U.type == UndulatorConfig::Type::Planar) ? 0.5 * k2 : k2) / g2;
    }
    A.slip = slip;
    A.front_now = grid.window() - (B.zeta_center + half);
    A.rear_end = B.zeta_center - half - slip;
    // smallest window that keeps both gaps
    A.window_need = A.gap_need + slip + L + A.gap_need;
    out.push_back(A);
  }
  return out;
}

// ============================================================================
Config Config::fromDeck(const Deck& deck, std::vector<std::string>* warnings) {
  std::vector<std::string> warn;
  deck.checkBlockNames({"constants", "grid", "laser", "flying_focus", "plasma", "ion", "bunch", "undulator", "quadrupole",
                        "plasma_lens", "solenoid", "chicane", "output", "numerics"});
  Config c;
  g_consts = &c.constants;  // the expressions of all blocks see the constants
  struct Reset { ~Reset() { g_consts = nullptr; } } reset;
  const Block* grid = deck.single("grid");
  if (!grid) throw DeckError(deck.source + ": the deck needs a 'grid' block");
  c.units.lambda0_um = readLambda0(deck, grid);

  if (const Block* b = deck.single("constants")) readConstants(*b, deck, c.units, c.constants);
  readGrid(*grid, deck, c.units, c.grid, warn);
  for (const Block* b : deck.all("laser")) {
    LaserConfig L;
    readLaser(*b, deck, c.units, c.grid, L, warn);
    c.lasers.push_back(L);
  }
  if (!c.lasers.empty() && c.grid.M < 1)
    throw DeckError(deck.where(grid->line) + ": a laser needs azimuthal mode 1: set grid.m_max >= 1");
  if (const Block* b = deck.single("flying_focus")) readFlyingFocus(*b, deck, c.units, c.grid, c.flying_focus, warn);
  if (const Block* b = deck.single("plasma")) readPlasma(*b, deck, c.units, c.grid, c.plasma, warn);
  double fsum = 0.0;
  for (const Block* b : deck.all("ion")) {
    if (!c.plasma.present)
      throw DeckError(deck.where(b->line) + ": an 'ion' block needs a 'plasma' block (ions follow its profile)");
    if (c.ions.size() >= 10) throw DeckError(deck.where(b->line) + ": at most 10 ion blocks are supported");
    IonConfig I;
    readIon(*b, deck, c.units, c.plasma, I, warn);
    fsum += I.fraction();
    c.ions.push_back(I);
  }
  if (!c.plasma.free_electrons && fsum > 0.0) {
    throw DeckError(deck.where(deck.single("plasma")->line) +
                    ": plasma.free_electrons = false, but ion blocks with charge_state > 0 carry positive charge "
                    "(it must be neutralized by free electrons)");
  }
  if (fsum > 1.0 + 1e-9) {
    const auto v = deck.all("ion");
    throw DeckError(deck.where(v.back()->line) + ": the ion fractions add up to " + std::to_string(fsum) +
                    " > 1 (fraction = charge_state x atom density / free-electron density)");
  }
  if (c.plasma.present && !c.plasma.free_electrons && c.ions.empty())
    warn.push_back("plasma.free_electrons = false and no ion blocks: the plasma block has no effect");
  for (const Block* b : deck.all("bunch")) {
    BunchConfig B;
    readBunch(*b, deck, c.units, c.grid, B, warn);
    c.bunches.push_back(B);
  }
  for (const Block* b : deck.all("undulator")) {
    UndulatorConfig U;
    readUndulator(*b, deck, c.units, U);
    // planar undulator, high current (benchmarks/harmonics: 4 kA, gamma 1000, K = 1.5): with the field step
    // lambda_u / 16 the beam blew up transversely (eps_n 0 -> 2.8 um in 25 periods; lambda_u / 24: 0.13 um,
    // / 32: 0.04 um, / 64: 0.02 um), and the gain dropped; the helical case is far less sensitive
    const double dzMax = c.grid.zTable ? c.grid.dz_max : c.grid.dz;
    double qmax = 0.0;
    for (const auto& B : c.bunches) qmax = std::max(qmax, B.charge_pC);
    if (U.type == UndulatorConfig::Type::Planar && qmax >= 1.0 && dzMax > U.period / 32.0 * 1.0001)
      warn.push_back(b->label() + ": planar undulator with dz = period / " +
                     std::to_string(int(std::lround(U.period / dzMax))) +
                     ": use dz <= period / 32 with intense beams (at period / 16 a 4 kA beam blew up transversely)");
    c.undulators.push_back(U);
  }
  for (const char* kind : {"quadrupole", "plasma_lens", "solenoid"})
    for (const Block* b : deck.all(kind)) {
      MagnetConfig G;
      readMagnet(*b, deck, c.units, G, warn);
      c.magnets.push_back(G);
    }
  for (const Block* b : deck.all("chicane")) {
    ChicaneConfig C;
    readChicane(*b, deck, c.units, C);
    c.chicanes.push_back(C);
  }
  std::sort(c.chicanes.begin(), c.chicanes.end(), [](const ChicaneConfig& a, const ChicaneConfig& b) { return a.z < b.z; });
  if (const Block* b = deck.single("output")) readOutput(*b, deck, c.units, c.output);
  // FEL diagnostics: on with undulators unless switched off; wavelength of the bunching
  {
    auto& O = c.output;
    O.fel = (O.fel_mode == 1) || (O.fel_mode == -1 && !c.undulators.empty() && !c.bunches.empty());
    if (O.fel && !(O.fel_wavelength > 0)) {
      if (c.undulators.empty() || c.bunches.empty())
        throw DeckError(deck.source + ": output.fel needs output.fel_wavelength (no undulator and bunch to take the "
                        "resonant wavelength from)");
      O.fel_wavelength = c.undulators[0].resonance(c.bunches[0].gamma);
    }
    if (!(O.fel_slice > 0)) O.fel_slice = 20.0 * O.fel_wavelength;
  }
  if (const Block* b = deck.single("numerics")) readNumerics(*b, deck, c.units, c.numerics);
  deck.checkAllUsed();
  // box size for bunches with v < c (Config::boxAdvice): warn if too small
  for (const auto& A : c.boxAdvice()) {
    if (A.rmax_now < A.rmax_need)
      warn.push_back("bunch '" + A.name + "': the radial boundary is closer than the reach of its fields: r_max = " +
                     fmtUm(A.rmax_now, c.units) + " < " + fmtUm(A.rmax_need, c.units) + " (" + A.rmax_why +
                     "); a wall closer than gamma L shields the space charge like a beam pipe");
    if (A.front_now < A.gap_need || A.rear_end < A.gap_need)
      warn.push_back("bunch '" + A.name + "': its space-charge field reaches the ends of the window (gap ahead " +
                     fmtUm(A.front_now, c.units) + ", behind at the end " + fmtUm(A.rear_end, c.units) +
                     "; needed " + fmtUm(A.gap_need, c.units) + "): E_z = 0 at the front acts as a uniform E_z");
  }
  // resolution of an undulator: the steps along z and the sub-steps of the
  // beam push must follow the wiggling (period / 16), the zeta grid the
  // radiation (resonant wavelength / 12)
  for (size_t k = 0; k < c.undulators.size(); ++k) {
    const auto& U = c.undulators[k];
    const std::string lab = "undulator[" + std::to_string(k) + "]";
    if (c.grid.dz_max > U.period / 16.0)
      warn.push_back(lab + ": grid.dz (up to " + std::to_string(c.grid.dz_max / U.period) +
                     " periods) should be <= period / 16 (the transverse current oscillates with the period)");
    if (c.numerics.beam_dz_max > U.period / 16.0)
      warn.push_back(lab + ": numerics.beam_dz_max should be <= period / 16 (sub-steps of the beam push)");
    for (const auto& B : c.bunches) {
      const double lr = U.resonance(B.gamma);
      if (c.grid.dzeta_max > lr / 12.0)
        warn.push_back(lab + ", bunch '" + B.name + "': the resonant wavelength is only " +
                       std::to_string(lr / c.grid.dzeta_max) + " zeta steps (use >= 12, better ~20)");
    }
  }
  // transport elements: the sub-steps of the beam push must resolve the element and the betatron phase
  for (const auto& G : c.magnets) {
    if (G.length < 4.0 * c.numerics.beam_dz_max)
      warn.push_back(G.label + ": shorter than 4 sub-steps of the beam push (numerics.beam_dz_max)");
    for (const auto& B : c.bunches) {
      const double k = (G.kind == MagnetConfig::Kind::Solenoid) ? 0.5 * std::fabs(G.field) / B.gamma
                                                                  : std::sqrt(std::fabs(G.gradient) / B.gamma);
      if (k * c.numerics.beam_dz_max > 0.1)
        warn.push_back(G.label + ", bunch '" + B.name + "': numerics.beam_dz_max exceeds 0.1 / k_beta (" +
                       fmtUm(0.1 / k, c.units) + ")");
    }
  }

  if (c.lasers.empty() && c.bunches.empty() && !c.flying_focus.present)
    warn.push_back("the deck defines neither a laser nor a bunch");
  if (!c.plasma.present) warn.push_back("no plasma block: vacuum propagation");
  if (warnings) *warnings = warn;
  return c;
}

}  // namespace gem
