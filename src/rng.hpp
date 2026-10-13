// ============================================================================
//  rng.hpp -- counter-based random numbers of the field ionization
//
//  A random number is a hash of its key, so it does not depend on the order in
//  which particles are processed (threads, GPU). The key of draw j of a
//  particle in pass `it` of slice n of step `step` is
//      ionKeyParticle(ionKeySlice(ionKeyStep(seed, step), n, it), id) + j
//  All backends use the same keys: serial, openmp and cuda draw identical
//  numbers for the same particle.
// ============================================================================
#pragma once
#include <cstdint>

#include "hd.hpp"

namespace gem {

// SplitMix64 finalizer: a bijective 64-bit mixing function
GEM_HD inline uint64_t mix64(uint64_t z) {
  z += 0x9e3779b97f4a7c15ULL;
  z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
  z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
  return z ^ (z >> 31);
}
// uniform in [0, 1)
GEM_HD inline double uniform01(uint64_t key) { return double(mix64(key) >> 11) * (1.0 / 9007199254740992.0); }

GEM_HD inline uint64_t ionKeyStep(uint64_t seed, int step) {
  return mix64(seed * 0x9E3779B97F4A7C15ULL + uint64_t(step));
}
GEM_HD inline uint64_t ionKeySlice(uint64_t keyStep, int n, int it) {
  return mix64(keyStep ^ ((uint64_t(n) << 8) | uint64_t(it)));
}
GEM_HD inline uint64_t ionKeyParticle(uint64_t keySlice, uint64_t id) { return mix64(keySlice ^ mix64(id)); }

}  // namespace gem
