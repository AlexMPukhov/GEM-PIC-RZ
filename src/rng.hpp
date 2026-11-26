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
#include <cmath>
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

// ----------------------------------------------------------------------------
// Quiet ionization (ion block key `split = K`): the K sibling ions of a loading
// slot share the key `id`; sibling j ionizes from charge state Z when its
// accumulated optical depth int w dt reaches the threshold -ln(1 - u) with
//     u = (pi_Z(j) + xi_Z) / K,
// pi_Z a random permutation of 0..K-1 and xi_Z uniform in [0, 1), both drawn per
// (slot, Z). The K thresholds of a level are stratified (one in each 1/K of the
// distribution), so the births of a slot follow the mean ionization history
// instead of K independent random times; for each sibling u is uniform and
// independent of its earlier levels, so the mean rates are unchanged.
// The permutation: rank of the hash of j among the K hashes of the slot.
// ----------------------------------------------------------------------------
GEM_HD inline double quietDepth(uint64_t id, int j, int Z, int K) {
  const uint64_t key = mix64(mix64(id) ^ (0x51ED2701ULL + uint64_t(Z)));
  const uint64_t hj = mix64(key + 1 + uint64_t(j));
  int rank = 0;
  for (int i = 0; i < K; ++i) {
    const uint64_t hi = mix64(key + 1 + uint64_t(i));
    rank += (hi < hj || (hi == hj && i < j)) ? 1 : 0;
  }
  const double u = (rank + uniform01(key)) / K;
  return -log(1.0 - u);
}

}  // namespace gem
