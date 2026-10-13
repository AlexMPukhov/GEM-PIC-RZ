// ============================================================================
//  portable.hpp -- minimal portability layer for the parallel engine
//
//  The kernels of engine.cpp are written once, as lambdas, and run either
//    * on an NVIDIA GPU (the file is compiled by nvcc: __CUDACC__), or
//    * on the CPU with OpenMP threads (or one thread without OpenMP).
//
//  Primitives
//    par::alloc / release / upload / download / copy / zero   memory ("device")
//    par::forN(n, f)              f(i) for 0 <= i < n  (n known on the host)
//    par::forDyn(a, b, cap, f)    f(i) for 0 <= i < *a + *b  (counts live in
//                                 device memory; cap = upper bound for the launch).
//                                 ordered = true: sequential on the CPU (keeps
//                                 the order of appends deterministic)
//    par::teams(nt, f)            f(team) for nt teams; a team is a CUDA block
//                                 (TEAM threads, barrier = __syncthreads) or a
//                                 single CPU thread
//    atomics, complex numbers (random numbers: rng.hpp)
//
//  GEM_EMULATE_GPU (CPU build with OpenMP, for testing): the GPU code paths
//  run on the CPU -- teams of GEM_EMULATE_TEAM threads with real barriers,
//  atomic deposition, unordered appends. Races or missing barriers in the
//  team kernels then show up as differences to the serial code.
//
//  Lambdas passed to the launchers are written GEM_LAMBDA(...) { ... } and
//  must capture only local variables (no `this`): they are copied to the GPU.
// ============================================================================
#pragma once
#include <math.h>

#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <stdexcept>
#include <string>
#include <utility>

#include "hd.hpp"

#if defined(__CUDACC__)
#include <cuda_runtime.h>
#define GEM_LAMBDA [=] __host__ __device__
#define GEM_BACKEND_CUDA 1
#else
#define GEM_LAMBDA [=]
#if defined(_OPENMP)
#include <omp.h>
#endif
#endif

namespace gem {

// ---------------------------------------------------------------- complex
// Same memory layout as std::complex<double> (host arrays are shared by cast).
struct cplx {
  double re, im;
  GEM_HD cplx() : re(0.0), im(0.0) {}
  GEM_HD cplx(double r, double i = 0.0) : re(r), im(i) {}
  GEM_HD cplx& operator+=(const cplx& b) { re += b.re; im += b.im; return *this; }
  GEM_HD cplx& operator-=(const cplx& b) { re -= b.re; im -= b.im; return *this; }
};
GEM_HD inline cplx operator+(const cplx& a, const cplx& b) { return cplx(a.re + b.re, a.im + b.im); }
GEM_HD inline cplx operator-(const cplx& a, const cplx& b) { return cplx(a.re - b.re, a.im - b.im); }
GEM_HD inline cplx operator*(const cplx& a, const cplx& b) {
  return cplx(a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re);
}
GEM_HD inline cplx operator*(double s, const cplx& a) { return cplx(s * a.re, s * a.im); }
GEM_HD inline cplx operator*(const cplx& a, double s) { return cplx(a.re * s, a.im * s); }
GEM_HD inline cplx operator/(const cplx& a, double s) { return cplx(a.re / s, a.im / s); }
GEM_HD inline double cabs(const cplx& a) { return hypot(a.re, a.im); }

namespace par {

constexpr int BLOCK = 256;  // threads per block for particle / grid loops
constexpr int TEAM = 128;   // threads per team (power of 2)

// ---------------------------------------------------------------- atomics
GEM_HD inline void atomicAddD(double* a, double v) {
#if defined(__CUDA_ARCH__)
  atomicAdd(a, v);
#elif defined(_OPENMP)
#pragma omp atomic
  *a += v;
#else
  *a += v;
#endif
}
GEM_HD inline void atomicAddC(cplx* a, const cplx& v) {
  atomicAddD(&a->re, v.re);
  atomicAddD(&a->im, v.im);
}
// returns the old value
GEM_HD inline int atomicFetchAdd(int* a, int v) {
#if defined(__CUDA_ARCH__)
  return atomicAdd(a, v);
#else
  int old;
#if defined(_OPENMP)
#pragma omp atomic capture
#endif
  { old = *a; *a += v; }
  return old;
#endif
}
GEM_HD inline void atomicAddU64(unsigned long long* a, unsigned long long v) {
#if defined(__CUDA_ARCH__)
  atomicAdd(a, v);
#elif defined(_OPENMP)
#pragma omp atomic
  *a += v;
#else
  *a += v;
#endif
}

// ---------------------------------------------------------------- teams
struct Team {
  int league, rank, size;
  double* scratch;  // size doubles, shared by the team
  GEM_HD void barrier() const {
#if defined(__CUDA_ARCH__)
    __syncthreads();
#elif defined(GEM_EMULATE_GPU) && defined(_OPENMP)
#pragma omp barrier
#endif
  }
  GEM_HD bool leader() const { return rank == 0; }
};
template <class F>
GEM_HD inline void teamFor(const Team& t, int n, F f) {
  for (int i = t.rank; i < n; i += t.size) f(i);
}
// maximum of v over the team (all members get the result)
GEM_HD inline double teamMax(const Team& t, double v) {
  t.scratch[t.rank] = v;
  t.barrier();
  for (int s = t.size / 2; s > 0; s >>= 1) {
    if (t.rank < s && t.scratch[t.rank + s] > t.scratch[t.rank]) t.scratch[t.rank] = t.scratch[t.rank + s];
    t.barrier();
  }
  const double r = t.scratch[0];
  t.barrier();
  return r;
}

#if defined(GEM_BACKEND_CUDA)
// ============================================================================
//  CUDA backend
// ============================================================================
inline void check(cudaError_t e, const char* what) {
  if (e == cudaErrorNoKernelImageForDevice || e == cudaErrorUnsupportedPtxVersion)
    throw std::runtime_error(std::string("CUDA error in ") + what + ": " + cudaGetErrorString(e) +
                             " -- the binary has no code for this GPU: rebuild with -DCMAKE_CUDA_ARCHITECTURES=native "
                             "(or the compute capability of the GPU, e.g. 80, 86, 89, 90)");
  if (e != cudaSuccess) throw std::runtime_error(std::string("CUDA error in ") + what + ": " + cudaGetErrorString(e));
}
template <class T> T* alloc(size_t n) {
  T* p = nullptr;
  if (n) check(cudaMalloc(reinterpret_cast<void**>(&p), n * sizeof(T)), "cudaMalloc");
  return p;
}
template <class T> void release(T* p) { if (p) cudaFree(p); }
template <class T> void upload(T* d, const T* h, size_t n) {
  if (n) check(cudaMemcpy(d, h, n * sizeof(T), cudaMemcpyHostToDevice), "upload");
}
template <class T> void download(T* h, const T* d, size_t n) {
  if (n) check(cudaMemcpy(h, d, n * sizeof(T), cudaMemcpyDeviceToHost), "download");
}
template <class T> void copy(T* dst, const T* src, size_t n) {
  if (n) check(cudaMemcpyAsync(dst, src, n * sizeof(T), cudaMemcpyDeviceToDevice), "copy");
}
template <class T> void zero(T* d, size_t n) {
  if (n) check(cudaMemsetAsync(d, 0, n * sizeof(T)), "memset");
}
inline void fence() { check(cudaDeviceSynchronize(), "kernel execution"); }

template <class F> __global__ void kFor(long n, F f) {
  for (long i = long(blockIdx.x) * blockDim.x + threadIdx.x; i < n; i += long(blockDim.x) * gridDim.x) f(i);
}
template <class F> __global__ void kForDyn(const int* a, const int* b, F f) {
  const long n = long(*a) + (b ? long(*b) : 0L);
  for (long i = long(blockIdx.x) * blockDim.x + threadIdx.x; i < n; i += long(blockDim.x) * gridDim.x) f(i);
}
template <class F> __global__ void kTeams(F f) {
  extern __shared__ double smem[];
  Team t{int(blockIdx.x), int(threadIdx.x), int(blockDim.x), smem};
  f(t);
}
inline int blocksFor(long n) {
  long b = (n + BLOCK - 1) / BLOCK;
  return int(b < 1 ? 1 : (b > 8192 ? 8192 : b));
}
template <class F> void forN(long n, F f) {
  if (n <= 0) return;
  kFor<<<blocksFor(n), BLOCK>>>(n, f);
  check(cudaGetLastError(), "launch");
}
template <class F> void forDyn(const int* a, const int* b, long cap, F f, bool /*ordered*/ = false) {
  if (cap <= 0) return;
  kForDyn<<<blocksFor(cap), BLOCK>>>(a, b, f);
  check(cudaGetLastError(), "launch");
}
template <class F> void teams(int nt, F f) {
  kTeams<<<nt, TEAM, TEAM * sizeof(double)>>>(f);
  check(cudaGetLastError(), "launch");
}
// check driver and device before the first allocation, with an explanation
inline void initDevice() {
  int drv = 0, rt = 0;
  cudaDriverGetVersion(&drv);
  cudaRuntimeGetVersion(&rt);
  char buf[512];
  if (drv == 0)
    throw std::runtime_error("no NVIDIA driver found: backend = cuda needs a machine with an NVIDIA GPU and its "
                             "driver (check with nvidia-smi); on this machine use backend = serial");
  if (drv < rt) {
    std::snprintf(buf, sizeof buf,
                  "the NVIDIA driver supports CUDA up to %d.%d, but this binary was built with CUDA %d.%d. Build with "
                  "a CUDA toolkit <= %d.%d (see docs/gpu.md) or update the driver",
                  drv / 1000, (drv % 1000) / 10, rt / 1000, (rt % 1000) / 10, drv / 1000, (drv % 1000) / 10);
    throw std::runtime_error(buf);
  }
  int n = 0;
  const cudaError_t e = cudaGetDeviceCount(&n);
  if (e != cudaSuccess || n == 0)
    throw std::runtime_error(std::string("no usable CUDA device: ") +
                             (e != cudaSuccess ? cudaGetErrorString(e) : "the driver reports 0 GPUs"));
  // kernels for this GPU in the binary? (a missing architecture shows up at the first launch)
}

inline std::string describe() {
  int dev = 0, n = 0;
  check(cudaGetDeviceCount(&n), "cudaGetDeviceCount");
  if (n == 0) throw std::runtime_error("no CUDA device found");
  check(cudaGetDevice(&dev), "cudaGetDevice");
  cudaDeviceProp p;
  check(cudaGetDeviceProperties(&p, dev), "cudaGetDeviceProperties");
  char buf[512];
  std::snprintf(buf, sizeof buf, "CUDA device %d: %s, %d SMs, %.1f GB, compute capability %d.%d", dev, p.name,
                p.multiProcessorCount, p.totalGlobalMem / 1073741824.0, p.major, p.minor);
  return buf;
}

#else
// ============================================================================
//  CPU backend (OpenMP threads if enabled)
// ============================================================================
template <class T> T* alloc(size_t n) { return n ? static_cast<T*>(::operator new(n * sizeof(T))) : nullptr; }
template <class T> void release(T* p) { ::operator delete(static_cast<void*>(p)); }
template <class T> void upload(T* d, const T* h, size_t n) { if (n) std::memcpy(d, h, n * sizeof(T)); }
template <class T> void download(T* h, const T* d, size_t n) { if (n) std::memcpy(h, d, n * sizeof(T)); }
template <class T> void copy(T* dst, const T* src, size_t n) { if (n) std::memmove(dst, src, n * sizeof(T)); }
template <class T> void zero(T* d, size_t n) { if (n) std::memset(static_cast<void*>(d), 0, n * sizeof(T)); }
inline void fence() {}

template <class F> void forN(long n, F f) {
#if defined(_OPENMP)
#pragma omp parallel for schedule(static) if (n > 512)
#endif
  for (long i = 0; i < n; ++i) f(i);
}
template <class F> void forDyn(const int* a, const int* b, long /*cap*/, F f, bool ordered = false) {
  const long n = long(*a) + (b ? long(*b) : 0L);
#if defined(GEM_EMULATE_GPU)
  ordered = false;  // as on the GPU
#endif
  if (ordered) {
    for (long i = 0; i < n; ++i) f(i);
    return;
  }
#if defined(_OPENMP)
#pragma omp parallel for schedule(dynamic, 64) if (n > 256)
#endif
  for (long i = 0; i < n; ++i) f(i);
}
#if defined(GEM_EMULATE_GPU) && defined(_OPENMP)
#ifndef GEM_EMULATE_TEAM
#define GEM_EMULATE_TEAM 8
#endif
template <class F> void teams(int nt, F f) {
  for (int t = 0; t < nt; ++t) {
    double scratch[GEM_EMULATE_TEAM];
#pragma omp parallel num_threads(GEM_EMULATE_TEAM)
    {
      if (omp_get_num_threads() != GEM_EMULATE_TEAM) throw std::runtime_error("team emulation: thread count");
      f(Team{t, omp_get_thread_num(), GEM_EMULATE_TEAM, scratch});
    }
  }
}
#else
template <class F> void teams(int nt, F f) {
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1) if (nt > 1)
#endif
  for (int t = 0; t < nt; ++t) {
    double s = 0.0;
    f(Team{t, 0, 1, &s});
  }
}
#endif
inline void initDevice() {}

inline std::string describe() {
#if defined(GEM_EMULATE_GPU)
  return "CPU emulation of the GPU code paths (teams of " + std::to_string(GEM_EMULATE_TEAM) + " threads)";
#elif defined(_OPENMP)
  return "CPU, OpenMP with " + std::to_string(omp_get_max_threads()) + " threads";
#else
  return "CPU, 1 thread (built without OpenMP)";
#endif
}
#endif

// ---------------------------------------------------------------- device array
template <class T>
class Array {
 public:
  Array() = default;
  explicit Array(size_t n) { reserve(n); }
  Array(const Array&) = delete;
  Array& operator=(const Array&) = delete;
  ~Array() { release(p_); }
  // grow to at least n elements; keep = preserve the old contents
  void reserve(size_t n, bool keep = false) {
    if (n <= cap_) return;
    const size_t nc = n > cap_ + cap_ / 2 ? n : cap_ + cap_ / 2;
    T* q = alloc<T>(nc);
    if (keep && p_) copy(q, p_, cap_);
    release(p_);
    p_ = q;
    cap_ = nc;
  }
  void assign(const T* host, size_t n) { reserve(n); upload(p_, host, n); }
  void clear(size_t n) { reserve(n); zero(p_, n); }
  T* data() { return p_; }
  const T* data() const { return p_; }
  size_t capacity() const { return cap_; }
  void swap(Array& o) { std::swap(p_, o.p_); std::swap(cap_, o.cap_); }

 private:
  T* p_ = nullptr;
  size_t cap_ = 0;
};

}  // namespace par
}  // namespace gem
