// ============================================================================
//  hd.hpp -- GEM_HD marks functions that run on the host and on a CUDA device
//  (empty for ordinary C++ compilers)
// ============================================================================
#pragma once
#if defined(__CUDACC__)
#define GEM_HD __host__ __device__
#else
#define GEM_HD
#endif
