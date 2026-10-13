// ============================================================================
//  eigen.hpp -- eigenvalues and eigenvectors of a real symmetric tridiagonal
//  matrix (implicit QL with Wilkinson shifts, after EISPACK tql2)
// ============================================================================
#pragma once
#include <vector>

namespace gem {

// d: diagonal (n), e: sub-diagonal (e[i] couples i and i+1, n-1 used).
// On return d holds the eigenvalues and Z (n*n, row-major, Z[i*n+k]) the
// orthonormal eigenvectors: column k belongs to d[k]. Returns false if the
// iteration did not converge.
bool symTridiagEigen(std::vector<double>& d, std::vector<double> e, std::vector<double>& Z);

}  // namespace gem
