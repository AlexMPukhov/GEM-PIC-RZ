// ============================================================================
//  compare_columns.cpp -- compare two numeric ASCII tables column by column
//
//  compare_columns <reference> <result> [rtol] [atol] [ncols]
//    lines starting with '#' and blank lines are skipped;
//    |a - b| <= atol + rtol * max(|a|, |b|) must hold for every value
//    (only the first ncols columns are compared if ncols > 0).
//  Exit code 0: files agree, 1: difference, 2: usage / file error.
// ============================================================================
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

static bool readTable(const char* path, std::vector<std::vector<double>>& rows) {
  std::ifstream in(path);
  if (!in) return false;
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#') continue;
    std::istringstream ss(line);
    std::vector<double> row;
    std::string tok;
    while (ss >> tok) {  // non-numeric tokens (e.g. bunch names) are skipped
      char* end = nullptr;
      const double x = std::strtod(tok.c_str(), &end);
      if (end != tok.c_str() && *end == '\0') row.push_back(x);
    }
    if (!row.empty()) rows.push_back(row);
  }
  return true;
}

int main(int argc, char** argv) {
  if (argc < 3) {
    std::fprintf(stderr, "usage: compare_columns <reference> <result> [rtol=1e-9] [atol=1e-12] [ncols=0]\n");
    return 2;
  }
  const double rtol = argc > 3 ? std::atof(argv[3]) : 1e-9;
  const double atol = argc > 4 ? std::atof(argv[4]) : 1e-12;
  const size_t ncols = argc > 5 ? size_t(std::atoi(argv[5])) : 0;
  std::vector<std::vector<double>> A, B;
  if (!readTable(argv[1], A)) { std::fprintf(stderr, "cannot read %s\n", argv[1]); return 2; }
  if (!readTable(argv[2], B)) { std::fprintf(stderr, "cannot read %s\n", argv[2]); return 2; }
  if (A.size() != B.size()) {
    std::printf("FAIL %s: %zu rows vs %zu rows\n", argv[2], A.size(), B.size());
    return 1;
  }
  double worst = 0.0;
  size_t wr = 0, wc = 0, bad = 0;
  for (size_t r = 0; r < A.size(); ++r) {
    size_t nc = std::min(A[r].size(), B[r].size());
    if (ncols > 0) nc = std::min(nc, ncols);
    for (size_t c = 0; c < nc; ++c) {
      const double a = A[r][c], b = B[r][c];
      const double d = std::fabs(a - b), lim = atol + rtol * std::max(std::fabs(a), std::fabs(b));
      if (d > lim) ++bad;
      const double rel = d / (lim > 0 ? lim : 1.0);
      if (rel > worst) { worst = rel; wr = r; wc = c; }
    }
  }
  if (bad) {
    std::printf("FAIL %s: %zu values differ; worst row %zu col %zu: %.9e vs %.9e\n", argv[2], bad, wr + 1, wc + 1,
                A[wr][wc], B[wr][wc]);
    return 1;
  }
  std::printf("OK   %s: %zu rows agree (worst |diff| / tolerance = %.3g)\n", argv[2], A.size(), worst);
  return 0;
}
