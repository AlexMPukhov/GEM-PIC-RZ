// ============================================================================
//  comm.cpp -- MPI layer (see comm.hpp)
// ============================================================================
#include "comm.hpp"

#include <cstdlib>
#include <list>

#if defined(GEMPIC_HAVE_MPI)
#define OMPI_SKIP_MPICXX 1  // only the C API is used: skip the deprecated C++ bindings
#define MPICH_SKIP_MPICXX 1
#include <mpi.h>
#endif

namespace gem {
namespace comm {

#if defined(GEMPIC_HAVE_MPI)

namespace {
struct Pending {
  Buffer buf;
  MPI_Request req;
};
std::list<Pending> pendingSends;
bool initialized = false;

void testPending() {
  for (auto it = pendingSends.begin(); it != pendingSends.end();) {
    int flag = 0;
    MPI_Test(&it->req, &flag, MPI_STATUS_IGNORE);
    if (flag) it = pendingSends.erase(it);
    else ++it;
  }
}
}  // namespace

void init(int* argc, char*** argv) {
  int done = 0;
  MPI_Initialized(&done);
  if (!done) MPI_Init(argc, argv);
  initialized = true;
}
void finalize() {
  if (!initialized) return;
  waitAll();
  MPI_Finalize();
  initialized = false;
}
int rank() {
  int r = 0;
  if (initialized) MPI_Comm_rank(MPI_COMM_WORLD, &r);
  return r;
}
int size() {
  int s = 1;
  if (initialized) MPI_Comm_size(MPI_COMM_WORLD, &s);
  return s;
}
bool available() { return true; }

void send(int dest, int tag, Buffer&& b) {
  if (b.data.size() > size_t(2147483647)) throw std::runtime_error("MPI message larger than 2 GB");
  pendingSends.push_back(Pending{std::move(b), MPI_REQUEST_NULL});
  Pending& p = pendingSends.back();
  MPI_Isend(p.buf.data.data(), int(p.buf.data.size()), MPI_BYTE, dest, tag, MPI_COMM_WORLD, &p.req);
  testPending();
}

Buffer recv(int src, int tag) {
  MPI_Status st;
  MPI_Probe(src, tag, MPI_COMM_WORLD, &st);
  int n = 0;
  MPI_Get_count(&st, MPI_BYTE, &n);
  Buffer b;
  b.data.resize(size_t(n));
  MPI_Recv(b.data.data(), n, MPI_BYTE, src, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  testPending();
  return b;
}

void waitAll() {
  for (auto& p : pendingSends) MPI_Wait(&p.req, MPI_STATUS_IGNORE);
  pendingSends.clear();
}

void abort(int code) {
  if (initialized) MPI_Abort(MPI_COMM_WORLD, code);
  std::exit(code);
}

#else  // no MPI: one rank

void init(int*, char***) {}
void finalize() {}
int rank() { return 0; }
int size() { return 1; }
bool available() { return false; }
void send(int, int, Buffer&&) { throw std::runtime_error("no MPI in this build"); }
Buffer recv(int, int) { throw std::runtime_error("no MPI in this build"); }
void waitAll() {}
void abort(int code) { std::exit(code); }

#endif

}  // namespace comm
}  // namespace gem
