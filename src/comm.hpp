// ============================================================================
//  comm.hpp -- thin MPI layer (point-to-point messages of packed byte buffers)
//
//  Without MPI (GEMPIC_HAVE_MPI undefined) there is one rank and no messages.
// ============================================================================
#pragma once
#include <cstddef>
#include <cstring>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

namespace gem {

// packing of trivially copyable data
class Buffer {
 public:
  template <class T> void put(const T& v) {
    static_assert(std::is_trivially_copyable<T>::value, "put: trivially copyable types only");
    putRaw(&v, sizeof(T));
  }
  template <class T> void putVec(const std::vector<T>& v) {
    put<unsigned long long>(v.size());
    if (!v.empty()) putRaw(v.data(), v.size() * sizeof(T));
  }
  void putRaw(const void* p, size_t n) {
    const size_t o = data.size();
    data.resize(o + n);
    if (n) std::memcpy(&data[o], p, n);
  }
  template <class T> T get() {
    T v;
    getRaw(&v, sizeof(T));
    return v;
  }
  template <class T> std::vector<T> getVec() {
    std::vector<T> v(size_t(get<unsigned long long>()));
    if (!v.empty()) getRaw(v.data(), v.size() * sizeof(T));
    return v;
  }
  void getRaw(void* p, size_t n) {
    if (pos + n > data.size()) throw std::runtime_error("message buffer: read beyond the end");
    if (n) std::memcpy(p, &data[pos], n);
    pos += n;
  }
  bool done() const { return pos == data.size(); }
  std::vector<char> data;
  size_t pos = 0;
};

namespace comm {
void init(int* argc, char*** argv);
void finalize();
int rank();
int size();
bool available();                         // built with MPI
void send(int dest, int tag, Buffer&& b);  // non-blocking; the buffer is kept until the send completes
Buffer recv(int src, int tag);            // blocking
void waitAll();                           // complete all pending sends
[[noreturn]] void abort(int code);
}  // namespace comm

}  // namespace gem
