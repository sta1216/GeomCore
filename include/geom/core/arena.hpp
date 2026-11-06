#pragma once

#include <cstddef>
#include <limits>
#include <memory_resource>
#include <span>
#include <stdexcept>

namespace geom {
inline namespace v1 {

// Monotonic scratch allocator for one algorithm call. Memory is released as a
// whole by reset(); individual deallocations are no-ops by design.
//
// alloc_array<T> returns *uninitialized* storage: T objects are neither
// constructed nor destroyed, so T must be an implicit-lifetime type with a
// trivial destructor (doubles, indices, plain structs). reset() invalidates
// every span previously handed out by this arena; using one afterwards is UB.
class ScratchArena {
 public:
  explicit ScratchArena(std::size_t initial_bytes = 64 * 1024)
      : resource_(initial_bytes, std::pmr::new_delete_resource()) {}

  std::pmr::memory_resource* resource() noexcept { return &resource_; }

  void reset() noexcept { resource_.release(); }

  // count == 0 returns an empty span without touching the arena. A request
  // whose byte size does not fit in std::size_t throws std::length_error
  // before any allocation: an overflowing count * sizeof(T) must never reach
  // the resource, which would hand back a span larger than the real storage.
  template <class T>
  std::span<T> alloc_array(std::size_t count) {
    if (count == 0) return {};
    if (count > (std::numeric_limits<std::size_t>::max)() / sizeof(T)) {
      throw std::length_error("ScratchArena::alloc_array: count * sizeof(T) overflows");
    }
    T* data = static_cast<T*>(resource_.allocate(count * sizeof(T), alignof(T)));
    return {data, count};
  }

 private:
  std::pmr::monotonic_buffer_resource resource_;
};

}  // namespace v1
}  // namespace geom
