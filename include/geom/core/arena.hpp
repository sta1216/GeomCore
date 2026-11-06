#pragma once

#include <cstddef>
#include <memory_resource>
#include <span>

namespace geom {
inline namespace v1 {

// Monotonic scratch allocator for one algorithm call. Memory is released as a
// whole by reset(); individual deallocations are no-ops by design.
class ScratchArena {
 public:
  explicit ScratchArena(std::size_t initial_bytes = 64 * 1024)
      : resource_(initial_bytes, std::pmr::new_delete_resource()) {}

  std::pmr::memory_resource* resource() noexcept { return &resource_; }

  void reset() noexcept { resource_.release(); }

  template <class T>
  std::span<T> alloc_array(std::size_t count) {
    T* data = static_cast<T*>(resource_.allocate(count * sizeof(T), alignof(T)));
    return {data, count};
  }

 private:
  std::pmr::monotonic_buffer_resource resource_;
};

}  // namespace v1
}  // namespace geom
