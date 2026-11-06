#pragma once

#include <cstdint>
#include <vector>

namespace geom {
inline namespace v1 {

// 32-bit packed handle: 24-bit index + 8-bit generation. Index 0xFFFFFF is
// reserved as the null handle.
class Handle {
 public:
  static constexpr std::uint32_t kIndexBits = 24;
  static constexpr std::uint32_t kIndexMask = (1u << kIndexBits) - 1u;
  static constexpr std::uint8_t kGenerationMask = 0xFFu;

  constexpr Handle() noexcept = default;

  static constexpr Handle make(std::uint32_t index, std::uint8_t generation) noexcept {
    return Handle((index << 8) | generation);
  }
  static constexpr Handle null() noexcept { return Handle(kIndexMask << 8); }

  constexpr std::uint32_t index() const noexcept { return packed_ >> 8; }
  constexpr std::uint8_t generation() const noexcept {
    return static_cast<std::uint8_t>(packed_ & kGenerationMask);
  }
  constexpr bool is_null() const noexcept { return index() == kIndexMask; }

  friend constexpr bool operator==(Handle a, Handle b) noexcept { return a.packed_ == b.packed_; }
  friend constexpr bool operator!=(Handle a, Handle b) noexcept { return !(a == b); }

 private:
  explicit constexpr Handle(std::uint32_t packed) noexcept : packed_(packed) {}
  std::uint32_t packed_ = kIndexMask << 8;
};

class HandleAllocator {
 public:
  Handle alloc() {
    if (!free_.empty()) {
      const std::uint32_t index = free_.back();
      free_.pop_back();
      return Handle::make(index, generations_[index]);
    }
    if (generations_.size() >= Handle::kIndexMask) return Handle::null();
    const std::uint32_t index = static_cast<std::uint32_t>(generations_.size());
    generations_.push_back(1u);
    return Handle::make(index, 1u);
  }

  void release(Handle h) {
    if (h.is_null() || !alive(h)) return;
    generations_[h.index()] = next_generation(h.index());
    free_.push_back(h.index());
  }

  bool alive(Handle h) const noexcept {
    if (h.is_null() || h.index() >= generations_.size()) return false;
    return generations_[h.index()] == h.generation();
  }

 private:
  std::uint8_t next_generation(std::uint32_t index) const noexcept {
    const std::uint8_t current = generations_[index];
    return current == 0xFFu ? 1u : static_cast<std::uint8_t>(current + 1u);
  }

  std::vector<std::uint8_t> generations_;
  std::vector<std::uint32_t> free_;
};

}  // namespace v1
}  // namespace geom
