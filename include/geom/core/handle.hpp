#pragma once

#include <cstdint>
#include <vector>

namespace geom {
inline namespace v1 {

// 32-bit packed handle: 24-bit index + 8-bit generation. Index 0xFFFFFF is the
// reserved null index and is never a valid live index.
//
// The generation counter is 8-bit and is bumped on every release, wrapping
// 0xFF -> 1. A handle that is 255 (or more) releases old can therefore compare
// alive again (ABA). Accepted for M0; widening the generation field is deferred.
class Handle {
 public:
  static constexpr std::uint32_t kIndexBits = 24;
  static constexpr std::uint32_t kNullIndex = 0xFFFFFFu;
  static constexpr std::uint32_t kIndexMask = (1u << kIndexBits) - 1u;
  static constexpr std::uint8_t kGenerationMask = 0xFFu;

  constexpr Handle() noexcept = default;

  // Precondition: index < kNullIndex. Out-of-range requests (including the
  // reserved null index) return null(); the index is never truncated.
  static constexpr Handle make(std::uint32_t index, std::uint8_t generation) noexcept {
    if (index >= kNullIndex) return null();
    return Handle((index << 8) | generation);
  }
  static constexpr Handle null() noexcept { return Handle(kIndexMask << 8); }

  constexpr std::uint32_t index() const noexcept { return packed_ >> 8; }
  constexpr std::uint8_t generation() const noexcept {
    return static_cast<std::uint8_t>(packed_ & kGenerationMask);
  }
  constexpr bool is_null() const noexcept { return index() == kNullIndex; }

  friend constexpr bool operator==(Handle a, Handle b) noexcept { return a.packed_ == b.packed_; }
  friend constexpr bool operator!=(Handle a, Handle b) noexcept { return !(a == b); }

 private:
  explicit constexpr Handle(std::uint32_t packed) noexcept : packed_(packed) {}
  std::uint32_t packed_ = kIndexMask << 8;
};

class HandleAllocator {
 public:
  // Returns Handle::null() when all kNullIndex usable slots are in use.
  Handle alloc() {
    if (!free_.empty()) {
      const std::uint32_t index = free_.back();
      free_.pop_back();
      return Handle::make(index, generations_[index]);
    }
    if (generations_.size() >= Handle::kNullIndex) return Handle::null();
    const std::uint32_t index = static_cast<std::uint32_t>(generations_.size());
    generations_.push_back(1u);
    return Handle::make(index, 1u);
  }

  void release(Handle h) {
    if (h.is_null() || !alive(h)) return;
    const std::uint32_t index = h.index();
    const std::uint8_t next = next_generation(index);
    // Enqueue before bumping: push_back can throw bad_alloc, and a slot that
    // was bumped but never enqueued would be leaked as permanently alive.
    free_.push_back(index);
    generations_[index] = next;
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
