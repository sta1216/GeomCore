#pragma once

#include <cstdint>

namespace geom {
inline namespace v1 {

enum class Errc : std::uint8_t {
  none = 0,
  degenerate,
  non_finite,
  index_out_of_range,
  numeric_failure,
  io_error,
  unsupported,
  cancelled,
};

struct Error {
  Errc code = Errc::none;
  const char* message = "";
  std::uint32_t context_a = 0;
  std::uint32_t context_b = 0;
};

constexpr const char* to_string(Errc code) noexcept {
  switch (code) {
    case Errc::none: return "none";
    case Errc::degenerate: return "degenerate";
    case Errc::non_finite: return "non_finite";
    case Errc::index_out_of_range: return "index_out_of_range";
    case Errc::numeric_failure: return "numeric_failure";
    case Errc::io_error: return "io_error";
    case Errc::unsupported: return "unsupported";
    case Errc::cancelled: return "cancelled";
  }
  return "unknown";
}

}  // namespace v1
}  // namespace geom
