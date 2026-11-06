#include <doctest/doctest.h>

#include <cstdint>
#include <limits>
#include <stdexcept>

#include <geom/core/arena.hpp>
#include <geom/core/handle.hpp>

using geom::v1::Handle;
using geom::v1::HandleAllocator;
using geom::v1::ScratchArena;

TEST_CASE("Handle packs index and generation") {
  const Handle h = Handle::make(1234u, 7u);
  CHECK(h.index() == 1234u);
  CHECK(h.generation() == 7u);
  CHECK_FALSE(h.is_null());
  CHECK(Handle::null().is_null());
  CHECK(h == Handle::make(1234u, 7u));
  CHECK_FALSE(h == Handle::make(1235u, 7u));
}

TEST_CASE("HandleAllocator invalidates stale handles on reuse") {
  HandleAllocator alloc;
  const Handle a = alloc.alloc();
  const Handle b = alloc.alloc();
  REQUIRE(alloc.alive(a));
  REQUIRE(alloc.alive(b));

  alloc.release(a);
  CHECK_FALSE(alloc.alive(a));

  const Handle c = alloc.alloc();
  CHECK(c.index() == a.index());
  CHECK(c.generation() != a.generation());
  CHECK(alloc.alive(c));
  CHECK_FALSE(alloc.alive(a));
}

TEST_CASE("ScratchArena allocates aligned arrays and resets") {
  ScratchArena arena;
  auto a = arena.alloc_array<double>(100);
  REQUIRE(a.size() == 100);
  a[0] = 1.5;
  a[99] = 2.5;
  CHECK(a[0] == 1.5);

  arena.reset();
  auto b = arena.alloc_array<double>(10);
  REQUIRE(b.size() == 10);
  b[0] = 7.0;
  CHECK(b[0] == 7.0);
}

TEST_CASE("HandleAllocator ignores null and stale releases") {
  HandleAllocator alloc;
  CHECK_FALSE(alloc.alive(Handle::null()));
  alloc.release(Handle::null());  // must be a no-op, not a free-list entry

  const Handle a = alloc.alloc();
  alloc.release(a);
  alloc.release(a);  // double release must not enqueue the slot twice
  const Handle b = alloc.alloc();
  CHECK(b.index() == a.index());
  CHECK(alloc.alive(b));

  alloc.release(a);  // stale generation must not free the live slot
  CHECK(alloc.alive(b));
}

TEST_CASE("ScratchArena returns distinct, aligned storage") {
  ScratchArena arena;
  auto x = arena.alloc_array<double>(8);
  auto y = arena.alloc_array<double>(8);
  CHECK(x.data() != y.data());
  CHECK(reinterpret_cast<std::uintptr_t>(x.data()) % alignof(double) == 0);
  CHECK(reinterpret_cast<std::uintptr_t>(y.data()) % alignof(double) == 0);
  x[0] = 1.0;
  y[0] = 2.0;
  CHECK(x[0] == 1.0);
  CHECK(y[0] == 2.0);
}

TEST_CASE("Handle::make rejects out-of-range indices") {
  // 0xFFFFFF is the reserved null index: make must canonicalize, not hand back
  // a handle that merely shares the index bits.
  CHECK(Handle::make(0xFFFFFFu, 1u) == Handle::null());
  // 0x1000000 would shift out of the index field and truncate to index 0.
  CHECK(Handle::make(0x1000000u, 3u).is_null());
  CHECK(Handle::make(0xFFFFFFFFu, 0u).is_null());
  // The largest usable index still packs both fields.
  const Handle edge = Handle::make(0xFFFFFEu, 7u);
  CHECK(edge.index() == 0xFFFFFEu);
  CHECK(edge.generation() == 7u);
  CHECK_FALSE(edge.is_null());
}

TEST_CASE("ScratchArena rejects overflowing array requests") {
  ScratchArena arena;
  const std::size_t count = (std::numeric_limits<std::size_t>::max)() / sizeof(double) + 2;
  CHECK_THROWS_AS(arena.alloc_array<double>(count), std::length_error);

  auto empty = arena.alloc_array<double>(0);
  CHECK(empty.empty());
  CHECK(empty.size() == 0);

  // A rejected request must not poison the arena.
  auto ok = arena.alloc_array<double>(4);
  CHECK(ok.size() == 4);
  ok[3] = 3.0;
  CHECK(ok[3] == 3.0);
}

TEST_CASE("HandleAllocator generation lifecycle is unchanged") {
  HandleAllocator alloc;
  const Handle first = alloc.alloc();
  std::uint8_t previous = first.generation();
  Handle h = first;

  for (int cycle = 0; cycle < 8; ++cycle) {
    alloc.release(h);
    CHECK_FALSE(alloc.alive(h));
    h = alloc.alloc();
    REQUIRE(h.index() == first.index());
    CHECK(h.generation() != previous);
    previous = h.generation();
  }

  // The 8-bit generation wraps after 255 releases; the slot must stay reusable.
  for (int cycle = 0; cycle < 300; ++cycle) {
    alloc.release(h);
    h = alloc.alloc();
    REQUIRE(h.index() == first.index());
    REQUIRE(alloc.alive(h));
  }
}

