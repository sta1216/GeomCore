#include <doctest/doctest.h>

#include <string_view>
#include <type_traits>
#include <variant>

#include <geom/brep/topology.hpp>
#include <geom/core/error.hpp>

TEST_CASE("empty brep has no live entities") {
  geom::Brep b;
  CHECK(b.vertex_count() == 0);
  CHECK(b.edge_count() == 0);
  CHECK(b.face_count() == 0);
  CHECK(b.solid_count() == 0);
  CHECK_FALSE(b.vertex_alive(geom::VertexId::null()));
  CHECK(b.geometry().curves.empty());
}

TEST_CASE("brep ids are generation-checked handles") {
  static_assert(std::is_same_v<geom::VertexId, geom::Handle>);
  const geom::VertexId v = geom::VertexId::make(3, 2);
  CHECK(v.index() == 3);
  CHECK(v.generation() == 2);
  CHECK_FALSE(v.is_null());
  CHECK(geom::VertexId::null().is_null());
}

TEST_CASE("geometry store preserves analytic type identity") {
  geom::GeometryStore g;
  g.surfaces.emplace_back(geom::Plane3{{0, 0, 1}, 0.0});
  CHECK(std::holds_alternative<geom::Plane3>(g.surfaces.back()));
  g.surfaces.emplace_back(geom::Cylinder3{});
  CHECK(std::holds_alternative<geom::Cylinder3>(g.surfaces.back()));
  g.curves.emplace_back(geom::CircleCurve3{});
  CHECK(std::holds_alternative<geom::CircleCurve3>(g.curves.back()));
}

TEST_CASE("new error codes stringify") {
  CHECK(std::string_view(geom::to_string(geom::Errc::invalid_topology)) == "invalid_topology");
  CHECK(std::string_view(geom::to_string(geom::Errc::tolerance_violation)) ==
        "tolerance_violation");
  CHECK(std::string_view(geom::to_string(geom::Errc::not_closed)) == "not_closed");
}

TEST_CASE("null and out-of-range handles are never alive") {
  geom::Brep b;
  CHECK_FALSE(b.edge_alive(geom::EdgeId::null()));
  CHECK_FALSE(b.coedge_alive(geom::CoedgeId::null()));
  CHECK_FALSE(b.loop_alive(geom::LoopId::null()));
  CHECK_FALSE(b.face_alive(geom::FaceId::null()));
  CHECK_FALSE(b.shell_alive(geom::ShellId::null()));
  CHECK_FALSE(b.solid_alive(geom::SolidId::null()));
  CHECK_FALSE(b.vertex_alive(geom::VertexId::make(7, 1)));
  CHECK_FALSE(b.edge_alive(geom::EdgeId::make(7, 1)));
  CHECK_FALSE(b.solid_alive(geom::SolidId::make(7, 1)));

  CHECK(b.vertices().empty());
  CHECK(b.edges().empty());
  CHECK(b.coedges().empty());
  CHECK(b.loops().empty());
  CHECK(b.faces().empty());
  CHECK(b.shells().empty());
  CHECK(b.solids().empty());
  CHECK(b.solid_ids().empty());
  CHECK(b.coedge_count() == 0);
  CHECK(b.loop_count() == 0);
  CHECK(b.shell_count() == 0);
  CHECK(b.geometry().surfaces.empty());
  CHECK(b.geometry().pcurves.empty());
}

TEST_CASE("brep is a copyable value type") {
  static_assert(std::is_copy_constructible_v<geom::Brep>);
  static_assert(std::is_copy_assignable_v<geom::Brep>);
  static_assert(std::is_move_constructible_v<geom::Brep>);
  geom::Brep a;
  geom::Brep b = a;
  CHECK(b.vertex_count() == 0);
  CHECK(b.geometry().curves.empty());
  a = b;
  CHECK(a.solid_count() == 0);
}

TEST_CASE("geometry slots default to kNoGeometry") {
  CHECK(geom::kNoGeometry == 0xFFFFFFFFu);
  CHECK(geom::Edge{}.curve == geom::kNoGeometry);
  CHECK(geom::Coedge{}.pcurve == geom::kNoGeometry);
  CHECK(geom::Face{}.surface == geom::kNoGeometry);
}
