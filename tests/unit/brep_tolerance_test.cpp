#include <doctest/doctest.h>

#include <array>
#include <cstdint>
#include <string_view>

#include <geom/brep/solid_editor.hpp>
#include <geom/brep/validate.hpp>

namespace {

geom::SolidEditor make_square_editor(double vertex_tol, double edge_tol,
                                     double face_tol, double solid_tol) {
  geom::SolidEditor ed;
  const geom::Point3d pts[4] = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}};
  std::array<geom::VertexId, 4> v{};
  for (int i = 0; i < 4; ++i) v[i] = ed.add_vertex(pts[i], vertex_tol).value();
  const std::uint32_t plane = ed.add_surface(geom::Plane3{{0, 0, 1}, 0.0});
  const std::uint32_t line = ed.add_curve(geom::LineCurve3{{0, 0, 0}, {1, 0, 0}});
  std::array<geom::EdgeId, 4> e{};
  for (int i = 0; i < 4; ++i) {
    e[i] = ed.add_edge(v[i], v[(i + 1) % 4], line, 0.0, 1.0, edge_tol).value();
  }
  const geom::FaceId f = ed.add_face(plane, false, face_tol).value();
  std::array<geom::CoedgeSpec, 4> specs{};
  for (int i = 0; i < 4; ++i) specs[i] = geom::CoedgeSpec{e[i], false};
  ed.add_loop(f, specs, true).value();
  const auto sh = ed.add_shell(std::array<geom::FaceId, 1>{f}, false).value();
  ed.add_solid(std::array<geom::ShellId, 1>{sh}, solid_tol).value();
  return ed;
}

}  // namespace

TEST_CASE("edge tolerance must dominate its vertices") {
  auto ed = make_square_editor(/*vertex=*/0.1, /*edge=*/0.05, /*face=*/0.2,
                               /*solid=*/0.2);
  const auto r = ed.commit();
  REQUIRE_FALSE(r.has_value());
  CHECK(r.error().code == geom::Errc::tolerance_violation);
}

TEST_CASE("solid tolerance must dominate its faces") {
  auto ed = make_square_editor(0.0, 0.0, /*face=*/1.0, /*solid=*/0.1);
  const auto r = ed.commit();
  REQUIRE_FALSE(r.has_value());
  CHECK(r.error().code == geom::Errc::tolerance_violation);
}

TEST_CASE("tolerance growth beyond ten times is a violation") {
  auto ed = make_square_editor(0.0, /*edge=*/0.1, /*face=*/2.0, /*solid=*/2.0);
  const auto r = ed.commit();
  REQUIRE_FALSE(r.has_value());
  CHECK(r.error().code == geom::Errc::tolerance_violation);
  CHECK(std::string_view(r.error().message).find("growth") != std::string_view::npos);
}

TEST_CASE("valid tolerance chain is accepted") {
  auto ed = make_square_editor(0.01, 0.02, 0.03, 0.03);
  const auto r = ed.commit();
  REQUIRE(r.has_value());
  CHECK(geom::validate_tolerances(r.value()).ok());
}

TEST_CASE("negative tolerance is rejected at construction") {
  auto ed = make_square_editor(0.0, 0.0, 0.0, 0.0);
  CHECK_FALSE(ed.add_vertex({0, 0, 0}, -1.0).has_value());
  CHECK_FALSE(ed.add_face(0, false, -0.5).has_value());
}
