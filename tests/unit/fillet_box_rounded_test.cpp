#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <span>
#include <string_view>
#include <utility>
#include <variant>
#include <vector>

#include <doctest/doctest.h>

#include <geom/brep/detail/curved_volume.hpp>
#include <geom/brep/fillet.hpp>
#include <geom/brep/mass_properties.hpp>
#include <geom/brep/primitives.hpp>
#include <geom/brep/solid_editor.hpp>
#include <geom/brep/tessellate.hpp>
#include <geom/brep/topology.hpp>
#include <geom/brep/validate.hpp>
#include <geom/core/error.hpp>
#include <geom/curves/analytic_curve.hpp>
#include <geom/surfaces/analytic.hpp>

namespace {
using geom::Brep;
using geom::EdgeId;
using geom::Point3d;

Brep box_2x3x4() {
  return geom::make_box(Point3d{0.0, 0.0, 0.0}, Point3d{2.0, 3.0, 4.0}, 1e-9).value();
}

std::vector<EdgeId> all_edges(const Brep& brep) {
  const auto ids = brep.edge_ids();
  return std::vector<EdgeId>(ids.begin(), ids.end());
}

double volume_of(const Brep& brep) {
  return geom::detail::curved_solid_volume(brep, brep.solid_ids().front(), 1e-9).value();
}

geom::Vec3d unit_vector(const geom::Vec3d& v) {
  return v * (1.0 / std::sqrt(geom::dot(v, v)));
}

// One spherical face bounded by three great-circle arcs a -> b -> c -> a. The
// `reversed` flag lets a test ask for the same trim with the face normal flipped.
Brep sphere_face_from_corners(const geom::Sphere3& sphere, const Point3d& a, const Point3d& b,
                              const Point3d& c, bool reversed) {
  geom::SolidEditor editor;
  const std::uint32_t surface = editor.add_surface(sphere);
  const geom::VertexId va = editor.add_vertex(a, 0.0).value();
  const geom::VertexId vb = editor.add_vertex(b, 0.0).value();
  const geom::VertexId vc = editor.add_vertex(c, 0.0).value();
  const auto add_arc = [&](const Point3d& from, const Point3d& to, geom::VertexId v0,
                           geom::VertexId v1) -> EdgeId {
    const geom::Vec3d x = unit_vector(from - sphere.center);
    const geom::Vec3d z = unit_vector(geom::cross(x, to - sphere.center));
    const geom::Vec3d y = geom::cross(z, x);
    const double angle = std::atan2(geom::dot(to - sphere.center, y),
                                    geom::dot(to - sphere.center, x));
    const std::uint32_t curve =
        editor.add_curve(geom::CircleCurve3{sphere.center, x, y, sphere.radius});
    return editor.add_edge(v0, v1, curve, 0.0, angle, 0.0).value();
  };
  const EdgeId ab = add_arc(a, b, va, vb);
  const EdgeId bc = add_arc(b, c, vb, vc);
  const EdgeId ca = add_arc(c, a, vc, va);

  const geom::FaceId face = editor.add_face(surface, reversed, 0.0).value();
  const std::vector<geom::CoedgeSpec> loop{geom::CoedgeSpec{ab}, geom::CoedgeSpec{bc},
                                           geom::CoedgeSpec{ca}};
  editor.add_loop(face, loop, true).value();
  const geom::ShellId shell =
      editor.add_shell(std::array<geom::FaceId, 1>{face}, false).value();
  editor.add_solid(std::array<geom::ShellId, 1>{shell}, 0.0).value();
  return editor.commit().value();
}

// A single spherical face bounded by three great-circle arcs through +x, +y and
// Sphere3::eval(pi/4, pi/3). One corner is nearer the pole than the other two and
// the equator corners span exactly pi/2 of longitude, but the two side arcs are
// not meridians: the trim is a triangular lune sector, not a spherical octant.
Brep non_octant_sphere_face() {
  constexpr double kPi = 3.14159265358979323846;
  const geom::Sphere3 sphere{
      Point3d{0.0, 0.0, 0.0}, geom::Vec3d{0.0, 0.0, 1.0}, geom::Vec3d{1.0, 0.0, 0.0}, 0.25};
  return sphere_face_from_corners(sphere, Point3d{0.25, 0.0, 0.0}, Point3d{0.0, 0.25, 0.0},
                                  sphere.eval(0.25 * kPi, kPi / 3.0), false);
}
}  // namespace

TEST_CASE("rounded box rounds all twelve edges with exact corner patches") {
  const Brep box = box_2x3x4();
  const std::vector<EdgeId> edges = all_edges(box);
  REQUIRE(edges.size() == 12u);

  const auto result = geom::fillet_edges(box, edges, geom::FilletOptions{0.25, 1e-9});
  REQUIRE(result.has_value());
  const Brep& rounded = result.value().brep;
  CHECK(result.value().filleted_edges == 12u);
  CHECK(rounded.solid_count() == 1u);
  CHECK(rounded.shell_count() == 1u);
  CHECK(rounded.vertex_count() == 24u);
  CHECK(rounded.edge_count() == 48u);
  CHECK(rounded.face_count() == 26u);
  CHECK(geom::validate_topology(rounded).ok());
  CHECK(volume_of(rounded) == doctest::Approx(23.53807140332).epsilon(1e-12));

  std::uint32_t planes = 0;
  std::uint32_t cylinders = 0;
  std::uint32_t spheres = 0;
  for (const geom::FaceId id : rounded.face_ids()) {
    const geom::SurfaceGeometry& surface =
        rounded.geometry().surfaces[rounded.face(id).surface];
    if (std::holds_alternative<geom::Plane3>(surface)) ++planes;
    if (std::holds_alternative<geom::Cylinder3>(surface)) ++cylinders;
    if (std::holds_alternative<geom::Sphere3>(surface)) ++spheres;
  }
  CHECK(planes == 6u);
  CHECK(cylinders == 12u);
  CHECK(spheres == 8u);
}

TEST_CASE("rounded box tessellates watertight") {
  const Brep box = box_2x3x4();
  const auto rounded = geom::fillet_edges(box, all_edges(box), geom::FilletOptions{0.25, 1e-9});
  REQUIRE(rounded.has_value());
  const auto mesh = geom::tessellate(rounded.value().brep, geom::TessellationOptions{});
  REQUIRE(mesh.has_value());
  CHECK(mesh.value().triangles.size() > 0u);
  std::map<std::pair<std::uint32_t, std::uint32_t>, std::uint32_t> uses;
  for (const auto& triangle : mesh.value().triangles) {
    for (int e = 0; e < 3; ++e) {
      std::uint32_t a = triangle[static_cast<std::size_t>(e)];
      std::uint32_t b = triangle[static_cast<std::size_t>((e + 1) % 3)];
      if (a > b) std::swap(a, b);
      ++uses[{a, b}];
    }
  }
  std::uint32_t boundary = 0;
  for (const auto& [edge, count] : uses) {
    if (count == 1) ++boundary;
  }
  CHECK(boundary == 0u);
  // No degenerate or duplicate triangles: a patch tessellator that fans from a
  // pole can emit collapsed rows if the pole column is handled twice.
  for (const auto& triangle : mesh.value().triangles) {
    const geom::Point3d& a = mesh.value().positions[triangle[0]];
    const geom::Point3d& b = mesh.value().positions[triangle[1]];
    const geom::Point3d& c = mesh.value().positions[triangle[2]];
    const geom::Vec3d n = geom::cross(b - a, c - a);
    CHECK(geom::dot(n, n) > 0.0);
  }
  // The corner patches must wind outward: a spherical triangle whose normal
  // points back at its sphere centre would be inside out.
  const Brep& solid = rounded.value().brep;
  std::vector<const geom::Sphere3*> face_sphere(solid.face_count(), nullptr);
  for (const geom::FaceId id : solid.face_ids()) {
    const geom::SurfaceGeometry& surface = solid.geometry().surfaces[solid.face(id).surface];
    if (const geom::Sphere3* sphere = std::get_if<geom::Sphere3>(&surface)) {
      face_sphere[id.index()] = sphere;
    }
  }
  bool saw_outward_sphere_triangle = false;
  for (std::size_t t = 0; t < mesh.value().triangles.size(); ++t) {
    const geom::Sphere3* sphere = face_sphere[mesh.value().triangle_face[t]];
    if (sphere == nullptr) continue;
    const auto& triangle = mesh.value().triangles[t];
    const geom::Point3d& a = mesh.value().positions[triangle[0]];
    const geom::Point3d& b = mesh.value().positions[triangle[1]];
    const geom::Point3d& c = mesh.value().positions[triangle[2]];
    if (geom::dot(geom::cross(b - a, c - a), a - sphere->center) > 0.0) {
      saw_outward_sphere_triangle = true;
      break;
    }
  }
  CHECK(saw_outward_sphere_triangle);
}

TEST_CASE("rounded box has exact bounds and corner patches") {
  const Brep box = box_2x3x4();
  const auto rounded = geom::fillet_edges(box, all_edges(box), geom::FilletOptions{0.25, 1e-9});
  REQUIRE(rounded.has_value());
  const Brep& solid = rounded.value().brep;
  double lo[3] = {1e300, 1e300, 1e300};
  double hi[3] = {-1e300, -1e300, -1e300};
  for (const geom::VertexId id : solid.vertex_ids()) {
    const geom::Point3d& p = solid.vertex(id).point;
    for (int axis = 0; axis < 3; ++axis) {
      lo[axis] = std::min(lo[axis], p[axis]);
      hi[axis] = std::max(hi[axis], p[axis]);
    }
  }
  CHECK(lo[0] == 0.0);
  CHECK(lo[1] == 0.0);
  CHECK(lo[2] == 0.0);
  CHECK(hi[0] == 2.0);
  CHECK(hi[1] == 3.0);
  CHECK(hi[2] == 4.0);
  bool saw_hi_plane = false;
  bool saw_corner_sphere = false;
  for (const geom::SurfaceGeometry& surface : solid.geometry().surfaces) {
    if (const geom::Plane3* plane = std::get_if<geom::Plane3>(&surface)) {
      if (plane->normal[0] == 1.0 && plane->offset == 2.0) saw_hi_plane = true;
    }
    if (const geom::Sphere3* sphere = std::get_if<geom::Sphere3>(&surface)) {
      if (sphere->center[0] == 1.75 && sphere->center[1] == 2.75 && sphere->center[2] == 3.75) {
        saw_corner_sphere = true;
      }
    }
  }
  CHECK(saw_hi_plane);
  CHECK(saw_corner_sphere);
}

TEST_CASE("rounded box mass properties match the Steiner volume") {
  const Brep box = box_2x3x4();
  const auto rounded = geom::fillet_edges(box, all_edges(box), geom::FilletOptions{0.25, 1e-9});
  REQUIRE(rounded.has_value());
  const auto mass = geom::mass_properties(rounded.value().brep);
  REQUIRE(mass.has_value());
  CHECK(mass.value().volume == doctest::Approx(23.53807140332).epsilon(1e-9));
  CHECK(mass.value().centroid[0] == doctest::Approx(1.0).epsilon(1e-9));
  CHECK(mass.value().centroid[1] == doctest::Approx(1.5).epsilon(1e-9));
  CHECK(mass.value().centroid[2] == doctest::Approx(2.0).epsilon(1e-9));
}

TEST_CASE("rounded box rejects a radius at half the smallest dimension") {
  const Brep box = box_2x3x4();
  const auto edges = all_edges(box);
  const auto at_limit = geom::fillet_edges(box, edges, geom::FilletOptions{1.0, 1e-9});
  REQUIRE_FALSE(at_limit.has_value());
  CHECK(at_limit.error().code == geom::Errc::unsupported);
  CHECK(std::string_view{at_limit.error().message} ==
        "fillet_edges: radius must be less than half the smallest box dimension");
  const auto near_limit = geom::fillet_edges(box, edges, geom::FilletOptions{1.0 - 1e-6, 1e-9});
  REQUIRE(near_limit.has_value());
  CHECK(near_limit.value().filleted_edges == 12u);
  CHECK(near_limit.value().brep.face_count() == 26u);
  CHECK(geom::validate_topology(near_limit.value().brep).ok());
}

TEST_CASE("rounded box rejects an incomplete or dead selection") {
  const Brep box = box_2x3x4();
  auto edges = all_edges(box);
  edges.pop_back();
  const auto eleven = geom::fillet_edges(box, edges, geom::FilletOptions{0.25, 1e-9});
  CHECK_FALSE(eleven.has_value());
  edges = all_edges(box);
  edges.push_back(edges.front());
  const auto duplicated = geom::fillet_edges(box, edges, geom::FilletOptions{0.25, 1e-9});
  REQUIRE(duplicated.has_value());
  CHECK(duplicated.value().filleted_edges == 12u);
  edges = all_edges(box);
  edges.push_back(geom::EdgeId::make(999u, 0u));
  const auto dead = geom::fillet_edges(box, edges, geom::FilletOptions{0.25, 1e-9});
  REQUIRE_FALSE(dead.has_value());
  CHECK(dead.error().code == geom::Errc::index_out_of_range);
  CHECK(std::string_view{dead.error().message} == "fillet_edges: edge does not belong to the solid");
  CHECK(dead.error().context_a == 999u);
  // A dead edge is the caller's mistake even when the rest of the selection is
  // incomplete: it must win over the family/parallel diagnostics.
  edges = all_edges(box);
  edges.pop_back();
  edges.push_back(geom::EdgeId::make(999u, 0u));
  const auto incomplete_dead = geom::fillet_edges(box, edges, geom::FilletOptions{0.25, 1e-9});
  REQUIRE_FALSE(incomplete_dead.has_value());
  CHECK(incomplete_dead.error().code == geom::Errc::index_out_of_range);
  CHECK(std::string_view{incomplete_dead.error().message} ==
        "fillet_edges: edge does not belong to the solid");
  CHECK(incomplete_dead.error().context_a == 999u);
}

TEST_CASE("sphere patch predicate matches only the eight corner faces") {
  const Brep box = box_2x3x4();
  const auto rounded = geom::fillet_edges(box, all_edges(box), geom::FilletOptions{0.25, 1e-9});
  REQUIRE(rounded.has_value());
  const Brep& solid = rounded.value().brep;
  std::size_t patches = 0;
  std::size_t others = 0;
  for (const geom::FaceId id : solid.face_ids()) {
    const geom::Face& face = solid.face(id);
    const bool sphere = face.surface < solid.geometry().surfaces.size() &&
                        std::holds_alternative<geom::Sphere3>(solid.geometry().surfaces[face.surface]);
    if (sphere) {
      CHECK(geom::detail::tessellate_face_is_sphere_patch(solid, id));
      ++patches;
    } else {
      CHECK_FALSE(geom::detail::tessellate_face_is_sphere_patch(solid, id));
      ++others;
    }
  }
  CHECK(patches == 8u);
  CHECK(others == 18u);
}

TEST_CASE("rounded box rejects a non-box shape") {
  const Brep cylinder =
      geom::make_cylinder(Point3d{0.0, 0.0, 0.0}, geom::Vec3d{0.0, 0.0, 1.0}, 1.0, 2.0, 1e-9)
          .value();
  const auto edges = all_edges(cylinder);
  const auto result = geom::fillet_edges(cylinder, edges, geom::FilletOptions{0.25, 1e-9});
  REQUIRE_FALSE(result.has_value());
  CHECK(result.error().code == geom::Errc::unsupported);
  CHECK(std::string_view{result.error().message} ==
        "fillet_edges: only a box or a convex prism is supported");
  // The strict box recogniser owns the box-only diagnostic; the public entry
  // point folds it into the box-or-prism rejection above.
  const auto bounds = geom::detail::fillet_box_bounds(cylinder, 1e-9);
  REQUIRE_FALSE(bounds.has_value());
  CHECK(bounds.error().code == geom::Errc::unsupported);
  CHECK(std::string_view{bounds.error().message} ==
        "fillet_edges: only an axis-aligned box is supported");
}

TEST_CASE("rounded box keeps the long axis of a non-cubic box") {
  const Brep box = box_2x3x4();
  const auto rounded = geom::fillet_edges(box, all_edges(box), geom::FilletOptions{0.25, 1e-9});
  REQUIRE(rounded.has_value());
  const Brep& solid = rounded.value().brep;
  CHECK(volume_of(solid) == doctest::Approx(23.53807140332).epsilon(1e-9));
  double lo[3] = {1e300, 1e300, 1e300};
  double hi[3] = {-1e300, -1e300, -1e300};
  for (const geom::VertexId id : solid.vertex_ids()) {
    const geom::Point3d& p = solid.vertex(id).point;
    for (int axis = 0; axis < 3; ++axis) {
      lo[axis] = std::min(lo[axis], p[axis]);
      hi[axis] = std::max(hi[axis], p[axis]);
    }
  }
  CHECK(lo[0] == 0.0);
  CHECK(lo[1] == 0.0);
  CHECK(lo[2] == 0.0);
  CHECK(hi[0] == 2.0);
  CHECK(hi[1] == 3.0);
  CHECK(hi[2] == 4.0);
}

TEST_CASE("sphere patch rejects a non-octant trim") {
  const Brep face = non_octant_sphere_face();
  const auto mesh = geom::tessellate(face, geom::TessellationOptions{});
  REQUIRE_FALSE(mesh.has_value());
  CHECK(mesh.error().code == geom::Errc::unsupported);
  CHECK(std::string_view{mesh.error().message} == "tessellate: sphere patch is not a single octant");
}

// A reversed octant winds inward if the patch builder trusts the uv order alone,
// so this case fails if the outward-orientation guard is removed.
TEST_CASE("sphere patch orients a reversed octant outward") {
  const geom::Sphere3 sphere{
      Point3d{0.0, 0.0, 0.0}, geom::Vec3d{0.0, 0.0, 1.0}, geom::Vec3d{1.0, 0.0, 0.0}, 0.25};
  const Brep face = sphere_face_from_corners(sphere, Point3d{0.25, 0.0, 0.0},
                                             Point3d{0.0, 0.25, 0.0}, Point3d{0.0, 0.0, 0.25},
                                             true);
  const auto mesh = geom::tessellate(face, geom::TessellationOptions{});
  REQUIRE(mesh.has_value());
  REQUIRE_FALSE(mesh.value().triangles.empty());
  std::size_t inward = 0;
  for (const auto& triangle : mesh.value().triangles) {
    const geom::Point3d& a = mesh.value().positions[triangle[0]];
    const geom::Point3d& b = mesh.value().positions[triangle[1]];
    const geom::Point3d& c = mesh.value().positions[triangle[2]];
    if (geom::dot(geom::cross(b - a, c - a), a - sphere.center) <= 0.0) ++inward;
  }
  CHECK(inward == 0u);
  CHECK(mesh.value().triangles.size() == 153u);
}
