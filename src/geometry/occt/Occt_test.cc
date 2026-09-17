#include <BRepAlgoAPI_Fuse.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepPrimAPI_MakeSphere.hxx>
#include <ShapeUpgrade_UnifySameDomain.hxx>
#include <TopExp_Explorer.hxx>
#include <catch2/catch_all.hpp>
#include <cmath>
#include <string>
#include <vector>

#include "geometry/occt/OcctBoolean.h"
#include "geometry/occt/OcctBridge.h"
#include "geometry/occt/OcctHull.h"
#include "geometry/occt/OcctMesh.h"
#include "geometry/occt/OcctRevolvedHull.h"

using Catch::Matchers::WithinRel;

namespace {

TopoDS_Shape box(double x, double y, double z, double dx = 10, double dy = 10, double dz = 10)
{
  return BRepPrimAPI_MakeBox(gp_Pnt(x, y, z), dx, dy, dz).Shape();
}

TopoDS_Shape sphere(double x, double y, double z, double r)
{
  return BRepPrimAPI_MakeSphere(gp_Pnt(x, y, z), r).Shape();
}

TopoDS_Shape cylinder(double x, double y, double z, double r, double h)
{
  return BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(x, y, z), gp_Dir(0, 0, 1)), r, h).Shape();
}

double volume(const TopoDS_Shape& s)
{
  return OcctBoolean::extent(s, 3);
}
size_t solids(const TopoDS_Shape& s)
{
  return OcctBoolean::piecesOf(s, 3).size();
}
size_t faces(const TopoDS_Shape& s)
{
  size_t n = 0;
  for (TopExp_Explorer it(s, TopAbs_FACE); it.More(); it.Next()) ++n;
  return n;
}

// OpenSCAD's polyhedron convention: faces wound clockwise seen from outside.
OcctBridge::MeshData cubeMesh()
{
  OcctBridge::MeshData m;
  m.vertices = {{0, 0, 0},  {10, 0, 0},  {10, 10, 0},  {0, 10, 0},
                {0, 0, 10}, {10, 0, 10}, {10, 10, 10}, {0, 10, 10}};
  m.faces = {{0, 1, 2, 3}, {4, 7, 6, 5}, {0, 4, 5, 1}, {1, 5, 6, 2}, {2, 6, 7, 3}, {3, 7, 4, 0}};
  return m;
}

}  // namespace

TEST_CASE("OcctBoolean fuse", "[occt][boolean]")
{
  SECTION("overlapping boxes merge into one solid with the shared volume counted once")
  {
    auto f = OcctBoolean::fuse({box(0, 0, 0), box(5, 5, 5)}, 3);
    CHECK(solids(f) == 1);
    CHECK_THAT(volume(f), WithinRel(2000.0 - 125.0, 1e-9));
  }
  SECTION("disjoint boxes stay two solids")
  {
    auto f = OcctBoolean::fuse({box(0, 0, 0), box(20, 0, 0)}, 3);
    CHECK(solids(f) == 2);
    CHECK_THAT(volume(f), WithinRel(2000.0, 1e-9));
  }
  SECTION("boxes sharing a face become one solid")
  {
    auto f = OcctBoolean::fuse({box(0, 0, 0), box(10, 0, 0)}, 3);
    CHECK(solids(f) == 1);
    CHECK_THAT(volume(f), WithinRel(2000.0, 1e-9));
  }
  SECTION("N-ary: twenty stacked cylinders")
  {
    std::vector<TopoDS_Shape> operands;
    for (int i = 0; i < 20; ++i) operands.push_back(cylinder(i * 0.5, 0, 0, 3, 5));
    auto f = OcctBoolean::fuse(operands, 3);
    CHECK(solids(f) == 1);
    CHECK(volume(f) > M_PI * 9 * 5);
    CHECK(volume(f) < 20 * M_PI * 9 * 5);
  }
}

TEST_CASE("OcctBoolean cut and common", "[occt][boolean]")
{
  const auto cube = box(-5, -5, -5);
  const auto ball = sphere(0, 0, 0, 6);
  const auto cut = OcctBoolean::cut({cube}, {ball}, 3);
  const auto common = OcctBoolean::common({cube}, {ball}, 3);
  SECTION("cut plus common is the whole")
  {
    CHECK_THAT(volume(cut) + volume(common), WithinRel(1000.0, 1e-7));
  }
  SECTION("cutting along a coincident surface leaves the body alone")
  {
    // The remainder abuts the sphere exactly; OCCT has answered this with
    // nothing. The result must be the remainder, not empty.
    auto again = OcctBoolean::cut({cut}, {ball}, 3);
    REQUIRE(!again.IsNull());
    CHECK_THAT(volume(again), WithinRel(volume(cut), 1e-6));
  }
  SECTION("common of disjoint shapes is empty")
  {
    CHECK(OcctBoolean::common({box(0, 0, 0)}, {box(20, 0, 0)}, 3).IsNull());
  }
  SECTION("cut that removes everything is empty")
  {
    CHECK(OcctBoolean::cut({box(0, 0, 0, 2, 2, 2)}, {box(-5, -5, -5, 20, 20, 20)}, 3).IsNull());
  }
}

TEST_CASE("OcctBoolean unify keeps its hands off a seam-crossed face", "[occt][boolean]")
{
  // UnifySameDomain drops the spherical cap crossed by the sphere's seam
  // here and returns an invalid solid; the guarded unify must notice.
  const auto fused = OcctBoolean::fuse({box(-5, -5, -5)}, 3);
  const auto raw = [] {
    BRepAlgoAPI_Fuse op(box(-5, -5, -5), sphere(0, 0, 0, 6));
    return op.Shape();
  }();
  ShapeUpgrade_UnifySameDomain unifier(raw, true, true, true);
  unifier.Build();
  const double unified = volume(unifier.Shape());
  const double before = volume(raw);
  INFO("raw " << before << " after UnifySameDomain " << unified);
  CHECK_THAT(volume(OcctBoolean::unify(raw, 3)), WithinRel(before, 1e-9));
  (void)fused;
}

TEST_CASE("OcctMesh builds solids from OpenSCAD-style faces", "[occt][mesh]")
{
  SECTION("a cube from six quads")
  {
    const auto m = cubeMesh();
    bool nonPlanar = false;
    int freeEdges = 0;
    auto s = OcctMesh::solidsFromFaces(m.vertices, m.faces, nonPlanar, freeEdges);
    REQUIRE(!s.IsNull());
    CHECK(!nonPlanar);
    CHECK(freeEdges == 0);
    CHECK_THAT(volume(s), WithinRel(1000.0, 1e-9));
    CHECK(faces(s) == 6);
  }
  SECTION("a missing face is reported as free edges, not a solid")
  {
    auto m = cubeMesh();
    m.faces.pop_back();
    bool nonPlanar = false;
    int freeEdges = 0;
    auto s = OcctMesh::solidsFromFaces(m.vertices, m.faces, nonPlanar, freeEdges);
    CHECK(s.IsNull());
    CHECK(freeEdges > 0);
  }
  SECTION("a non-planar face is reported")
  {
    auto m = cubeMesh();
    m.vertices[6] = {10, 10, 12};  // lifts one corner of three faces
    bool nonPlanar = false;
    int freeEdges = 0;
    auto s = OcctMesh::solidsFromFaces(m.vertices, m.faces, nonPlanar, freeEdges);
    CHECK(s.IsNull());
    CHECK(nonPlanar);
  }
}

TEST_CASE("OcctMesh repairs collinear zero-area triangles", "[occt][mesh]")
{
  // A triangulated cube where one triangle's edge carries a midpoint
  // vertex on one side only, filled by a zero-area triangle: the pattern
  // Manifold's hull emits. Index-manifold, yet the degenerate triangle
  // cannot become a face.
  OcctBridge::MeshData m;
  m.vertices = {{0, 0, 0},   {10, 0, 0},   {10, 10, 0}, {0, 10, 0}, {0, 0, 10},
                {10, 0, 10}, {10, 10, 10}, {0, 10, 10}, {5, 0, 0}};  // 8: midpoint of edge 0-1
  m.faces = {
    {0, 1, 2}, {0, 2, 3},             // bottom (clockwise from outside)
    {4, 7, 6}, {4, 6, 5},             // top
    {0, 4, 5}, {0, 5, 8}, {8, 5, 1},  // front, with the split edge 0-8-1
    {1, 0, 8},                        // zero-area filler along 0-1
    {1, 5, 6}, {1, 6, 2},             // right
    {2, 6, 7}, {2, 7, 3},             // back
    {3, 7, 4}, {3, 4, 0},             // left
  };
  std::string why;
  auto s = OcctMesh::solidsFromMesh(m, why);
  INFO(why);
  REQUIRE(!s.IsNull());
  CHECK_THAT(volume(s), WithinRel(1000.0, 1e-9));
  CHECK(faces(s) == 6);
}

TEST_CASE("OcctHull closed forms", "[occt][hull]")
{
  std::string rung;
  SECTION("two equal spheres on a line: a capsule")
  {
    auto h = OcctHull::hull({sphere(0, 0, 0, 3), sphere(20, 0, 0, 3)}, 3, rung);
    REQUIRE(!h.IsNull());
    CHECK(rung == "hull of equal spheres");
    CHECK_THAT(volume(h), WithinRel(M_PI * 9 * 20 + 4.0 / 3.0 * M_PI * 27, 1e-9));
  }
  SECTION("two unequal spheres: caps and a tangent cone")
  {
    const double ra = 4, rb = 9, d = std::sqrt(20.0 * 20 + 5 * 5 + 3 * 3);
    auto h = OcctHull::hull({sphere(0, 0, 0, ra), sphere(20, 5, 3, rb)}, 3, rung);
    REQUIRE(!h.IsNull());
    CHECK(rung == "hull of two spheres");
    const double sinA = (rb - ra) / d, cosA = std::sqrt(1 - sinA * sinA);
    const double h1 = ra * (1 - sinA), h2 = rb * (1 + sinA), length = d * cosA * cosA;
    const double rho1 = ra * cosA, rho2 = rb * cosA;
    const double exact = M_PI * h1 * h1 * (3 * ra - h1) / 3 +
                         M_PI * length / 3 * (rho1 * rho1 + rho1 * rho2 + rho2 * rho2) +
                         M_PI * h2 * h2 * (3 * rb - h2) / 3;
    CHECK_THAT(volume(h), WithinRel(exact, 1e-6));
    CHECK(faces(h) == 3);
  }
  SECTION("four corner posts: a rounded box")
  {
    auto h = OcctHull::hull({cylinder(0, 0, 0, 3, 8), cylinder(30, 0, 0, 3, 8), cylinder(0, 20, 0, 3, 8),
                             cylinder(30, 20, 0, 3, 8)},
                            3, rung);
    REQUIRE(!h.IsNull());
    CHECK(rung == "hull of parallel cylinders");
    CHECK_THAT(volume(h), WithinRel((30 * 20 + 2 * 3 * (30 + 20) + M_PI * 9) * 8, 1e-9));
  }
  SECTION("two boxes: the convex hull of their vertices")
  {
    auto h = OcctHull::hull({box(0, 0, 0), box(0, 0, 20)}, 3, rung);
    REQUIRE(!h.IsNull());
    CHECK(rung == "hull of polyhedra");
    CHECK_THAT(volume(h), WithinRel(3000.0, 1e-9));
    CHECK(faces(h) == 6);
  }
  SECTION("a sphere merely translated is exactly an offset hull")
  {
    auto h = OcctHull::hull(
      {sphere(0, 0, 0, 3), sphere(20, 0, 0, 3), sphere(0, 15, 0, 3), sphere(20, 15, 0, 3)}, 3, rung);
    REQUIRE(!h.IsNull());
    // Coplanar centers: the equal-spheres rung declines and the revolved
    // rung builds the rounded coin.
    CHECK(rung == "hull of translated revolution solids");
    const double A = 300, P = 70, r = 3;
    CHECK_THAT(volume(h),
               WithinRel(2 * r * A + M_PI * r * r / 2 * P + 4.0 / 3.0 * M_PI * r * r * r, 1e-9));
  }
}

TEST_CASE("OcctHull minkowski with a ball is an offset", "[occt][hull]")
{
  std::string rung;
  auto m = OcctHull::minkowski({box(0, 0, 0), sphere(0, 0, 0, 2)}, 3, rung);
  REQUIRE(!m.IsNull());
  CHECK(rung == "minkowski with a ball");
  // Steiner: V + A r + (r^2/2) sum(L theta) + 4/3 pi r^3, for a cube with
  // 12 edges of length 10 and exterior angle pi/2.
  const double r = 2;
  const double exact = 1000 + 600 * r + r * r / 2 * 12 * 10 * M_PI / 2 + 4.0 / 3.0 * M_PI * r * r * r;
  CHECK_THAT(volume(m), WithinRel(exact, 1e-9));
}
