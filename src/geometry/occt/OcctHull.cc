#include "geometry/occt/OcctHull.h"

#include <BRepAdaptor_Curve.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepBuilderAPI_Copy.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakePolygon.hxx>
#include <BRepBuilderAPI_MakeSolid.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepBuilderAPI_Sewing.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepBndLib.hxx>
#include <BRepGProp.hxx>
#include <BRepMesh_IncrementalMesh.hxx>
#include <Bnd_Box.hxx>
#include <Poly_Triangulation.hxx>
#include <TopLoc_Location.hxx>
#include <BRepOffsetAPI_MakeOffset.hxx>
#include <BRepOffsetAPI_MakeOffsetShape.hxx>
#include <BRepPrimAPI_MakeCone.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepPrimAPI_MakePrism.hxx>
#include <BRepPrimAPI_MakeSphere.hxx>
#include <BRep_Tool.hxx>
#include <GC_MakeArcOfCircle.hxx>
#include <GProp_GProps.hxx>
#include <GeomAbs_JoinType.hxx>
#include <GeomAbs_SurfaceType.hxx>
#include <Geom_TrimmedCurve.hxx>
#include <Standard_Failure.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Shell.hxx>
#include <TopoDS_Solid.hxx>
#include <TopoDS_Vertex.hxx>
#include <TopoDS_Wire.hxx>
#include <gp_Ax2.hxx>
#include <gp_Ax3.hxx>
#include <gp_Circ.hxx>
#include <gp_Pnt.hxx>
#include <gp_Trsf.hxx>
#include <gp_Vec.hxx>
#include <algorithm>
#include <array>
#include <cmath>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "geometry/occt/OcctBoolean.h"
#include "geometry/occt/OcctMesh.h"
#include "geometry/occt/OcctRevolvedHull.h"

#ifdef ENABLE_MANIFOLD
#include <manifold/manifold.h>
#endif

namespace OcctHull {

namespace {

constexpr double RADIUS_REL_TOL = 1e-6;
constexpr double COLLINEAR_REL_TOL = 1e-6;
constexpr double PARALLEL_TOL = 1e-6;
constexpr double SPAN_TOL = 1e-6;
constexpr double ORIGIN_TOL = 1e-6;
constexpr double OFFSET_TOL = 1e-4;
constexpr double SELF_CHECK_RTOL = 1e-6;
constexpr size_t MIN_MESH_BALL_VERTICES = 24;
constexpr double MESH_BALL_REL_TOL = 1e-3;

using P3 = std::array<double, 3>;
using P2 = std::array<double, 2>;

double dist(const P3& a, const P3& b)
{
  return std::sqrt((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) +
                   (a[2] - b[2]) * (a[2] - b[2]));
}

double dist2(const P2& a, const P2& b)
{
  return std::hypot(a[0] - b[0], a[1] - b[1]);
}

bool isValid(const TopoDS_Shape& shape)
{
  if (shape.IsNull()) return false;
  try {
    return BRepCheck_Analyzer(shape).IsValid();
  } catch (const Standard_Failure&) {
    return false;
  }
}

double volume(const TopoDS_Shape& shape)
{
  return OcctBoolean::extent(shape, 3);
}
double area(const TopoDS_Shape& shape)
{
  return OcctBoolean::extent(shape, 2);
}

std::vector<TopoDS_Shape> facesOf(const TopoDS_Shape& shape)
{
  std::vector<TopoDS_Shape> faces;
  for (TopExp_Explorer it(shape, TopAbs_FACE); it.More(); it.Next()) faces.push_back(it.Current());
  return faces;
}

std::vector<TopoDS_Shape> edgesOf(const TopoDS_Shape& shape)
{
  std::vector<TopoDS_Shape> edges;
  for (TopExp_Explorer it(shape, TopAbs_EDGE); it.More(); it.Next()) edges.push_back(it.Current());
  return edges;
}

size_t solidCount(const TopoDS_Shape& shape)
{
  size_t n = 0;
  for (TopExp_Explorer it(shape, TopAbs_SOLID); it.More(); it.Next()) ++n;
  return n;
}

std::vector<P3> verticesOf(const TopoDS_Shape& shape)
{
  std::vector<P3> points;
  for (TopExp_Explorer it(shape, TopAbs_VERTEX); it.More(); it.Next()) {
    const auto p = BRep_Tool::Pnt(TopoDS::Vertex(it.Current()));
    const P3 q{p.X(), p.Y(), p.Z()};
    bool seen = false;
    for (const auto& existing : points) {
      if (dist(existing, q) < 1e-9) {
        seen = true;
        break;
      }
    }
    if (!seen) points.push_back(q);
  }
  return points;
}

GeomAbs_SurfaceType surfaceType(const TopoDS_Shape& face)
{
  return BRepAdaptor_Surface(TopoDS::Face(face)).GetType();
}

P3 faceCenter(const TopoDS_Shape& face)
{
  GProp_GProps props;
  BRepGProp::SurfaceProperties(face, props);
  const auto c = props.CentreOfMass();
  return {c.X(), c.Y(), c.Z()};
}

// --- classifiers -----------------------------------------------------------

struct SphereInfo {
  P3 center;
  double radius;
};

std::optional<SphereInfo> asSphere(const TopoDS_Shape& shape)
{
  const auto faces = facesOf(shape);
  if (solidCount(shape) != 1 || faces.size() != 1) return std::nullopt;
  BRepAdaptor_Surface surface(TopoDS::Face(faces.front()));
  if (surface.GetType() != GeomAbs_Sphere) return std::nullopt;
  const auto sphere = surface.Sphere();
  const auto c = sphere.Location();
  return SphereInfo{{c.X(), c.Y(), c.Z()}, sphere.Radius()};
}

struct CylinderInfo {
  P3 a, b;  // cap centers
  double radius;
};

// A plain cylinder: one cylindrical face between two planar caps. A cone
// (r1 != r2) has a conical face and is excluded by the type check.
std::optional<CylinderInfo> asCylinder(const TopoDS_Shape& shape)
{
  const auto faces = facesOf(shape);
  if (solidCount(shape) != 1 || faces.size() != 3) return std::nullopt;
  std::vector<TopoDS_Shape> cylinders, planes;
  for (const auto& f : faces) {
    const auto type = surfaceType(f);
    if (type == GeomAbs_Cylinder) cylinders.push_back(f);
    else if (type == GeomAbs_Plane) planes.push_back(f);
  }
  if (cylinders.size() != 1 || planes.size() != 2) return std::nullopt;
  BRepAdaptor_Surface surface(TopoDS::Face(cylinders.front()));
  return CylinderInfo{faceCenter(planes[0]), faceCenter(planes[1]), surface.Cylinder().Radius()};
}

struct CircleInfo {
  P2 center;
  double radius;
};

std::optional<CircleInfo> asDisc(const TopoDS_Shape& shape)
{
  if (solidCount(shape) != 0) return std::nullopt;
  const auto faces = facesOf(shape);
  if (faces.size() != 1) return std::nullopt;
  const auto edges = edgesOf(faces.front());
  if (edges.size() != 1) return std::nullopt;
  BRepAdaptor_Curve curve(TopoDS::Edge(edges.front()));
  if (curve.GetType() != GeomAbs_Circle) return std::nullopt;
  const auto circle = curve.Circle();
  const auto c = circle.Location();
  if (std::fabs(c.Z()) > 1e-9) return std::nullopt;
  return CircleInfo{{c.X(), c.Y()}, circle.Radius()};
}

// --- convex hulls ----------------------------------------------------------

double cross2(const P2& o, const P2& a, const P2& b)
{
  return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0]);
}

// Andrew's monotone chain; counterclockwise, no collinear points.
std::vector<P2> convexHull2d(std::vector<P2> points)
{
  std::sort(points.begin(), points.end());
  points.erase(std::unique(points.begin(), points.end(),
                           [](const P2& a, const P2& b) { return dist2(a, b) < 1e-9; }),
               points.end());
  if (points.size() < 3) return points;
  std::vector<P2> hull(2 * points.size());
  size_t k = 0;
  for (const auto& p : points) {
    while (k >= 2 && cross2(hull[k - 2], hull[k - 1], p) <= 1e-12) --k;
    hull[k++] = p;
  }
  for (size_t i = points.size() - 1, t = k + 1; i-- > 0;) {
    while (k >= t && cross2(hull[k - 2], hull[k - 1], points[i]) <= 1e-12) --k;
    hull[k++] = points[i];
  }
  hull.resize(k - 1);
  return hull;
}

// The convex hull of a 3D point set as a solid with merged planar faces.
TopoDS_Shape convexHull3d(const std::vector<P3>& points)
{
#ifdef ENABLE_MANIFOLD
  std::vector<manifold::vec3> pts;
  pts.reserve(points.size());
  for (const auto& p : points) pts.emplace_back(p[0], p[1], p[2]);
  manifold::Manifold hull = manifold::Manifold::Hull(pts);
  if (hull.IsEmpty()) return {};
  const auto mesh = hull.GetMeshGL64();
  OcctBridge::MeshData data;
  const size_t numVert = mesh.vertProperties.size() / mesh.numProp;
  data.vertices.reserve(numVert);
  for (size_t i = 0; i < numVert; ++i) {
    data.vertices.emplace_back(mesh.vertProperties[i * mesh.numProp],
                               mesh.vertProperties[i * mesh.numProp + 1],
                               mesh.vertProperties[i * mesh.numProp + 2]);
  }
  for (size_t i = 0; i + 2 < mesh.triVerts.size(); i += 3) {
    // Manifold winds triangles counterclockwise from outside; the mesh
    // builder expects OpenSCAD's clockwise convention.
    data.faces.push_back({static_cast<size_t>(mesh.triVerts[i + 2]),
                          static_cast<size_t>(mesh.triVerts[i + 1]),
                          static_cast<size_t>(mesh.triVerts[i])});
  }
  bool nonPlanar = false;
  int freeEdges = 0;
  return OcctMesh::solidsFromFaces(data.vertices, data.faces, nonPlanar, freeEdges);
#else
  (void)points;
  return {};
#endif
}

// --- constructions ---------------------------------------------------------

TopoDS_Shape offset3d(const TopoDS_Shape& shape, double r)
{
  try {
    BRepOffsetAPI_MakeOffsetShape op;
    op.PerformByJoin(shape, r, OFFSET_TOL, BRepOffset_Skin, true, false, GeomAbs_Arc, true);
    if (!op.IsDone()) return {};
    auto result = op.Shape();
    if (!isValid(result)) return {};
    if (volume(result) <= volume(shape)) return {};
    return result;
  } catch (const Standard_Failure&) {
    return {};
  }
}

TopoDS_Shape faceFromWire(const TopoDS_Wire& wire)
{
  BRepBuilderAPI_MakeFace maker(wire, true);
  if (!maker.IsDone()) return {};
  return maker.Face();
}

TopoDS_Shape offset2d(const TopoDS_Shape& face, double r)
{
  try {
    TopoDS_Wire outer;
    int wires = 0;
    for (TopExp_Explorer it(face, TopAbs_WIRE); it.More(); it.Next()) {
      outer = TopoDS::Wire(it.Current());
      ++wires;
    }
    if (wires != 1) return {};
    for (const double sign : {1.0, -1.0}) {
      BRepOffsetAPI_MakeOffset op;
      op.Init(GeomAbs_Arc);
      op.AddWire(outer);
      op.Perform(sign * r);
      if (!op.IsDone()) continue;
      TopoDS_Wire wire;
      for (TopExp_Explorer it(op.Shape(), TopAbs_WIRE); it.More(); it.Next()) {
        wire = TopoDS::Wire(it.Current());
      }
      if (wire.IsNull()) continue;
      auto result = faceFromWire(wire);
      if (!isValid(result)) continue;
      if (area(result) > area(face)) return result;
    }
    return {};
  } catch (const Standard_Failure&) {
    return {};
  }
}

TopoDS_Shape sphereAt(const P3& c, double r)
{
  return BRepPrimAPI_MakeSphere(gp_Pnt(c[0], c[1], c[2]), r).Shape();
}

TopoDS_Shape circleAt(const P2& c, double r)
{
  const gp_Circ circ(gp_Ax2(gp_Pnt(c[0], c[1], 0), gp_Dir(0, 0, 1)), r);
  return faceFromWire(BRepBuilderAPI_MakeWire(BRepBuilderAPI_MakeEdge(circ).Edge()).Wire());
}

TopoDS_Shape capsule(const P3& a, const P3& b, double r)
{
  const double length = dist(a, b);
  if (length < 1e-9) return sphereAt(a, r);
  const gp_Dir dir((b[0] - a[0]) / length, (b[1] - a[1]) / length, (b[2] - a[2]) / length);
  auto cylinder = BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(a[0], a[1], a[2]), dir), r, length).Shape();
  return OcctBoolean::fuse({cylinder, sphereAt(a, r), sphereAt(b, r)}, 3);
}

// A 2D capsule as one closed wire (two lines, two tangent arcs), so the
// extruded solid has four lateral faces rather than the twelve a fused
// rectangle-plus-circles leaves.
TopoDS_Shape stadium(const P2& a, const P2& b, double r)
{
  const double length = dist2(a, b);
  if (length < 1e-9) return circleAt(a, r);
  const double ux = (b[0] - a[0]) / length, uy = (b[1] - a[1]) / length;
  const double nx = -uy, ny = ux;
  const gp_Pnt c1(a[0] + nx * r, a[1] + ny * r, 0), c2(b[0] + nx * r, b[1] + ny * r, 0);
  const gp_Pnt c3(b[0] - nx * r, b[1] - ny * r, 0), c4(a[0] - nx * r, a[1] - ny * r, 0);
  try {
    BRepBuilderAPI_MakeWire wire;
    wire.Add(BRepBuilderAPI_MakeEdge(c1, c2).Edge());
    wire.Add(BRepBuilderAPI_MakeEdge(GC_MakeArcOfCircle(c2, gp_Vec(ux, uy, 0), c3).Value()).Edge());
    wire.Add(BRepBuilderAPI_MakeEdge(c3, c4).Edge());
    wire.Add(BRepBuilderAPI_MakeEdge(GC_MakeArcOfCircle(c4, gp_Vec(-ux, -uy, 0), c1).Value()).Edge());
    if (!wire.IsDone()) return {};
    return faceFromWire(wire.Wire());
  } catch (const Standard_Failure&) {
    return {};
  }
}

struct Extremes {
  double span;
  P3 a, b;
};

Extremes pairwiseExtremes(const std::vector<P3>& points)
{
  Extremes best{0.0, points.front(), points.front()};
  for (size_t i = 0; i < points.size(); ++i) {
    for (size_t j = i + 1; j < points.size(); ++j) {
      const double d = dist(points[i], points[j]);
      if (d > best.span) best = {d, points[i], points[j]};
    }
  }
  return best;
}

bool collinear(const std::vector<P3>& points, const Extremes& ext)
{
  if (ext.span < 1e-12) return true;
  P3 dir{(ext.b[0] - ext.a[0]) / ext.span, (ext.b[1] - ext.a[1]) / ext.span,
         (ext.b[2] - ext.a[2]) / ext.span};
  for (const auto& p : points) {
    P3 rel{p[0] - ext.a[0], p[1] - ext.a[1], p[2] - ext.a[2]};
    const double t = rel[0] * dir[0] + rel[1] * dir[1] + rel[2] * dir[2];
    P3 foot{ext.a[0] + t * dir[0], ext.a[1] + t * dir[1], ext.a[2] + t * dir[2]};
    if (dist(p, foot) > COLLINEAR_REL_TOL * ext.span) return false;
  }
  return true;
}

TopoDS_Shape hull3dOffset(const std::vector<P3>& points, double r)
{
  const auto ext = pairwiseExtremes(points);
  if (ext.span < 1e-9) return sphereAt(points.front(), r);
  if (collinear(points, ext)) return capsule(ext.a, ext.b, r);
  auto polytope = convexHull3d(points);
  if (polytope.IsNull()) return {};
  return offset3d(polytope, r);
}

TopoDS_Shape polygonFace(const std::vector<P2>& ring)
{
  if (ring.size() < 3) return {};
  BRepBuilderAPI_MakePolygon polygon;
  for (const auto& p : ring) polygon.Add(gp_Pnt(p[0], p[1], 0));
  polygon.Close();
  if (!polygon.IsDone()) return {};
  return faceFromWire(polygon.Wire());
}

TopoDS_Shape hull2dOffset(const std::vector<P2>& points, double r)
{
  std::vector<P3> lifted;
  for (const auto& p : points) lifted.push_back({p[0], p[1], 0.0});
  const auto ext = pairwiseExtremes(lifted);
  if (ext.span < 1e-9) return circleAt(points.front(), r);
  if (collinear(lifted, ext)) return stadium({ext.a[0], ext.a[1]}, {ext.b[0], ext.b[1]}, r);
  auto hull = convexHull2d(points);
  auto face = polygonFace(hull);
  if (face.IsNull()) return {};
  return offset2d(face, r);
}

// Exact hull of two spheres with ra < rb, neither containing the other:
// two spherical caps sewn to the external tangent cone. With d the
// center distance and sin(a) = (rb - ra)/d, the cone grazes each sphere
// at latitude -a (a circle of radius r*cos(a) at axial offset -r*sin(a)
// toward the small side), exact for any non-contained spacing.
TopoDS_Shape twoSphereHull(const P3& ca, double ra, const P3& cb, double rb)
{
  const double d = dist(ca, cb);
  const double sinA = (rb - ra) / d;
  const double cosA = std::sqrt(1.0 - sinA * sinA);
  const double alpha = std::asin(sinA);
  const gp_Dir u((cb[0] - ca[0]) / d, (cb[1] - ca[1]) / d, (cb[2] - ca[2]) / d);

  auto sphereCap = [&](const P3& c, double r, double latLo, double latHi) -> TopoDS_Shape {
    BRepPrimAPI_MakeSphere maker(gp_Ax2(gp_Pnt(c[0], c[1], c[2]), u), r, latLo, latHi);
    for (const auto& f : facesOf(maker.Shape())) {
      if (surfaceType(f) == GeomAbs_Sphere) return f;
    }
    return {};
  };

  const double x1 = -ra * sinA;
  const double x2 = d - rb * sinA;
  const double length = x2 - x1;
  if (length <= 0) return {};
  const gp_Pnt base(ca[0] + u.X() * x1, ca[1] + u.Y() * x1, ca[2] + u.Z() * x1);
  TopoDS_Shape coneFace;
  try {
    BRepPrimAPI_MakeCone cone(gp_Ax2(base, u), ra * cosA, rb * cosA, length);
    for (const auto& f : facesOf(cone.Shape())) {
      if (surfaceType(f) == GeomAbs_Cone) coneFace = f;
    }
  } catch (const Standard_Failure&) {
    return {};
  }
  const auto capA = sphereCap(ca, ra, -M_PI / 2, -alpha);
  const auto capB = sphereCap(cb, rb, -alpha, M_PI / 2);
  if (capA.IsNull() || capB.IsNull() || coneFace.IsNull()) return {};

  TopoDS_Shape result;
  try {
    BRepBuilderAPI_Sewing sewing(1e-6 * std::max(d, rb));
    sewing.Add(capA);
    sewing.Add(coneFace);
    sewing.Add(capB);
    sewing.Perform();
    TopoDS_Shell shell;
    for (TopExp_Explorer it(sewing.SewedShape(), TopAbs_SHELL); it.More(); it.Next()) {
      shell = TopoDS::Shell(it.Current());
    }
    if (shell.IsNull()) return {};
    BRepBuilderAPI_MakeSolid maker(shell);
    if (!maker.IsDone()) return {};
    result = maker.Solid();
  } catch (const Standard_Failure&) {
    return {};
  }
  if (!isValid(result)) return {};

  const double h1 = ra * (1 - sinA), h2 = rb * (1 + sinA);
  const double rho1 = ra * cosA, rho2 = rb * cosA;
  const double exact = M_PI * h1 * h1 * (3 * ra - h1) / 3 +
                       M_PI * length / 3 * (rho1 * rho1 + rho1 * rho2 + rho2 * rho2) +
                       M_PI * h2 * h2 * (3 * rb - h2) / 3;
  const double got = volume(result);
  if (std::fabs(got - exact) > SELF_CHECK_RTOL * exact) return {};
  GProp_GProps props;
  BRepGProp::VolumeProperties(result, props);
  if (props.Mass() < 0) result = result.Reversed();
  return result;
}

// --- rungs -----------------------------------------------------------------

TopoDS_Shape hullOfSpheres(const std::vector<TopoDS_Shape>& shapes, std::string& rung)
{
  std::vector<SphereInfo> spheres;
  for (const auto& s : shapes) {
    auto info = asSphere(s);
    if (!info) return {};
    spheres.push_back(*info);
  }
  const double r0 = spheres.front().radius;
  bool equal = true;
  for (const auto& s : spheres) {
    if (std::fabs(s.radius - r0) > RADIUS_REL_TOL * r0) equal = false;
  }
  if (!equal) {
    if (spheres.size() != 2) return {};
    const size_t small = spheres[0].radius <= spheres[1].radius ? 0 : 1;
    const auto& a = spheres[small];
    const auto& b = spheres[1 - small];
    rung = "hull of two spheres";
    if (dist(a.center, b.center) + a.radius <= b.radius * (1 + RADIUS_REL_TOL)) {
      return shapes[1 - small];
    }
    return twoSphereHull(a.center, a.radius, b.center, b.radius);
  }
  std::vector<P3> centers;
  for (const auto& s : spheres) centers.push_back(s.center);
  rung = "hull of equal spheres";
  return hull3dOffset(centers, r0);
}

TopoDS_Shape hullOfCylinders(const std::vector<TopoDS_Shape>& shapes, std::string& rung)
{
  std::vector<CylinderInfo> cylinders;
  for (const auto& s : shapes) {
    auto info = asCylinder(s);
    if (!info) return {};
    cylinders.push_back(*info);
  }
  const double r0 = cylinders.front().radius;
  for (const auto& c : cylinders) {
    if (std::fabs(c.radius - r0) > RADIUS_REL_TOL * r0) return {};
  }
  const auto& first = cylinders.front();
  const double length0 = dist(first.a, first.b);
  if (length0 < 1e-9) return {};
  const P3 dir0{(first.b[0] - first.a[0]) / length0, (first.b[1] - first.a[1]) / length0,
                (first.b[2] - first.a[2]) / length0};

  std::optional<std::pair<double, double>> span;
  std::vector<P3> axisPoints;
  for (auto c : cylinders) {
    const double length = dist(c.a, c.b);
    if (length < 1e-9) return {};
    const P3 dir{(c.b[0] - c.a[0]) / length, (c.b[1] - c.a[1]) / length, (c.b[2] - c.a[2]) / length};
    const double dot = dir[0] * dir0[0] + dir[1] * dir0[1] + dir[2] * dir0[2];
    if (std::fabs(std::fabs(dot) - 1.0) > PARALLEL_TOL) return {};
    if (dot < 0) std::swap(c.a, c.b);
    auto along = [&](const P3& p) {
      return (p[0] - first.a[0]) * dir0[0] + (p[1] - first.a[1]) * dir0[1] +
             (p[2] - first.a[2]) * dir0[2];
    };
    const double ta = along(c.a), tb = along(c.b);
    const std::pair<double, double> thisSpan{std::min(ta, tb), std::max(ta, tb)};
    if (!span) {
      span = thisSpan;
    } else if (std::fabs(thisSpan.first - span->first) > SPAN_TOL * length0 ||
               std::fabs(thisSpan.second - span->second) > SPAN_TOL * length0) {
      return {};
    }
    axisPoints.push_back({(c.a[0] + c.b[0]) / 2, (c.a[1] + c.b[1]) / 2, (c.a[2] + c.b[2]) / 2});
  }

  const gp_Ax3 frame(gp_Pnt(first.a[0], first.a[1], first.a[2]), gp_Dir(dir0[0], dir0[1], dir0[2]));
  std::vector<P2> points;
  for (const auto& p : axisPoints) {
    const gp_Vec rel(p[0] - first.a[0], p[1] - first.a[1], p[2] - first.a[2]);
    points.push_back({rel.Dot(gp_Vec(frame.XDirection())), rel.Dot(gp_Vec(frame.YDirection()))});
  }
  auto face = hull2dOffset(points, r0);
  if (face.IsNull()) return {};
  rung = "hull of parallel cylinders";
  try {
    gp_Trsf toWorld;
    toWorld.SetTransformation(frame);
    toWorld.Invert();
    auto placed = BRepBuilderAPI_Transform(face, toWorld, true).Shape();
    const double height = span->second - span->first;
    BRepPrimAPI_MakePrism prism(placed, gp_Vec(dir0[0] * height, dir0[1] * height, dir0[2] * height));
    if (!prism.IsDone()) return {};
    auto result = prism.Shape();
    return isValid(result) ? result : TopoDS_Shape();
  } catch (const Standard_Failure&) {
    return {};
  }
}

TopoDS_Shape hullOfPolyhedra(const std::vector<TopoDS_Shape>& shapes, std::string& rung)
{
  std::vector<P3> points;
  for (const auto& s : shapes) {
    if (solidCount(s) == 0) return {};
    for (const auto& f : facesOf(s)) {
      if (surfaceType(f) != GeomAbs_Plane) return {};
    }
    for (const auto& p : verticesOf(s)) points.push_back(p);
  }
  if (points.size() < 4) return {};
  rung = "hull of polyhedra";
  auto result = convexHull3d(points);
  return isValid(result) ? result : TopoDS_Shape();
}

TopoDS_Shape hullOfTwoDiscs(const std::vector<TopoDS_Shape>& shapes, std::string& rung)
{
  if (shapes.size() != 2) return {};
  auto a = asDisc(shapes[0]);
  auto b = asDisc(shapes[1]);
  if (!a || !b) return {};
  const size_t small = a->radius <= b->radius ? 0 : 1;
  const CircleInfo& s = small == 0 ? *a : *b;
  const CircleInfo& l = small == 0 ? *b : *a;
  rung = "hull of two discs";
  const double d = dist2(s.center, l.center);
  if (d + s.radius <= l.radius * (1 + RADIUS_REL_TOL)) return shapes[1 - small];
  const double sinA = (l.radius - s.radius) / d;
  const double cosA = std::sqrt(1.0 - sinA * sinA);
  const double ux = (l.center[0] - s.center[0]) / d, uy = (l.center[1] - s.center[1]) / d;
  const double nx = -uy, ny = ux;
  auto graze = [&](const CircleInfo& c, double side) {
    return gp_Pnt(c.center[0] + c.radius * (-sinA * ux + side * cosA * nx),
                  c.center[1] + c.radius * (-sinA * uy + side * cosA * ny), 0);
  };
  const auto p1p = graze(s, 1), p1m = graze(s, -1), p2p = graze(l, 1), p2m = graze(l, -1);
  const gp_Pnt back(s.center[0] - s.radius * ux, s.center[1] - s.radius * uy, 0);
  const gp_Pnt front(l.center[0] + l.radius * ux, l.center[1] + l.radius * uy, 0);
  TopoDS_Shape face;
  try {
    BRepBuilderAPI_MakeWire wire;
    wire.Add(BRepBuilderAPI_MakeEdge(GC_MakeArcOfCircle(p1p, back, p1m).Value()).Edge());
    wire.Add(BRepBuilderAPI_MakeEdge(p1m, p2m).Edge());
    wire.Add(BRepBuilderAPI_MakeEdge(GC_MakeArcOfCircle(p2m, front, p2p).Value()).Edge());
    wire.Add(BRepBuilderAPI_MakeEdge(p2p, p1p).Edge());
    if (!wire.IsDone()) return {};
    face = faceFromWire(wire.Wire());
  } catch (const Standard_Failure&) {
    return {};
  }
  if (!isValid(face)) return {};
  const double alpha = std::asin(sinA);
  const double sin2a = std::sin(2 * alpha);
  const double exact = 0.5 * s.radius * s.radius * (M_PI - 2 * alpha - sin2a) +
                       0.5 * l.radius * l.radius * (M_PI + 2 * alpha + sin2a) +
                       (s.radius + l.radius) * d * cosA * cosA * cosA;
  if (std::fabs(area(face) - exact) > SELF_CHECK_RTOL * exact) return {};
  return face;
}

TopoDS_Shape hullOfPolygons(const std::vector<TopoDS_Shape>& shapes, std::string& rung)
{
  std::vector<P2> points;
  for (const auto& s : shapes) {
    if (solidCount(s) != 0 || facesOf(s).empty()) return {};
    for (const auto& e : edgesOf(s)) {
      if (BRepAdaptor_Curve(TopoDS::Edge(e)).GetType() != GeomAbs_Line) return {};
    }
    for (const auto& v : verticesOf(s)) {
      if (std::fabs(v[2]) > 1e-9) return {};
      points.push_back({v[0], v[1]});
    }
  }
  auto hull = convexHull2d(points);
  if (hull.size() < 3) return {};
  rung = "hull of polygons";
  auto face = polygonFace(hull);
  return isValid(face) ? face : TopoDS_Shape();
}

// --- minkowski -------------------------------------------------------------

std::optional<double> meshBallRadius(const TopoDS_Shape& shape)
{
  const auto points = verticesOf(shape);
  if (points.size() < MIN_MESH_BALL_VERTICES) return std::nullopt;
  P3 c{0, 0, 0};
  for (const auto& p : points) {
    c[0] += p[0] / points.size();
    c[1] += p[1] / points.size();
    c[2] += p[2] / points.size();
  }
  double lo = 1e300, hi = 0, sum = 0;
  for (const auto& p : points) {
    const double d = dist(p, c);
    lo = std::min(lo, d);
    hi = std::max(hi, d);
    sum += d;
  }
  if (hi < 1e-9) return std::nullopt;
  if (hi - lo > MESH_BALL_REL_TOL * hi) return std::nullopt;
  if (dist(c, {0, 0, 0}) > MESH_BALL_REL_TOL * hi) return std::nullopt;
  return sum / points.size();
}

std::optional<double> ballRadius(const TopoDS_Shape& shape, unsigned int dim)
{
  if (dim == 3) {
    if (auto sphere = asSphere(shape)) {
      if (dist(sphere->center, {0, 0, 0}) > ORIGIN_TOL) return std::nullopt;
      return sphere->radius;
    }
    return meshBallRadius(shape);
  }
  if (auto disc = asDisc(shape)) {
    if (dist2(disc->center, {0, 0}) > ORIGIN_TOL) return std::nullopt;
    return disc->radius;
  }
  return meshBallRadius(shape);
}

}  // namespace

std::vector<TopoDS_Shape> components(const std::vector<TopoDS_Shape>& children, unsigned int dim)
{
  std::vector<TopoDS_Shape> out;
  for (const auto& child : children) {
    for (const auto& piece : OcctBoolean::piecesOf(child, dim)) out.push_back(piece);
  }
  return out;
}

TopoDS_Shape hullOfTessellation(const std::vector<TopoDS_Shape>& comps, int segments)
{
  if (comps.empty() || segments < 3) return {};
  try {
    Bnd_Box box;
    for (const auto& c : comps) BRepBndLib::Add(c, box);
    if (box.IsVoid()) return {};
    double x0, y0, z0, x1, y1, z1;
    box.Get(x0, y0, z0, x1, y1, z1);
    const double diag = std::sqrt((x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0) + (z1 - z0) * (z1 - z0));
    // The angle sets the segment count on every curve; the length only
    // keeps a huge flat face from being split further than it needs.
    const double angular = 2.0 * M_PI / segments;
    const double linear = std::max(diag * 1e-3, 1e-6);
    std::vector<P3> points;
    for (const auto& c : comps) {
      // Meshing writes the triangulation into the shape, and a shape that
      // already carries one is not re-meshed, so mesh a copy: the caller's
      // shape stays untouched and the fineness asked for is the fineness
      // delivered.
      const TopoDS_Shape copy = BRepBuilderAPI_Copy(c).Shape();
      BRepMesh_IncrementalMesh mesher(copy, linear, false, angular, false);
      for (TopExp_Explorer it(copy, TopAbs_FACE); it.More(); it.Next()) {
        TopLoc_Location loc;
        const auto tri = BRep_Tool::Triangulation(TopoDS::Face(it.Current()), loc);
        if (tri.IsNull()) continue;
        const gp_Trsf trsf = loc.Transformation();
        for (int i = 1; i <= tri->NbNodes(); ++i) {
          const gp_Pnt p = tri->Node(i).Transformed(trsf);
          points.push_back({p.X(), p.Y(), p.Z()});
        }
      }
    }
    // Faces share their boundary nodes, so most points arrive several
    // times; and the hull's output depends on its input order, so sort for
    // a result that is the same from one run to the next.
    std::sort(points.begin(), points.end());
    points.erase(std::unique(points.begin(), points.end()), points.end());
    if (points.size() < 4) return {};
    return convexHull3d(points);
  } catch (const Standard_Failure&) {
    return {};
  }
}

TopoDS_Shape hull(const std::vector<TopoDS_Shape>& comps, unsigned int dim, std::string& rung)
{
  if (comps.empty()) return {};
  try {
    if (dim == 3) {
      if (auto r = hullOfSpheres(comps, rung); !r.IsNull()) return r;
      if (auto r = hullOfCylinders(comps, rung); !r.IsNull()) return r;
      if (auto r = OcctRevolvedHull::hull(comps); !r.IsNull()) {
        rung = "hull of translated revolution solids";
        return r;
      }
      if (auto r = hullOfPolyhedra(comps, rung); !r.IsNull()) return r;
    } else if (dim == 2) {
      if (auto r = hullOfTwoDiscs(comps, rung); !r.IsNull()) return r;
      if (auto r = hullOfPolygons(comps, rung); !r.IsNull()) return r;
    }
  } catch (const Standard_Failure&) {
    return {};
  }
  rung.clear();
  return {};
}

TopoDS_Shape minkowski(const std::vector<TopoDS_Shape>& children, unsigned int dim, std::string& rung)
{
  if (children.size() != 2) return {};
  try {
    for (const size_t ballIndex : {size_t{1}, size_t{0}}) {
      const auto radius = ballRadius(children[ballIndex], dim);
      if (!radius) continue;
      const auto& target = children[1 - ballIndex];
      TopoDS_Shape result;
      if (dim == 3) {
        std::vector<TopoDS_Shape> parts;
        for (const auto& solid : OcctBoolean::piecesOf(target, 3)) {
          auto part = offset3d(solid, *radius);
          if (part.IsNull()) {
            parts.clear();
            break;
          }
          parts.push_back(part);
        }
        if (!parts.empty()) result = OcctBoolean::fuse(parts, 3);
      } else {
        std::vector<TopoDS_Shape> parts;
        for (const auto& face : OcctBoolean::piecesOf(target, 2)) {
          auto part = offset2d(face, *radius);
          if (part.IsNull()) {
            parts.clear();
            break;
          }
          parts.push_back(part);
        }
        if (!parts.empty()) result = OcctBoolean::fuse(parts, 2);
      }
      if (!result.IsNull()) {
        rung = "minkowski with a ball";
        return result;
      }
    }
  } catch (const Standard_Failure&) {
    return {};
  }
  return {};
}

}  // namespace OcctHull
