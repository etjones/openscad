#include "geometry/occt/OcctMesh.h"

#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakePolygon.hxx>
#include <BRepBuilderAPI_MakeSolid.hxx>
#include <BRepBuilderAPI_Sewing.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepClass3d_SolidClassifier.hxx>
#include <BRepGProp.hxx>
#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <GProp_GProps.hxx>
#include <ShapeAnalysis_FreeBounds.hxx>
#include <Standard_Failure.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopTools_IndexedDataMapOfShapeListOfShape.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Compound.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Shell.hxx>
#include <TopoDS_Solid.hxx>
#include <TopoDS_Vertex.hxx>
#include <TopoDS_Wire.hxx>
#include <algorithm>
#include <cmath>
#include <memory>
#include <vector>

#include "geometry/ClipperUtils.h"
#include "geometry/Polygon2d.h"
#include "geometry/occt/OcctBoolean.h"

namespace OcctMesh {

namespace {

double signedArea(const VectorOfVector2d& ring)
{
  double area = 0;
  for (size_t i = 0, n = ring.size(); i < n; ++i) {
    const auto& a = ring[i];
    const auto& b = ring[(i + 1) % n];
    area += a[0] * b[1] - b[0] * a[1];
  }
  return area / 2;
}

bool pointInRing(const Vector2d& p, const VectorOfVector2d& ring)
{
  bool in = false;
  for (size_t i = 0, j = ring.size() - 1; i < ring.size(); j = i++) {
    const auto& a = ring[i];
    const auto& b = ring[j];
    if ((a[1] > p[1]) != (b[1] > p[1]) && p[0] < (b[0] - a[0]) * (p[1] - a[1]) / (b[1] - a[1]) + a[0]) {
      in = !in;
    }
  }
  return in;
}

// Drop consecutive duplicates (and a closing repeat of the first point).
VectorOfVector2d dedupe(const VectorOfVector2d& ring)
{
  VectorOfVector2d out;
  for (const auto& v : ring) {
    if (out.empty() || (out.back() - v).norm() > 1e-12) out.push_back(v);
  }
  if (out.size() > 1 && (out.front() - out.back()).norm() <= 1e-12) out.pop_back();
  return out;
}

TopoDS_Wire wireFromRing2d(const VectorOfVector2d& ring, bool counterClockwise)
{
  auto pts = dedupe(ring);
  if (pts.size() < 3) return {};
  if ((signedArea(pts) > 0) != counterClockwise) std::reverse(pts.begin(), pts.end());
  BRepBuilderAPI_MakePolygon polygon;
  for (const auto& v : pts) polygon.Add(gp_Pnt(v[0], v[1], 0));
  polygon.Close();
  if (!polygon.IsDone()) return {};
  return polygon.Wire();
}

struct Outer {
  const Outline2d *outline;
  double area;
  std::vector<const Outline2d *> holes;
};

}  // namespace

TopoDS_Shape facesFromPolygon2d(const Polygon2d& polygon)
{
  std::unique_ptr<Polygon2d> sanitized;
  const Polygon2d *poly = &polygon;
  if (!polygon.isSanitized()) {
    sanitized = ClipperUtils::sanitize(polygon);
    poly = sanitized.get();
  }
  if (!poly || poly->isEmpty()) return {};

  std::vector<Outer> outers;
  std::vector<const Outline2d *> holes;
  for (const auto& outline : poly->outlines()) {
    if (outline.vertices.size() < 3) continue;
    if (outline.positive) {
      outers.push_back({&outline, std::fabs(signedArea(outline.vertices)), {}});
    } else {
      holes.push_back(&outline);
    }
  }
  // A hole belongs to the smallest positive outline that contains it.
  for (const auto *hole : holes) {
    Outer *best = nullptr;
    for (auto& outer : outers) {
      if (pointInRing(hole->vertices.front(), outer.outline->vertices) &&
          (!best || outer.area < best->area)) {
        best = &outer;
      }
    }
    if (best) best->holes.push_back(hole);
  }

  std::vector<TopoDS_Shape> faces;
  for (const auto& outer : outers) {
    auto wire = wireFromRing2d(outer.outline->vertices, true);
    if (wire.IsNull()) continue;
    BRepBuilderAPI_MakeFace maker(wire, true);
    if (!maker.IsDone()) continue;
    for (const auto *hole : outer.holes) {
      auto holeWire = wireFromRing2d(hole->vertices, false);
      if (!holeWire.IsNull()) maker.Add(holeWire);
    }
    if (maker.IsDone()) faces.push_back(maker.Face());
  }
  if (faces.empty()) return {};
  return faces.size() == 1 ? faces.front() : OcctBoolean::makeCompound(faces);
}

TopoDS_Shape faceFromRing(const std::vector<Vector3d>& ring)
{
  std::vector<Vector3d> pts;
  for (const auto& v : ring) {
    if (pts.empty() || (pts.back() - v).norm() > 1e-12) pts.push_back(v);
  }
  if (pts.size() > 1 && (pts.front() - pts.back()).norm() <= 1e-12) pts.pop_back();
  if (pts.size() < 3) return {};
  try {
    BRepBuilderAPI_MakePolygon polygon;
    for (const auto& v : pts) polygon.Add(gp_Pnt(v[0], v[1], v[2]));
    polygon.Close();
    if (!polygon.IsDone()) return {};
    BRepBuilderAPI_MakeFace maker(polygon.Wire(), true);
    if (!maker.IsDone()) return {};
    return maker.Face();
  } catch (const Standard_Failure&) {
    return {};
  }
}

namespace {

bool isInsideOut(const TopoDS_Solid& solid)
{
  BRepClass3d_SolidClassifier classifier(solid);
  classifier.PerformInfinitePoint(1e-7);
  return classifier.State() == TopAbs_IN;
}

double volumeOf(const TopoDS_Shape& shape)
{
  GProp_GProps props;
  BRepGProp::VolumeProperties(shape, props);
  return props.Mass();
}

gp_Pnt anyVertex(const TopoDS_Shape& shape)
{
  TopExp_Explorer it(shape, TopAbs_VERTEX);
  return BRep_Tool::Pnt(TopoDS::Vertex(it.Current()));
}

int countFreeEdges(const TopoDS_Shape& shape)
{
  TopTools_IndexedDataMapOfShapeListOfShape edgeFaces;
  TopExp::MapShapesAndAncestors(shape, TopAbs_EDGE, TopAbs_FACE, edgeFaces);
  int free = 0;
  for (int i = 1; i <= edgeFaces.Extent(); ++i) {
    if (edgeFaces.FindFromIndex(i).Extent() < 2) ++free;
  }
  return free;
}

// Sew faces into shells, then solids with cavities. Null if nothing closes.
TopoDS_Shape solidsFromFaceList(const std::vector<TopoDS_Shape>& faces, double tolerance, int& freeEdges)
{
  freeEdges = 0;
  if (faces.empty()) return {};
  BRepBuilderAPI_Sewing sewing(tolerance);
  for (const auto& f : faces) sewing.Add(f);
  sewing.Perform();
  const auto sewn = sewing.SewedShape();
  if (sewn.IsNull()) return {};
  freeEdges = countFreeEdges(sewn);
  if (freeEdges > 0) return {};

  struct Candidate {
    TopoDS_Shell shell;
    TopoDS_Solid solid;
    double volume;
  };
  std::vector<Candidate> candidates;
  for (TopExp_Explorer it(sewn, TopAbs_SHELL); it.More(); it.Next()) {
    auto shell = TopoDS::Shell(it.Current());
    BRepBuilderAPI_MakeSolid maker(shell);
    if (!maker.IsDone()) continue;
    auto solid = maker.Solid();
    if (isInsideOut(solid)) {
      solid = TopoDS::Solid(solid.Reversed());
      shell = TopoDS::Shell(shell.Reversed());
    }
    candidates.push_back({shell, solid, std::fabs(volumeOf(solid))});
  }
  if (candidates.empty()) return {};
  std::sort(candidates.begin(), candidates.end(),
            [](const Candidate& a, const Candidate& b) { return a.volume > b.volume; });

  // A shell lying inside a larger one is that solid's cavity.
  std::vector<std::vector<size_t>> cavities(candidates.size());
  std::vector<bool> isCavity(candidates.size(), false);
  for (size_t i = 1; i < candidates.size(); ++i) {
    const auto point = anyVertex(candidates[i].shell);
    for (size_t j = 0; j < i; ++j) {
      if (isCavity[j]) continue;
      BRepClass3d_SolidClassifier classifier(candidates[j].solid, point, 1e-7);
      if (classifier.State() == TopAbs_IN) {
        cavities[j].push_back(i);
        isCavity[i] = true;
        break;
      }
    }
  }

  std::vector<TopoDS_Shape> solids;
  for (size_t j = 0; j < candidates.size(); ++j) {
    if (isCavity[j]) continue;
    if (cavities[j].empty()) {
      solids.push_back(candidates[j].solid);
      continue;
    }
    BRepBuilderAPI_MakeSolid maker(candidates[j].shell);
    for (const size_t i : cavities[j]) {
      maker.Add(TopoDS::Shell(candidates[i].shell.Reversed()));
    }
    solids.push_back(maker.IsDone() ? TopoDS_Shape(maker.Solid()) : candidates[j].solid);
  }
  if (solids.empty()) return {};
  auto result = solids.size() == 1 ? solids.front() : OcctBoolean::makeCompound(solids);
  return OcctBoolean::unify(result, 3);
}

}  // namespace

TopoDS_Shape solidsFromMesh(const OcctBridge::MeshData& mesh)
{
  std::vector<TopoDS_Shape> faces;
  double scale = 1.0;
  for (const auto& v : mesh.vertices) scale = std::max(scale, v.norm());
  for (const auto& polygon : mesh.faces) {
    std::vector<Vector3d> ring;
    ring.reserve(polygon.size());
    for (const auto index : polygon) {
      if (index < mesh.vertices.size()) ring.push_back(mesh.vertices[index]);
    }
    auto face = faceFromRing(ring);
    if (!face.IsNull()) faces.push_back(face);
  }
  int freeEdges = 0;
  return solidsFromFaceList(faces, 1e-6 * scale, freeEdges);
}

TopoDS_Shape solidsFromFaces(const std::vector<Vector3d>& points,
                             const std::vector<std::vector<size_t>>& faces, bool& nonPlanar,
                             int& freeEdges)
{
  nonPlanar = false;
  freeEdges = 0;
  std::vector<TopoDS_Shape> built;
  double scale = 1.0;
  for (const auto& p : points) scale = std::max(scale, p.norm());
  for (const auto& face : faces) {
    std::vector<Vector3d> ring;
    // OpenSCAD orders polyhedron faces clockwise from outside; OCCT wants
    // counterclockwise for an outward normal.
    for (auto it = face.rbegin(); it != face.rend(); ++it) {
      if (*it < points.size()) ring.push_back(points[*it]);
    }
    if (ring.size() < 3) continue;
    auto f = faceFromRing(ring);
    if (f.IsNull()) {
      nonPlanar = true;
      return {};
    }
    built.push_back(f);
  }
  return solidsFromFaceList(built, 1e-7 * scale, freeEdges);
}

}  // namespace OcctMesh
