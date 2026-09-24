#include "geometry/occt/OcctMesh.h"

#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <functional>
#include <unordered_map>
#include <cstdint>
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
#include <string>
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

// Turn closed shells into solids, nesting any shell that lies inside
// another as that solid's cavity. Null if nothing closes.
TopoDS_Shape solidsFromShells(const std::vector<TopoDS_Shell>& shells, std::string& why)
{
  struct Candidate {
    TopoDS_Shell shell;
    TopoDS_Solid solid;
    double volume;
  };
  std::vector<Candidate> candidates;
  for (auto shell : shells) {
    BRepBuilderAPI_MakeSolid maker(shell);
    if (!maker.IsDone()) {
      why = "a shell could not be made into a solid";
      continue;
    }
    auto solid = maker.Solid();
    if (isInsideOut(solid)) {
      solid = TopoDS::Solid(solid.Reversed());
      shell = TopoDS::Shell(shell.Reversed());
    }
    candidates.push_back({shell, solid, std::fabs(volumeOf(solid))});
  }
  if (candidates.empty()) {
    if (why.empty()) why = "no shell in the sewn result";
    return {};
  }
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

// Builds shells straight from the mesh's own topology.
//
// The mesh already says which faces share which vertices, so the edges can
// be shared by index. Sewing throws that away and rediscovers it by
// comparing edges geometrically, which is what it is for when the faces
// come from unrelated sources, and which dominates a large mesh: two
// corpus models spent their whole export inside
// BRepBuilderAPI_Sewing::FindCandidates.
//
// Returns nothing when the mesh is not a clean closed manifold -- an edge
// used by other than exactly two faces in opposite directions, a repeated
// or degenerate vertex, a wire or face OpenCASCADE will not build. The
// caller then sews, which copes with meshes this cannot.
std::vector<TopoDS_Shell> shellsFromIndexedMesh(const OcctBridge::MeshData& mesh, std::string& why)
{
  if (mesh.faces.empty() || mesh.vertices.empty()) return {};
  if (mesh.vertices.size() > (1ULL << 31)) return {};

  struct Use {
    TopoDS_Edge edge;
    int forward = 0;
    int reverse = 0;
    size_t face = 0;
  };
  const auto key = [](size_t a, size_t b) {
    return a < b ? (static_cast<uint64_t>(a) << 32) | b : (static_cast<uint64_t>(b) << 32) | a;
  };

  std::vector<TopoDS_Vertex> vertices(mesh.vertices.size());
  std::unordered_map<uint64_t, Use> edges;
  edges.reserve(mesh.faces.size() * 2);

  try {
    for (size_t f = 0; f < mesh.faces.size(); ++f) {
      const auto& face = mesh.faces[f];
      if (face.size() < 3) return {};
      for (size_t i = 0; i < face.size(); ++i) {
        const size_t a = face[i];
        const size_t b = face[(i + 1) % face.size()];
        if (a == b || a >= mesh.vertices.size() || b >= mesh.vertices.size()) return {};
        auto& use = edges[key(a, b)];
        if (use.edge.IsNull()) {
          for (const size_t v : {std::min(a, b), std::max(a, b)}) {
            if (vertices[v].IsNull()) {
              const auto& p = mesh.vertices[v];
              vertices[v] = BRepBuilderAPI_MakeVertex(gp_Pnt(p[0], p[1], p[2])).Vertex();
            }
          }
          BRepBuilderAPI_MakeEdge maker(vertices[std::min(a, b)], vertices[std::max(a, b)]);
          if (!maker.IsDone()) return {};
          use.edge = maker.Edge();
          use.face = f;
        }
        ++(a < b ? use.forward : use.reverse);
      }
    }

    // Every edge of a closed, consistently wound surface is walked once in
    // each direction. Anything else is a hole, a seam or a flipped face,
    // and not something to guess at.
    for (const auto& [ignored, use] : edges) {
      if (use.forward != 1 || use.reverse != 1) {
        why = "the mesh is not a closed manifold";
        return {};
      }
    }

    // Faces that share an edge belong to the same shell.
    std::vector<size_t> parent(mesh.faces.size());
    for (size_t i = 0; i < parent.size(); ++i) parent[i] = i;
    const std::function<size_t(size_t)> root = [&](size_t i) {
      while (parent[i] != i) i = parent[i] = parent[parent[i]];
      return i;
    };

    std::vector<TopoDS_Face> built(mesh.faces.size());
    for (size_t f = 0; f < mesh.faces.size(); ++f) {
      const auto& face = mesh.faces[f];
      BRepBuilderAPI_MakeWire wire;
      for (size_t i = 0; i < face.size(); ++i) {
        const size_t a = face[i];
        const size_t b = face[(i + 1) % face.size()];
        auto& use = edges[key(a, b)];
        wire.Add(TopoDS::Edge(a < b ? use.edge : use.edge.Reversed()));
        if (!wire.IsDone()) return {};
        const size_t other = root(use.face);
        if (other != root(f)) parent[other] = root(f);
      }
      BRepBuilderAPI_MakeFace maker(wire.Wire(), true);
      if (!maker.IsDone()) return {};
      built[f] = maker.Face();
    }

    std::unordered_map<size_t, TopoDS_Shell> shells;
    BRep_Builder builder;
    for (size_t f = 0; f < built.size(); ++f) {
      auto& shell = shells[root(f)];
      if (shell.IsNull()) builder.MakeShell(shell);
      builder.Add(shell, built[f]);
    }
    std::vector<TopoDS_Shell> out;
    out.reserve(shells.size());
    for (auto& [ignored, shell] : shells) {
      shell.Closed(BRep_Tool::IsClosed(shell));
      out.push_back(shell);
    }
    return out;
  } catch (const Standard_Failure&) {
    return {};
  }
}

// Sew faces into shells, then solids. Sewing matches edges by comparing
// them geometrically, which is what it is for when the faces come from
// unrelated sources, and what makes it slow on a large mesh.
TopoDS_Shape solidsFromFaceList(const std::vector<TopoDS_Shape>& faces, double tolerance, int& freeEdges,
                                std::string& why)
{
  freeEdges = 0;
  if (faces.empty()) {
    why = "no faces";
    return {};
  }
  BRepBuilderAPI_Sewing sewing(tolerance);
  for (const auto& f : faces) sewing.Add(f);
  // Deliberately not given the export's time budget: sewing a mesh is the
  // fallback the budget falls back *to*, so cutting it short would leave
  // the node with nothing at all.
  sewing.Perform();
  const auto sewn = sewing.SewedShape();
  if (sewn.IsNull()) {
    why = "sewing produced nothing";
    return {};
  }
  freeEdges = countFreeEdges(sewn);
  if (freeEdges > 0) {
    why = std::to_string(freeEdges) + " free edge(s) after sewing " + std::to_string(faces.size()) +
          " face(s)";
    return {};
  }
  std::vector<TopoDS_Shell> shells;
  for (TopExp_Explorer it(sewn, TopAbs_SHELL); it.More(); it.Next()) {
    shells.push_back(TopoDS::Shell(it.Current()));
  }
  return solidsFromShells(shells, why);
}

}  // namespace

namespace {

// Remove zero-area triangles (three collinear vertices) by dropping the
// triangle and inserting its middle vertex into the neighbor that shares
// the long edge, so every edge still has exactly two faces. Manifold's
// hull emits such triangles, and a dropped one leaves a T-junction that
// sewing does not always close.
OcctBridge::MeshData repairCollinear(const OcctBridge::MeshData& input)
{
  OcctBridge::MeshData mesh = input;
  const auto& V = mesh.vertices;
  auto middleOf = [&](const std::vector<size_t>& tri, size_t& p, size_t& m, size_t& q) {
    // The middle vertex projects between the other two along their line.
    for (int k = 0; k < 3; ++k) {
      const size_t a = tri[k], b = tri[(k + 1) % 3], c = tri[(k + 2) % 3];
      const Vector3d ac = V[c] - V[a];
      const double len2 = ac.squaredNorm();
      if (len2 < 1e-24) continue;
      const double t = (V[b] - V[a]).dot(ac) / len2;
      if (t > 0 && t < 1) {
        p = a;
        m = b;
        q = c;
        return true;
      }
    }
    return false;
  };
  for (int pass = 0; pass < 8; ++pass) {
    bool changed = false;
    for (size_t i = 0; i < mesh.faces.size(); ++i) {
      const auto& tri = mesh.faces[i];
      if (tri.size() != 3) continue;
      const Vector3d e1 = V[tri[1]] - V[tri[0]], e2 = V[tri[2]] - V[tri[0]];
      const double longest = std::max({e1.norm(), e2.norm(), (V[tri[2]] - V[tri[1]]).norm()});
      if (longest < 1e-12) continue;
      if (e1.cross(e2).norm() > 1e-9 * longest * longest) continue;
      size_t p, m, q;
      if (!middleOf(tri, p, m, q)) continue;
      // The neighbor across edge (p, q) gets m inserted between them.
      for (size_t j = 0; j < mesh.faces.size() && !changed; ++j) {
        if (j == i) continue;
        auto& ring = mesh.faces[j];
        for (size_t k = 0; k < ring.size(); ++k) {
          const size_t a = ring[k], b = ring[(k + 1) % ring.size()];
          if ((a == p && b == q) || (a == q && b == p)) {
            ring.insert(ring.begin() + static_cast<long>(k) + 1, m);
            changed = true;
            break;
          }
        }
      }
      if (changed) {
        mesh.faces.erase(mesh.faces.begin() + static_cast<long>(i));
        break;
      }
    }
    if (!changed) break;
    --pass;  // keep going while triangles are still being repaired
    if (mesh.faces.size() < 4) break;
  }
  return mesh;
}

}  // namespace

TopoDS_Shape solidsFromMesh(const OcctBridge::MeshData& input, std::string& why)
{
  const auto mesh = repairCollinear(input);
  std::vector<TopoDS_Shape> faces;
  size_t dropped = 0;
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
    else ++dropped;
  }
  // The indexed build is exact and cheap when the mesh is clean; sewing is
  // the general case and is what handles everything else.
  std::string indexedWhy;
  const auto shells = shellsFromIndexedMesh(mesh, indexedWhy);
  if (!shells.empty()) {
    std::string shellWhy;
    auto direct = solidsFromShells(shells, shellWhy);
    if (!direct.IsNull()) return direct;
  }

  int freeEdges = 0;
  auto result = solidsFromFaceList(faces, 1e-6 * scale, freeEdges, why);
  if (result.IsNull() && dropped > 0) {
    why += " (" + std::to_string(dropped) + " of " + std::to_string(mesh.faces.size()) +
           " polygons could not be made into faces)";
  }
  return result;
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
  std::string why;
  return solidsFromFaceList(built, 1e-7 * scale, freeEdges, why);
}

}  // namespace OcctMesh
