#include "geometry/occt/OcctBuilder.h"

#include <BRepBuilderAPI_GTransform.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <BRepBndLib.hxx>
#include <BRepGProp.hxx>
#include <BRepLib.hxx>
#include <BRepOffsetAPI_ThruSections.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeCone.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepPrimAPI_MakePrism.hxx>
#include <BRepPrimAPI_MakeRevol.hxx>
#include <BRepPrimAPI_MakeSphere.hxx>
#include <BRepTools_WireExplorer.hxx>
#include <Bnd_Box.hxx>
#include <GProp_GProps.hxx>
#include <ShapeAnalysis.hxx>
#include <Standard_Failure.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Solid.hxx>
#include <TopoDS_Wire.hxx>
#include <gp_Ax1.hxx>
#include <gp_Ax2.hxx>
#include <gp_Circ.hxx>
#include <gp_GTrsf.hxx>
#include <gp_Trsf.hxx>
#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "core/CgalAdvNode.h"
#include "core/ColorNode.h"
#include "core/CsgOpNode.h"
#include "core/CurveDiscretizer.h"
#include "core/ImportNode.h"
#include "core/LinearExtrudeNode.h"
#include "core/ModuleInstantiation.h"
#include "core/RotateExtrudeNode.h"
#include "core/SurfaceNode.h"
#include "core/TextNode.h"
#include "core/TransformNode.h"
#include "core/enums.h"
#include "core/primitives.h"
#include "geometry/ClipperUtils.h"
#include "geometry/Geometry.h"
#include "geometry/GeometryEvaluator.h"
#include "geometry/Polygon2d.h"
#include "geometry/occt/OcctBoolean.h"
#include "geometry/occt/OcctMesh.h"
#include "utils/degree_trig.h"
#include "geometry/occt/OcctBridge.h"
#include "geometry/occt/OcctHull.h"

namespace {

constexpr double SINGULAR_EPS = 1e-12;
// A matrix within this of a uniform scale times a rotation is treated as
// exactly that. OpenSCAD's .csg text carries six significant figures, so a
// rotation read back from it is orthonormal only to ~5e-7; treating that
// as shear would send an exact cylinder through a B-spline approximation.
constexpr double RIGID_TOL = 1e-6;
constexpr double EXTENT_EPS = 1e-9;

OcctGeometry single(const TopoDS_Shape& shape, unsigned int dim, const Color4f& color)
{
  OcctGeometry geometry;
  if (shape.IsNull()) return geometry;
  geometry.dim = dim;
  geometry.bodies.push_back({shape, color});
  return geometry;
}

std::vector<TopoDS_Shape> shapesOf(const OcctGeometry& geometry)
{
  std::vector<TopoDS_Shape> shapes;
  for (const auto& body : geometry.bodies) shapes.push_back(body.shape);
  return shapes;
}

bool sameColor(const Color4f& a, const Color4f& b)
{
  if (!a.isValid() && !b.isValid()) return true;
  if (a.isValid() != b.isValid()) return false;
  return a == b;
}

double signedVolume(const TopoDS_Shape& shape)
{
  GProp_GProps props;
  BRepGProp::VolumeProperties(shape, props);
  return props.Mass();
}

// Every solid oriented outward. A reflected or reversed solid encloses
// negative volume; a REVERSED-flagged one is dropped by OCCT's N-ary fuse
// wherever it overlaps another, so orient the faces properly instead.
TopoDS_Shape orientedOutward(const TopoDS_Shape& shape)
{
  if (shape.IsNull()) return shape;
  std::vector<TopoDS_Shape> solids;
  bool changed = false;
  for (TopExp_Explorer it(shape, TopAbs_SOLID); it.More(); it.Next()) {
    auto solid = TopoDS::Solid(it.Current());
    if (signedVolume(solid) < 0) {
      TopoDS_Solid forward = TopoDS::Solid(solid.Oriented(TopAbs_FORWARD));
      if (!BRepLib::OrientClosedSolid(forward)) forward = TopoDS::Solid(solid.Reversed());
      solid = forward;
      changed = true;
    }
    solids.push_back(solid);
  }
  if (!changed) return shape;
  return solids.size() == 1 ? solids.front() : OcctBoolean::makeCompound(solids);
}

TopoDS_Shape translated(const TopoDS_Shape& shape, double x, double y, double z)
{
  gp_Trsf trsf;
  trsf.SetTranslation(gp_Vec(x, y, z));
  return BRepBuilderAPI_Transform(shape, trsf, false).Shape();
}

std::vector<Vector2d> ngon(double radius, int count)
{
  std::vector<Vector2d> points;
  points.reserve(count);
  for (int i = 0; i < count; ++i) {
    const double phi = (360.0 * i) / count;
    points.emplace_back(radius * cos_degrees(phi), radius * sin_degrees(phi));
  }
  return points;
}

struct Decomposed {
  double scale;
  double rotation[3][3];
};

// Split a 3x3 into a uniform scale and an orthonormal matrix (a rotation,
// or a reflection when the determinant is negative). Nullopt when the
// matrix carries non-uniform scale or shear.
std::optional<Decomposed> decompose(const Transform3d& m)
{
  double cols[3][3];
  double norms[3];
  for (int c = 0; c < 3; ++c) {
    double n = 0;
    for (int r = 0; r < 3; ++r) {
      cols[c][r] = m(r, c);
      n += cols[c][r] * cols[c][r];
    }
    norms[c] = std::sqrt(n);
    if (norms[c] < SINGULAR_EPS) return std::nullopt;
  }
  const double scale = (norms[0] + norms[1] + norms[2]) / 3;
  for (const double n : norms) {
    if (std::fabs(n - scale) > RIGID_TOL * std::max(scale, 1.0)) return std::nullopt;
  }
  // Gram-Schmidt on the unit columns.
  double u[3][3];
  for (int c = 0; c < 3; ++c) {
    for (int r = 0; r < 3; ++r) u[c][r] = cols[c][r] / norms[c];
    for (int p = 0; p < c; ++p) {
      double dot = 0;
      for (int r = 0; r < 3; ++r) dot += u[c][r] * u[p][r];
      if (std::fabs(dot) > 1e-6) return std::nullopt;  // genuinely sheared
      for (int r = 0; r < 3; ++r) u[c][r] -= dot * u[p][r];
    }
    double n = 0;
    for (int r = 0; r < 3; ++r) n += u[c][r] * u[c][r];
    n = std::sqrt(n);
    for (int r = 0; r < 3; ++r) u[c][r] /= n;
  }
  Decomposed d{scale, {}};
  for (int r = 0; r < 3; ++r) {
    for (int c = 0; c < 3; ++c) d.rotation[r][c] = u[c][r];
  }
  return d;
}

}  // namespace

OcctBuilder::OcctBuilder(const Tree& tree, int facetThreshold)
  : tree_(tree), facetThreshold_(facetThreshold)
{
}

OcctBuilder::~OcctBuilder() = default;

OcctGeometry OcctBuilder::build(const AbstractNode& node)
{
  return buildNode(node, Color4f());
}

void OcctBuilder::warn(const AbstractNode& node, const std::string& message) const
{
  OcctBridge::warn(node, tree_, message);
}

std::optional<int> OcctBuilder::facetCount(const CurveDiscretizer& discretizer, double r) const
{
  if (facetThreshold_ <= 0) return std::nullopt;
  const auto count = discretizer.explicitSegmentCount(r);
  if (!count || *count >= facetThreshold_) return std::nullopt;
  return count;
}

std::vector<OcctGeometry> OcctBuilder::buildChildren(const AbstractNode& node, const Color4f& inherited)
{
  std::vector<OcctGeometry> children;
  for (const auto& child : node.getChildren()) {
    if (child->modinst->isBackground()) continue;
    children.push_back(buildNode(*child, inherited));
  }
  return children;
}

OcctGeometry OcctBuilder::buildNode(const AbstractNode& node, const Color4f& inherited)
{
  // The same subtree turns up many times in a model that scatters a shape:
  // one corpus pot places 810 copies of a single vent, each of which would
  // otherwise be built, rendered and sewn from scratch. The key is the
  // node's id string, which the tree memoises, plus the colour inherited
  // from above, since that is carried on the bodies rather than the shape.
  std::string key = OcctBridge::idString(tree_, node);
  if (inherited.isValid()) {
    float r = 0, g = 0, b = 0, a = 0;
    inherited.getRgba(r, g, b, a);
    key += "|" + std::to_string(r) + "," + std::to_string(g) + "," + std::to_string(b) + "," +
           std::to_string(a);
  }
  if (const auto it = cache_.find(key); it != cache_.end()) return it->second;
  auto built = buildNodeUncached(node, inherited);
  cache_.emplace(std::move(key), built);
  return built;
}

OcctGeometry OcctBuilder::buildNodeUncached(const AbstractNode& node, const Color4f& inherited)
{
  if (const auto *n = dynamic_cast<const CubeNode *>(&node)) return cube(*n, inherited);
  if (const auto *n = dynamic_cast<const SphereNode *>(&node)) return sphere(*n, inherited);
  if (const auto *n = dynamic_cast<const CylinderNode *>(&node)) return cylinder(*n, inherited);
  if (const auto *n = dynamic_cast<const PolyhedronNode *>(&node)) return polyhedron(*n, inherited);
  if (const auto *n = dynamic_cast<const SquareNode *>(&node)) return square(*n, inherited);
  if (const auto *n = dynamic_cast<const CircleNode *>(&node)) return circle(*n, inherited);
  if (const auto *n = dynamic_cast<const PolygonNode *>(&node)) return polygon(*n, inherited);
  if (const auto *n = dynamic_cast<const TextNode *>(&node)) return text(*n, inherited);
  if (const auto *n = dynamic_cast<const TransformNode *>(&node)) return transform(*n, inherited);
  if (const auto *n = dynamic_cast<const ColorNode *>(&node)) {
    // The preview keeps the outermost color: an inner color() never
    // overrides one already set above it.
    const Color4f color = inherited.isValid() ? inherited : n->color;
    return unionOf(node, buildChildren(node, color));
  }
  if (const auto *n = dynamic_cast<const CsgOpNode *>(&node)) {
    switch (n->type) {
    case OpenSCADOperator::DIFFERENCE:   return differenceOf(node, buildChildren(node, inherited));
    case OpenSCADOperator::INTERSECTION: return intersectionOf(node, buildChildren(node, inherited));
    default:                             return unionOf(node, buildChildren(node, inherited));
    }
  }
  if (const auto *n = dynamic_cast<const LinearExtrudeNode *>(&node))
    return linearExtrude(*n, inherited);
  if (const auto *n = dynamic_cast<const RotateExtrudeNode *>(&node))
    return rotateExtrude(*n, inherited);
  if (const auto *n = dynamic_cast<const CgalAdvNode *>(&node)) {
    if (n->type == CgalAdvType::HULL || n->type == CgalAdvType::MINKOWSKI) {
      return hullOrMinkowski(*n, inherited);
    }
  }
  // Every remaining leaf (import, surface) and polygon-producing node
  // (projection, offset, roof) is rendered by OpenSCAD's own evaluator.
  if (dynamic_cast<const CgalAdvNode *>(&node) || dynamic_cast<const AbstractPolyNode *>(&node)) {
    return meshFallback(node, inherited, node.name() + "() has no B-rep equivalent");
  }
  // Groups, lists, render(), the root, and anything else: implicit union.
  return unionOf(node, buildChildren(node, inherited));
}

// --- combining -------------------------------------------------------------

namespace {

// Keep the non-empty children of one dimension (the first seen), warning
// about the others, as OpenSCAD does.
std::vector<OcctGeometry> sameDimension(const std::vector<OcctGeometry>& children, unsigned int& dim,
                                        bool& mixed)
{
  std::vector<OcctGeometry> kept;
  dim = 0;
  mixed = false;
  for (const auto& child : children) {
    if (child.isEmpty()) continue;
    if (dim == 0) dim = child.dim;
    if (child.dim != dim) {
      mixed = true;
      continue;
    }
    kept.push_back(child);
  }
  return kept;
}

}  // namespace

TopoDS_Shape OcctBuilder::checkedFuse(const std::vector<TopoDS_Shape>& operands, unsigned int dim)
{
  // A union that cannot be verified hands back its operands unjoined,
  // which is the right material in the right places -- only unmerged. It
  // is not worth trading the whole node's exact surfaces for a mesh, so
  // unlike a cut it does not mark the node unverified.
  return OcctBoolean::fuse(operands, dim);
}

TopoDS_Shape OcctBuilder::checkedCut(const std::vector<TopoDS_Shape>& args,
                                     const std::vector<TopoDS_Shape>& tools, unsigned int dim)
{
  bool ok = true;
  auto shape = OcctBoolean::cut(args, tools, dim, &ok);
  if (!ok) unverified_ = true;
  return shape;
}

TopoDS_Shape OcctBuilder::checkedCommon(const std::vector<TopoDS_Shape>& args,
                                        const std::vector<TopoDS_Shape>& tools, unsigned int dim)
{
  bool ok = true;
  auto shape = OcctBoolean::common(args, tools, dim, &ok);
  if (!ok) unverified_ = true;
  return shape;
}

namespace {

// A hull no closed form covers is normally OpenSCAD's own render of it,
// so the region matches the STL exactly. Past this many triangles that
// render is instead replaced by the hull of the children tessellated at
// HULL_FALLBACK_SEGMENTS: OpenCASCADE's booleans are badly superlinear in
// facet count, and a hull rendered at $fn=180 carried 17k to 33k triangles
// that 2k describe to 0.1% of its volume, and never finished. Below the
// ceiling the model's fineness is honoured, since a coarser model rendered
// coarser is also cheaper; above it the exporter's fineness wins.
constexpr size_t HULL_MESH_BUDGET = 10000;
constexpr int HULL_FALLBACK_SEGMENTS = 60;

// Restores the builder's "unverified" flag on scope exit, so a node that
// falls back to a mesh does not taint its ancestors.
struct VerificationScope {
  bool& flag;
  bool saved;
  explicit VerificationScope(bool& f) : flag(f), saved(f) { flag = false; }
  ~VerificationScope() { flag = saved; }
  [[nodiscard]] bool failed() const { return flag; }
};

const char *UNVERIFIED = "a boolean below it could not be verified in B-rep (see the warning above)";

}  // namespace

OcctGeometry OcctBuilder::unionOf(const AbstractNode& node, std::vector<OcctGeometry> children)
{
  const VerificationScope scope(unverified_);
  auto result = unionOfChecked(node, std::move(children));
  if (!scope.failed()) return result;
  const Color4f color = result.isEmpty() ? Color4f() : result.bodies.front().color;
  return meshFallback(node, color, UNVERIFIED);
}

// The union itself is not what marks a node unverified (see checkedFuse);
// this catches a cut or intersection that failed inside a child of it.
OcctGeometry OcctBuilder::unionOfChecked(const AbstractNode& node, std::vector<OcctGeometry> children)
{
  unsigned int dim = 0;
  bool mixed = false;
  auto kept = sameDimension(children, dim, mixed);
  if (mixed) warn(node, "Mixing 2D and 3D objects is not supported");
  if (kept.empty()) return {};
  if (kept.size() == 1) return kept.front();

  std::vector<OcctBody> bodies;
  for (auto& child : kept) {
    for (auto& body : child.bodies) bodies.push_back(std::move(body));
  }
  bool oneColor = true;
  for (size_t i = 1; i < bodies.size(); ++i) {
    if (!sameColor(bodies[i].color, bodies[0].color)) {
      oneColor = false;
      break;
    }
  }
  if (oneColor) {
    std::vector<TopoDS_Shape> shapes;
    for (const auto& body : bodies) shapes.push_back(body.shape);
    return single(checkedFuse(shapes, dim), dim, bodies[0].color);
  }
  OcctGeometry flat;
  flat.dim = dim;
  flat.bodies = std::move(bodies);
  return partitionedUnion({flat}, dim);
}

// Union of bodies carrying different colors: fuse runs of one color, and
// let higher-priority bodies claim contested material. An assigned color
// wins over uncolored material; between two assigned colors the later
// operand wins. Bodies that share no material stay as they are.
OcctGeometry OcctBuilder::partitionedUnion(std::vector<OcctGeometry> children, unsigned int dim)
{
  std::vector<OcctBody> bodies;
  for (auto& child : children) {
    for (auto& body : child.bodies) bodies.push_back(std::move(body));
  }
  std::vector<TopoDS_Shape> all;
  for (const auto& body : bodies) all.push_back(body.shape);
  const auto fused = checkedFuse(all, dim);
  const double fusedExtent = OcctBoolean::extent(fused, dim);
  double naive = 0;
  for (const auto& s : all) naive += OcctBoolean::extent(s, dim);

  OcctGeometry result;
  result.dim = dim;
  // Mass properties of curved faces carry integration noise around 1e-8
  // relative; a stricter tolerance sent a set of merely touching bodies
  // into the partition below, where a cut against a coincident surface
  // can come back empty.
  if (std::fabs(fusedExtent - naive) <= 1e-6 * std::max(naive, 1.0)) {
    // Disjoint, or merely touching: nothing to partition.
    result.bodies = std::move(bodies);
    return result;
  }

  // Ascending priority: uncolored in source order, then colored in source order.
  std::vector<size_t> order;
  for (size_t i = 0; i < bodies.size(); ++i) {
    if (!bodies[i].isColored()) order.push_back(i);
  }
  for (size_t i = 0; i < bodies.size(); ++i) {
    if (bodies[i].isColored()) order.push_back(i);
  }
  struct Run {
    size_t first;
    Color4f color;
    std::vector<TopoDS_Shape> shapes;
    TopoDS_Shape fused;
  };
  std::vector<Run> runs;
  for (const size_t i : order) {
    if (!runs.empty() && sameColor(runs.back().color, bodies[i].color)) {
      runs.back().shapes.push_back(bodies[i].shape);
    } else {
      runs.push_back({i, bodies[i].color, {bodies[i].shape}, {}});
    }
  }
  for (auto& run : runs) run.fused = checkedFuse(run.shapes, dim);

  // Each run's surviving piece, by run index; a run whose cut fails is
  // merged with the runs it collided with into one uncolored body, so
  // the color loss stays local to that collision.
  struct Kept {
    size_t first;
    OcctBody body;
    std::vector<size_t> runs;  // runs this body accounts for
  };
  std::vector<Kept> kept;
  std::vector<size_t> higher;
  for (size_t r = runs.size(); r-- > 0;) {
    const auto& run = runs[r];
    // Only bodies that actually share material get cut away. A body that
    // merely touches a higher one (a part sitting in the cavity cut for
    // it, a difference next to the shape it was cut with) keeps its
    // color untouched, and OCCT is spared a cut along a coincident
    // surface, which it can answer with nothing.
    const double own = OcctBoolean::extent(run.fused, dim);
    std::vector<size_t> claimants;
    std::vector<TopoDS_Shape> tools;
    for (const size_t h : higher) {
      auto shared = checkedCommon({run.fused}, {runs[h].fused}, dim);
      if (!shared.IsNull() && OcctBoolean::extent(shared, dim) > 1e-6 * std::max(own, 1.0)) {
        claimants.push_back(h);
        tools.push_back(runs[h].fused);
      }
    }
    higher.push_back(r);
    if (claimants.empty()) {
      kept.push_back({run.first, {run.fused, run.color}, {r}});
      continue;
    }
    // What the cut must leave: the run minus everything it shares.
    double expected = own;
    if (auto shared = checkedCommon({run.fused}, tools, dim); !shared.IsNull()) {
      expected = own - OcctBoolean::extent(shared, dim);
    }
    auto piece = checkedCut({run.fused}, tools, dim);
    const double got = piece.IsNull() ? 0.0 : OcctBoolean::extent(piece, dim);
    if (std::fabs(got - expected) <= 1e-5 * std::max(own, 1.0)) {
      if (got > EXTENT_EPS) kept.push_back({run.first, {piece, run.color}, {r}});
      continue;
    }
    // The cut is not the answer. Merge this run with the runs it collided
    // with into one uncolored body; every other body keeps its color.
    std::vector<TopoDS_Shape> group{run.fused};
    std::vector<size_t> members{r};
    size_t firstIndex = run.first;
    for (auto it = kept.begin(); it != kept.end();) {
      bool collides = false;
      for (const size_t m : it->runs) {
        if (std::find(claimants.begin(), claimants.end(), m) != claimants.end()) collides = true;
      }
      if (collides) {
        for (const size_t m : it->runs) {
          group.push_back(runs[m].fused);
          members.push_back(m);
        }
        firstIndex = std::min(firstIndex, it->first);
        it = kept.erase(it);
      } else {
        ++it;
      }
    }
    OcctBridge::warn("STEP export: a boolean between colored bodies did not come back right; " +
                     std::to_string(members.size()) +
                     " bodies are merged and exported without color, the rest keep theirs");
    kept.push_back({firstIndex, {checkedFuse(group, dim), Color4f()}, members});
  }
  std::sort(kept.begin(), kept.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
  double total = 0;
  for (const auto& k : kept) total += OcctBoolean::extent(k.body.shape, dim);
  if (std::getenv("OPENSCAD_OCCT_DEBUG")) {
    std::string msg = "partition: fused=" + std::to_string(fusedExtent) +
                      " naive=" + std::to_string(naive) + " total=" + std::to_string(total) + " pieces:";
    for (const auto& k : kept) {
      msg += " " + std::to_string(OcctBoolean::extent(k.body.shape, dim)) + "(" +
             std::to_string(OcctBoolean::piecesOf(k.body.shape, dim).size()) + ")";
    }
    msg += " runs:";
    for (const auto& run : runs) msg += " " + std::to_string(OcctBoolean::extent(run.fused, dim));
    OcctBridge::warn(msg);
  }
  if (std::fabs(total - fusedExtent) > 1e-6 * std::max(fusedExtent, 1.0)) {
    OcctBridge::warn(
      "STEP export: partitioning a union by color lost material; the union is exported without "
      "its colors instead");
    return single(fused, dim, Color4f());
  }
  for (auto& k : kept) result.bodies.push_back(std::move(k.body));
  return result;
}

OcctGeometry OcctBuilder::differenceOf(const AbstractNode& node, std::vector<OcctGeometry> children)
{
  const VerificationScope scope(unverified_);
  auto result = differenceOfChecked(node, std::move(children));
  if (!scope.failed()) return result;
  const Color4f color = result.isEmpty() ? Color4f() : result.bodies.front().color;
  return meshFallback(node, color, UNVERIFIED);
}

OcctGeometry OcctBuilder::differenceOfChecked(const AbstractNode& node,
                                              std::vector<OcctGeometry> children)
{
  if (children.empty()) return {};
  // The first child is the minuend even when empty: with it gone there is
  // nothing to cut from, and the first subtrahend must not be promoted.
  OcctGeometry minuend = std::move(children.front());
  if (minuend.isEmpty()) return {};
  const unsigned int dim = minuend.dim;
  std::vector<TopoDS_Shape> tools;
  bool mixed = false;
  for (size_t i = 1; i < children.size(); ++i) {
    if (children[i].isEmpty()) continue;
    if (children[i].dim != dim) {
      mixed = true;
      continue;
    }
    for (const auto& body : children[i].bodies) tools.push_back(body.shape);
  }
  if (mixed) warn(node, "Mixing 2D and 3D objects is not supported");
  if (tools.empty()) return minuend;

  OcctGeometry result;
  result.dim = dim;
  bool oneColor = true;
  for (size_t i = 1; i < minuend.bodies.size(); ++i) {
    if (!sameColor(minuend.bodies[i].color, minuend.bodies[0].color)) oneColor = false;
  }
  if (oneColor) {
    auto shape = checkedCut(shapesOf(minuend), tools, dim);
    return single(shape, dim, minuend.bodies[0].color);
  }
  for (const auto& body : minuend.bodies) {
    auto shape = checkedCut({body.shape}, tools, dim);
    if (!shape.IsNull()) result.bodies.push_back({shape, body.color});
  }
  if (result.bodies.empty()) return {};
  return result;
}

OcctGeometry OcctBuilder::intersectionOf(const AbstractNode& node, std::vector<OcctGeometry> children)
{
  const VerificationScope scope(unverified_);
  auto result = intersectionOfChecked(node, std::move(children));
  if (!scope.failed()) return result;
  const Color4f color = result.isEmpty() ? Color4f() : result.bodies.front().color;
  return meshFallback(node, color, UNVERIFIED);
}

OcctGeometry OcctBuilder::intersectionOfChecked(const AbstractNode& node,
                                                std::vector<OcctGeometry> children)
{
  if (children.empty()) return {};
  // Any empty operand empties the intersection.
  for (const auto& child : children) {
    if (child.isEmpty()) return {};
  }
  const unsigned int dim = children.front().dim;
  for (const auto& child : children) {
    if (child.dim != dim) {
      warn(node, "Mixing 2D and 3D objects is not supported");
      return {};
    }
  }
  if (children.size() == 1) return children.front();

  OcctGeometry result;
  result.dim = dim;
  for (const auto& body : children.front().bodies) {
    TopoDS_Shape current = body.shape;
    for (size_t i = 1; i < children.size() && !current.IsNull(); ++i) {
      current = checkedCommon({current}, shapesOf(children[i]), dim);
    }
    if (!current.IsNull()) result.bodies.push_back({current, body.color});
  }
  if (result.bodies.empty()) return {};
  return result;
}

// --- leaves ----------------------------------------------------------------

OcctGeometry OcctBuilder::cube(const CubeNode& node, const Color4f& color)
{
  if (!(node.x > 0 && node.y > 0 && node.z > 0)) return {};
  TopoDS_Shape box = BRepPrimAPI_MakeBox(node.x, node.y, node.z).Shape();
  if (node.center) box = translated(box, -node.x / 2, -node.y / 2, -node.z / 2);
  return single(box, 3, color);
}

OcctGeometry OcctBuilder::sphere(const SphereNode& node, const Color4f& color)
{
  if (!(node.r > 0)) return {};
  if (const auto count = facetCount(node.discretizer, node.r)) {
    // OpenSCAD's sphere is a stack of rings; let it build that mesh.
    return meshFallback(node, color, "faceted sphere, " + std::to_string(*count) + " segments");
  }
  return single(BRepPrimAPI_MakeSphere(node.r).Shape(), 3, color);
}

OcctGeometry OcctBuilder::cylinder(const CylinderNode& node, const Color4f& color)
{
  if (!(node.h > 0) || (!(node.r1 > 0) && !(node.r2 > 0))) return {};
  const double r1 = std::max(node.r1, 0.0);
  const double r2 = std::max(node.r2, 0.0);
  const double z0 = node.center ? -node.h / 2 : 0.0;
  TopoDS_Shape shape;
  if (const auto faceted = facetCount(node.discretizer, std::max(r1, r2))) {
    const int count = *faceted;
    std::vector<Vector3d> points;
    std::vector<std::vector<size_t>> faces;
    const auto bottom = r1 > 0 ? ngon(r1, count) : std::vector<Vector2d>{};
    const auto top = r2 > 0 ? ngon(r2, count) : std::vector<Vector2d>{};
    if (bottom.empty()) {
      points.emplace_back(0, 0, z0);
    } else {
      for (const auto& p : bottom) points.emplace_back(p[0], p[1], z0);
    }
    const size_t nb = points.size();
    if (top.empty()) {
      points.emplace_back(0, 0, z0 + node.h);
    } else {
      for (const auto& p : top) points.emplace_back(p[0], p[1], z0 + node.h);
    }
    const size_t nt = points.size() - nb;
    if (nb > 1) {
      std::vector<size_t> face;
      for (size_t i = nb; i-- > 0;) face.push_back(i);
      faces.push_back(face);
    }
    if (nt > 1) {
      std::vector<size_t> face;
      for (size_t i = 0; i < nt; ++i) face.push_back(nb + i);
      faces.push_back(face);
    }
    for (size_t i = 0; i < static_cast<size_t>(count); ++i) {
      const size_t j = (i + 1) % count;
      if (nb > 1 && nt > 1) {
        faces.push_back({i, j, nb + j, nb + i});
      } else if (nb > 1) {
        faces.push_back({i, j, nb});
      } else {
        faces.push_back({0, 1 + j, 1 + i});
      }
    }
    bool nonPlanar = false;
    int freeEdges = 0;
    shape = OcctMesh::solidsFromFaces(points, faces, nonPlanar, freeEdges);
  } else if (std::fabs(r1 - r2) < 1e-12) {
    shape = BRepPrimAPI_MakeCylinder(r1, node.h).Shape();
    if (node.center) shape = translated(shape, 0, 0, z0);
  } else {
    shape = BRepPrimAPI_MakeCone(r1, r2, node.h).Shape();
    if (node.center) shape = translated(shape, 0, 0, z0);
  }
  return single(shape, 3, color);
}

OcctGeometry OcctBuilder::polyhedron(const PolyhedronNode& node, const Color4f& color)
{
  if (node.points.empty() || node.faces.empty()) return {};
  std::vector<std::vector<size_t>> faces;
  faces.reserve(node.faces.size());
  for (const auto& face : node.faces) {
    std::vector<size_t> indices;
    for (const auto index : face) {
      if (index >= 0) indices.push_back(static_cast<size_t>(index));
    }
    faces.push_back(std::move(indices));
  }
  bool nonPlanar = false;
  int freeEdges = 0;
  auto shape = OcctMesh::solidsFromFaces(node.points, faces, nonPlanar, freeEdges);
  if (nonPlanar) {
    return meshFallback(node, color, "polyhedron() has a non-planar face, which OpenSCAD tessellates");
  }
  if (freeEdges > 0) {
    warn(node, "polyhedron() leaves " + std::to_string(freeEdges) +
                 " free edge(s), so it is not a valid 2-manifold and contributes nothing");
    return {};
  }
  if (shape.IsNull()) {
    return meshFallback(node, color, "polyhedron() could not be sewn into a solid");
  }
  return single(orientedOutward(shape), 3, color);
}

OcctGeometry OcctBuilder::square(const SquareNode& node, const Color4f& color)
{
  if (!(node.x > 0 && node.y > 0)) return {};
  const double x0 = node.center ? -node.x / 2 : 0.0;
  const double y0 = node.center ? -node.y / 2 : 0.0;
  const std::vector<Vector3d> ring{
    {x0, y0, 0}, {x0 + node.x, y0, 0}, {x0 + node.x, y0 + node.y, 0}, {x0, y0 + node.y, 0}};
  return single(OcctMesh::faceFromRing(ring), 2, color);
}

OcctGeometry OcctBuilder::circle(const CircleNode& node, const Color4f& color)
{
  if (!(node.r > 0)) return {};
  if (const auto count = facetCount(node.discretizer, node.r)) {
    std::vector<Vector3d> ring;
    for (const auto& p : ngon(node.r, *count)) ring.emplace_back(p[0], p[1], 0);
    return single(OcctMesh::faceFromRing(ring), 2, color);
  }
  const gp_Circ circ(gp_Ax2(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1)), node.r);
  const auto edge = BRepBuilderAPI_MakeEdge(circ).Edge();
  const auto wire = BRepBuilderAPI_MakeWire(edge).Wire();
  return single(BRepBuilderAPI_MakeFace(wire, true).Face(), 2, color);
}

OcctGeometry OcctBuilder::polygon(const PolygonNode& node, const Color4f& color)
{
  auto geometry = node.createGeometry();
  const auto *poly = dynamic_cast<const Polygon2d *>(geometry.get());
  if (!poly || poly->isEmpty()) return {};
  return single(OcctMesh::facesFromPolygon2d(*poly), 2, color);
}

OcctGeometry OcctBuilder::text(const TextNode& node, const Color4f& color)
{
  auto polygons = node.createPolygonList();
  auto merged = ClipperUtils::apply(polygons, Clipper2Lib::ClipType::Union);
  if (!merged || merged->isEmpty()) return {};
  return single(OcctMesh::facesFromPolygon2d(*merged), 2, color);
}

// --- transforms and extrusions ---------------------------------------------

OcctGeometry OcctBuilder::transform(const TransformNode& node, const Color4f& inherited)
{
  auto geometry = unionOf(node, buildChildren(node, inherited));
  if (geometry.isEmpty()) return {};
  if (matrix_contains_infinity(node.matrix) || matrix_contains_nan(node.matrix)) {
    warn(node, "Transformation matrix contains Not-a-Number and/or Infinity - removing object.");
    return {};
  }
  Transform3d m = node.matrix;
  if (geometry.dim == 2) {
    // 2D geometry has no z: it feels only the in-plane part of the matrix
    // and never leaves the XY plane.
    Transform3d flat = Transform3d::Identity();
    flat(0, 0) = m(0, 0);
    flat(0, 1) = m(0, 1);
    flat(0, 3) = m(0, 3);
    flat(1, 0) = m(1, 0);
    flat(1, 1) = m(1, 1);
    flat(1, 3) = m(1, 3);
    m = flat;
    if (std::fabs(m(0, 0) * m(1, 1) - m(0, 1) * m(1, 0)) < SINGULAR_EPS) {
      warn(node, "Transformation flattens the object - removing it.");
      return {};
    }
  } else if (std::fabs(m.matrix().block<3, 3>(0, 0).determinant()) < SINGULAR_EPS) {
    warn(node, "Transformation flattens the object - removing it.");
    return {};
  }

  const auto decomposed = decompose(m);
  OcctGeometry result;
  result.dim = geometry.dim;
  for (const auto& body : geometry.bodies) {
    TopoDS_Shape shape;
    try {
      if (decomposed) {
        const auto& r = decomposed->rotation;
        gp_Trsf rigid;
        rigid.SetValues(r[0][0], r[0][1], r[0][2], m(0, 3), r[1][0], r[1][1], r[1][2], m(1, 3), r[2][0],
                        r[2][1], r[2][2], m(2, 3));
        gp_Trsf trsf = rigid;
        if (std::fabs(decomposed->scale - 1.0) > RIGID_TOL) {
          gp_Trsf scaling;
          scaling.SetScale(gp_Pnt(0, 0, 0), decomposed->scale);
          trsf = rigid * scaling;
        }
        BRepBuilderAPI_Transform op(body.shape, trsf, true);
        if (op.IsDone()) shape = op.Shape();
      } else {
        gp_GTrsf gtrsf;
        gtrsf.SetVectorialPart(
          gp_Mat(m(0, 0), m(0, 1), m(0, 2), m(1, 0), m(1, 1), m(1, 2), m(2, 0), m(2, 1), m(2, 2)));
        gtrsf.SetTranslationPart(gp_XYZ(m(0, 3), m(1, 3), m(2, 3)));
        BRepBuilderAPI_GTransform op(body.shape, gtrsf, true);
        if (op.IsDone()) shape = op.Shape();
      }
    } catch (const Standard_Failure&) {
      shape.Nullify();
    }
    if (shape.IsNull()) {
      return meshFallback(node, inherited, "the transformation could not be applied to B-rep geometry");
    }
    if (geometry.dim == 3) shape = orientedOutward(shape);
    result.bodies.push_back({shape, body.color});
  }
  return result;
}

OcctGeometry OcctBuilder::flatChildren(const AbstractNode& node, const Color4f& inherited,
                                       unsigned int wantDim)
{
  auto geometry = unionOf(node, buildChildren(node, inherited));
  if (geometry.isEmpty()) return {};
  if (geometry.dim != wantDim) {
    warn(node, wantDim == 2 ? "Ignoring 3D child object for 2D operation"
                            : "Ignoring 2D child object for 3D operation");
    return {};
  }
  return geometry;
}

OcctGeometry OcctBuilder::linearExtrude(const LinearExtrudeNode& node, const Color4f& inherited)
{
  auto profile = flatChildren(node, inherited, 2);
  if (profile.isEmpty()) return {};
  if (node.height[2] <= 0) return {};
  if (node.has_twist && node.twist != 0) {
    return meshFallback(node, inherited, "linear_extrude() with twist has no exact B-rep form");
  }
  const Vector3d height = node.height;
  const Vector3d offset = node.center ? Vector3d(-height / 2.0) : Vector3d::Zero();
  const bool scaled = node.scale_x != 1.0 || node.scale_y != 1.0;

  OcctGeometry result;
  result.dim = 3;
  for (const auto& body : profile.bodies) {
    TopoDS_Shape shape;
    try {
      if (!scaled) {
        TopoDS_Shape base = translated(body.shape, offset[0], offset[1], offset[2]);
        BRepPrimAPI_MakePrism prism(base, gp_Vec(height[0], height[1], height[2]));
        if (prism.IsDone()) shape = prism.Shape();
      } else {
        if (!(node.scale_x > 0 && node.scale_y > 0)) {
          return meshFallback(node, inherited, "linear_extrude() scaling to zero");
        }
        // Loft every wire of every face to its scaled copy at the top:
        // outer wires make solids, hole wires are cut away again.
        std::vector<TopoDS_Shape> outers, holes;
        for (TopExp_Explorer faces(body.shape, TopAbs_FACE); faces.More(); faces.Next()) {
          const auto face = TopoDS::Face(faces.Current());
          const auto outer = ShapeAnalysis::OuterWire(face);
          for (TopExp_Explorer wires(face, TopAbs_WIRE); wires.More(); wires.Next()) {
            const auto wire = TopoDS::Wire(wires.Current());
            gp_GTrsf scale;
            scale.SetVectorialPart(gp_Mat(node.scale_x, 0, 0, 0, node.scale_y, 0, 0, 0, 1));
            BRepBuilderAPI_GTransform scaledOp(wire, scale, true);
            if (!scaledOp.IsDone()) return meshFallback(node, inherited, "linear_extrude() scale");
            auto top = translated(scaledOp.Shape(), height[0], height[1], height[2]);
            auto bottom = TopoDS_Shape(wire);
            BRepOffsetAPI_ThruSections loft(true, true);
            loft.AddWire(TopoDS::Wire(bottom));
            loft.AddWire(TopoDS::Wire(top));
            loft.Build();
            if (!loft.IsDone())
              return meshFallback(node, inherited, "linear_extrude() scale loft failed");
            (wire.IsSame(outer) ? outers : holes).push_back(loft.Shape());
          }
        }
        bool ok = true;
        shape = OcctBoolean::fuse(outers, 3, &ok);
        if (ok && !holes.empty()) shape = OcctBoolean::cut({shape}, holes, 3, &ok);
        if (!ok) shape.Nullify();
        if (!shape.IsNull()) shape = translated(shape, offset[0], offset[1], offset[2]);
      }
    } catch (const Standard_Failure&) {
      shape.Nullify();
    }
    if (shape.IsNull()) return meshFallback(node, inherited, "linear_extrude() failed in B-rep");
    result.bodies.push_back({orientedOutward(shape), body.color});
  }
  return result;
}

OcctGeometry OcctBuilder::rotateExtrude(const RotateExtrudeNode& node, const Color4f& inherited)
{
  auto profile = flatChildren(node, inherited, 2);
  if (profile.isEmpty()) return {};
  if (node.angle == 0) return {};

  // A tight box: the default one is padded by the shape's tolerance, which
  // reads as "crosses the axis" for a profile that merely touches it.
  Bnd_Box box;
  for (const auto& body : profile.bodies) BRepBndLib::AddOptimal(body.shape, box, false, false);
  double xmin, ymin, zmin, xmax, ymax, zmax;
  box.Get(xmin, ymin, zmin, xmax, ymax, zmax);
  const double axisTol = 1e-6 * std::max({std::fabs(xmin), std::fabs(xmax), 1.0});
  if (xmax > axisTol && xmin < -axisTol) {
    OcctBridge::error(node, tree_,
                      "Children of rotate_extrude() may not lie across the Y axis (Range of X coords "
                      "for all children [" +
                        std::to_string(xmin) + " : " + std::to_string(xmax) + "])");
    return {};
  }

  const double angle = std::clamp(node.angle, -360.0, 360.0);
  const bool full = std::fabs(angle) >= 360.0 - 1e-9;
  const double startDeg = angle < 0 ? node.start + angle : node.start;

  OcctGeometry result;
  result.dim = 3;
  for (const auto& body : profile.bodies) {
    TopoDS_Shape shape;
    try {
      // The XY profile becomes the XZ profile: (x, y) -> (x, 0, y).
      gp_Trsf toXZ;
      toXZ.SetRotation(gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(1, 0, 0)), M_PI / 2);
      gp_Trsf toStart;
      toStart.SetRotation(gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1)), startDeg * M_PI / 180.0);
      TopoDS_Shape placed = BRepBuilderAPI_Transform(body.shape, toStart * toXZ, true).Shape();
      const gp_Ax1 axis(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1));
      if (full) {
        BRepPrimAPI_MakeRevol revol(placed, axis);
        if (revol.IsDone()) shape = revol.Shape();
      } else {
        BRepPrimAPI_MakeRevol revol(placed, axis, std::fabs(angle) * M_PI / 180.0);
        if (revol.IsDone()) shape = revol.Shape();
      }
    } catch (const Standard_Failure&) {
      shape.Nullify();
    }
    if (shape.IsNull()) return meshFallback(node, inherited, "rotate_extrude() failed in B-rep");
    result.bodies.push_back({orientedOutward(shape), body.color});
  }
  return result;
}

// --- hull and minkowski ----------------------------------------------------

OcctGeometry OcctBuilder::hullOrMinkowski(const CgalAdvNode& node, const Color4f& inherited)
{
  const bool isHull = node.type == CgalAdvType::HULL;
  auto children = buildChildren(node, inherited);
  unsigned int dim = 0;
  bool mixed = false;
  auto kept = sameDimension(children, dim, mixed);
  if (mixed) warn(node, "Mixing 2D and 3D objects is not supported");
  // hull() and minkowski() of nothing are nothing, and asking OpenSCAD to
  // mesh an empty subtree is an error rather than a fallback.
  if (kept.empty()) return {};

  std::vector<TopoDS_Shape> shapes;
  Color4f color = inherited;
  for (const auto& child : kept) {
    for (const auto& body : child.bodies) {
      shapes.push_back(body.shape);
      if (!color.isValid() && body.isColored()) color = body.color;
    }
  }
  std::string rung;
  TopoDS_Shape result;
  if (isHull) {
    result = OcctHull::hull(OcctHull::components(shapes, dim), dim, rung);
  } else {
    result = OcctHull::minkowski(shapes, dim, rung);
  }
  if (!result.IsNull()) {
    if (dim == 3) result = orientedOutward(result);
    return single(result, dim, color);
  }
  if (isHull && dim == 3) {
    if (!evaluator_) evaluator_ = OcctBridge::makeEvaluator(tree_);
    const auto rendered = OcctBridge::render(*evaluator_, node);
    size_t triangles = 0;
    for (const auto& mesh : rendered.meshes) triangles += mesh.faces.size();
    if (triangles <= HULL_MESH_BUDGET) {
      return meshFallback(node, inherited, "children fit none of the closed-form hull cases", rendered);
    }
    result = OcctHull::hullOfTessellation(OcctHull::components(shapes, dim), HULL_FALLBACK_SEGMENTS);
    if (!result.IsNull()) {
      const std::string note =
        node.name() +
        "(): children fit none of the closed-form hull cases; built as the hull of "
        "the children tessellated at " +
        std::to_string(HULL_FALLBACK_SEGMENTS) + " segments. Surfaces in this region are polyhedral.";
      fallbacks_.push_back(note);
      warn(node, "STEP export: " + note);
      return single(orientedOutward(result), dim, color);
    }
    // The tessellation could not be built; the render already in hand is
    // still the right fallback, whatever its size.
    return meshFallback(node, inherited, "children fit none of the closed-form hull cases", rendered);
  }
  return meshFallback(node, inherited,
                      isHull ? "children fit none of the closed-form hull cases"
                             : "no operand is a sphere or circle at the origin");
}

// --- mesh fallback ---------------------------------------------------------

OcctGeometry OcctBuilder::meshFallback(const AbstractNode& node, const Color4f& color,
                                       const std::string& why)
{
  if (!evaluator_) evaluator_ = OcctBridge::makeEvaluator(tree_);
  return meshFallback(node, color, why, OcctBridge::render(*evaluator_, node));
}

OcctGeometry OcctBuilder::meshFallback(const AbstractNode& node, const Color4f& color,
                                       const std::string& why, const OcctBridge::Rendered& rendered)
{
  if (rendered.dim == 0) return {};

  OcctGeometry result;
  result.dim = rendered.dim;
  std::string meshWhy;
  for (const auto& mesh : rendered.meshes) {
    auto shape = OcctMesh::solidsFromMesh(mesh, meshWhy);
    if (!shape.IsNull()) result.bodies.push_back({shape, color});
  }
  for (const auto& polygon : rendered.polygons) {
    Polygon2d poly;
    for (const auto& outline : polygon.outlines) {
      Outline2d o;
      o.vertices = outline.vertices;
      o.positive = outline.positive;
      poly.addOutline(o);
    }
    poly.setSanitized(true);
    auto shape = OcctMesh::facesFromPolygon2d(poly);
    if (!shape.IsNull()) result.bodies.push_back({shape, color});
  }
  if (result.isEmpty()) {
    warn(node, node.name() + "() rendered to a mesh that could not be made into a solid: " + meshWhy +
                 "; it contributes nothing (" + why + ")");
    return {};
  }
  const std::string note = node.name() + "(): " + why +
                           "; rendered to a mesh via OpenSCAD's own evaluator. Surfaces in this "
                           "region are not analytic.";
  fallbacks_.push_back(note);
  warn(node, "STEP export: " + note);
  return result;
}
