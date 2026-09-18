#pragma once

#include <optional>
#include <string>
#include <vector>

#include "core/node.h"
#include "geometry/linalg.h"
#include "geometry/occt/OcctGeometry.h"

class Tree;
class GeometryEvaluator;
class Polygon2d;
class CubeNode;
class SphereNode;
class CylinderNode;
class PolyhedronNode;
class SquareNode;
class CircleNode;
class PolygonNode;
class TextNode;
class TransformNode;
class ColorNode;
class CsgOpNode;
class LinearExtrudeNode;
class RotateExtrudeNode;
class CgalAdvNode;

// Evaluates an OpenSCAD node tree into OCCT B-rep geometry.
//
// Like the CSG exporter, this works from the evaluated node tree rather
// than from rendered meshes: the tree still holds analytic primitives
// (a cylinder is a cylinder, not a ring of quads), which is what makes a
// real STEP file possible. Nodes with no B-rep equivalent (hull,
// minkowski, projection, surface, mesh imports, offset, twisted or
// hollow scaled extrusions) are rendered by the regular GeometryEvaluator
// and their mesh is sewn into a solid, so nothing fails to export; that
// region merely loses its analytic surfaces.
class OcctBuilder
{
public:
  // facetThreshold: a $fn below this is honored as real polygonal
  // geometry (circle(r, $fn=6) is a hexagon); at or above it, the curve
  // is exact. 0 disables faceting.
  OcctBuilder(const Tree& tree, int facetThreshold = 20);
  ~OcctBuilder();

  OcctGeometry build(const AbstractNode& node);

  // Human-readable notes about subtrees that took the mesh path.
  [[nodiscard]] const std::vector<std::string>& fallbacks() const { return fallbacks_; }

private:
  OcctGeometry buildNode(const AbstractNode& node, const Color4f& inherited);
  std::vector<OcctGeometry> buildChildren(const AbstractNode& node, const Color4f& inherited);

  OcctGeometry unionOf(const AbstractNode& node, std::vector<OcctGeometry> children);
  OcctGeometry differenceOf(const AbstractNode& node, std::vector<OcctGeometry> children);
  OcctGeometry intersectionOf(const AbstractNode& node, std::vector<OcctGeometry> children);
  OcctGeometry partitionedUnion(std::vector<OcctGeometry> children, unsigned int dim);
  OcctGeometry unionOfChecked(const AbstractNode& node, std::vector<OcctGeometry> children);
  OcctGeometry differenceOfChecked(const AbstractNode& node, std::vector<OcctGeometry> children);
  OcctGeometry intersectionOfChecked(const AbstractNode& node, std::vector<OcctGeometry> children);

  OcctGeometry cube(const CubeNode& node, const Color4f& color);
  OcctGeometry sphere(const SphereNode& node, const Color4f& color);
  OcctGeometry cylinder(const CylinderNode& node, const Color4f& color);
  OcctGeometry polyhedron(const PolyhedronNode& node, const Color4f& color);
  OcctGeometry square(const SquareNode& node, const Color4f& color);
  OcctGeometry circle(const CircleNode& node, const Color4f& color);
  OcctGeometry polygon(const PolygonNode& node, const Color4f& color);
  OcctGeometry text(const TextNode& node, const Color4f& color);
  OcctGeometry transform(const TransformNode& node, const Color4f& inherited);
  OcctGeometry linearExtrude(const LinearExtrudeNode& node, const Color4f& inherited);
  OcctGeometry rotateExtrude(const RotateExtrudeNode& node, const Color4f& inherited);
  OcctGeometry hullOrMinkowski(const CgalAdvNode& node, const Color4f& inherited);

  // Render the subtree with OpenSCAD's own evaluator and sew the mesh.
  OcctGeometry meshFallback(const AbstractNode& node, const Color4f& color, const std::string& why);

  // The union of a node's 2D children as one Polygon2d-derived compound of faces.
  OcctGeometry flatChildren(const AbstractNode& node, const Color4f& inherited, unsigned int wantDim);

  [[nodiscard]] std::optional<int> facetCount(const class CurveDiscretizer& discretizer, double r) const;
  void warn(const AbstractNode& node, const std::string& message) const;

  const Tree& tree_;
  int facetThreshold_;
  std::unique_ptr<GeometryEvaluator> evaluator_;
  std::vector<std::string> fallbacks_;
  // Set when a boolean below the node being built could not be verified;
  // the node is then rendered to a mesh instead of trusting the result.
  bool unverified_ = false;

  TopoDS_Shape checkedFuse(const std::vector<TopoDS_Shape>& operands, unsigned int dim);
  TopoDS_Shape checkedCut(const std::vector<TopoDS_Shape>& args, const std::vector<TopoDS_Shape>& tools,
                          unsigned int dim);
  TopoDS_Shape checkedCommon(const std::vector<TopoDS_Shape>& args,
                             const std::vector<TopoDS_Shape>& tools, unsigned int dim);
};
