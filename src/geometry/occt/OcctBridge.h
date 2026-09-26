#pragma once

#include <memory>
#include <string>
#include <vector>

#include "geometry/linalg.h"

// The parts of OpenSCAD the OCCT evaluator needs that cannot be included
// next to OpenCASCADE headers: OCCT declares a global class `Message`,
// and so does utils/printutils.h, which core/Tree.h and geometry/PolySet.h
// pull in. Everything here is implemented in a translation unit that
// includes no OCCT header.
class AbstractNode;
class Tree;
class GeometryEvaluator;

namespace OcctBridge {

struct MeshData {
  std::vector<Vector3d> vertices;
  std::vector<std::vector<size_t>> faces;
};

struct Outline {
  VectorOfVector2d vertices;
  bool positive = true;
};

struct PolygonData {
  std::vector<Outline> outlines;
};

// What OpenSCAD's own evaluator renders for a subtree, as plain data.
struct Rendered {
  unsigned int dim = 0;
  std::vector<MeshData> meshes;       // dim 3
  std::vector<PolygonData> polygons;  // dim 2
};

std::unique_ptr<GeometryEvaluator> makeEvaluator(const Tree& tree);

// The tree's own identifier for a subtree, memoised by the tree. Equal
// strings mean equal geometry, which is how repeated subtrees are built
// once.
std::string idString(const Tree& tree, const AbstractNode& node);
Rendered render(GeometryEvaluator& evaluator, const AbstractNode& node);

void warn(const std::string& message);
void warn(const AbstractNode& node, const Tree& tree, const std::string& message);
void error(const std::string& message);
void error(const AbstractNode& node, const Tree& tree, const std::string& message);

}  // namespace OcctBridge
