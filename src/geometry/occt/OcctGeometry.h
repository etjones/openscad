#pragma once

#include <TopoDS_Shape.hxx>
#include <vector>

#include "geometry/linalg.h"

// One body of B-rep geometry: a solid (3D) or a face (2D), or a compound
// of several, together with the color() it carries. An invalid Color4f
// means "no color was assigned".
struct OcctBody {
  TopoDS_Shape shape;
  Color4f color;

  [[nodiscard]] bool isColored() const { return color.isValid(); }
};

// What the OCCT evaluator produces for a node: nothing (dim 0), 2D faces
// or 3D solids. Bodies are kept apart instead of fused into one compound
// so a color assignment survives every operation above it.
struct OcctGeometry {
  unsigned int dim = 0;
  std::vector<OcctBody> bodies;

  [[nodiscard]] bool isEmpty() const { return bodies.empty(); }
  [[nodiscard]] bool hasColor() const
  {
    for (const auto& body : bodies) {
      if (body.isColored()) return true;
    }
    return false;
  }
};
