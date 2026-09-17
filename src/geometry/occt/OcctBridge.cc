#include "geometry/occt/OcctBridge.h"

#include <memory>
#include <string>
#include <vector>

#include "core/ModuleInstantiation.h"
#include "core/Tree.h"
#include "core/node.h"
#include "geometry/Geometry.h"
#include "geometry/GeometryEvaluator.h"
#include "geometry/PolySet.h"
#include "geometry/Polygon2d.h"
#include "utils/printutils.h"

namespace OcctBridge {

std::unique_ptr<GeometryEvaluator> makeEvaluator(const Tree& tree)
{
  return std::make_unique<GeometryEvaluator>(tree);
}

namespace {

void collect(const std::shared_ptr<const Geometry>& geometry, Rendered& out)
{
  if (!geometry || geometry->isEmpty()) return;
  if (const auto list = std::dynamic_pointer_cast<const GeometryList>(geometry)) {
    for (const auto& item : list->flatten()) collect(item.second, out);
    return;
  }
  if (const auto ps = std::dynamic_pointer_cast<const PolySet>(geometry)) {
    if (ps->getDimension() != 3) return;
    MeshData mesh;
    mesh.vertices = ps->vertices;
    for (const auto& polygon : ps->indices) {
      std::vector<size_t> face;
      for (const auto index : polygon) face.push_back(static_cast<size_t>(index));
      mesh.faces.push_back(std::move(face));
    }
    out.dim = 3;
    out.meshes.push_back(std::move(mesh));
    return;
  }
  if (const auto poly = std::dynamic_pointer_cast<const Polygon2d>(geometry)) {
    PolygonData data;
    for (const auto& outline : poly->outlines()) {
      data.outlines.push_back({outline.vertices, outline.positive});
    }
    out.dim = 2;
    out.polygons.push_back(std::move(data));
  }
}

}  // namespace

Rendered render(GeometryEvaluator& evaluator, const AbstractNode& node)
{
  Rendered out;
  collect(evaluator.evaluateGeometry(node, false), out);
  return out;
}

void warn(const std::string& message)
{
  LOG(message_group::Warning, "%1$s", message);
}

void warn(const AbstractNode& node, const Tree& tree, const std::string& message)
{
  LOG(message_group::Warning, node.modinst->location(), tree.getDocumentPath(), "%1$s", message);
}

void error(const std::string& message)
{
  LOG(message_group::Error, "%1$s", message);
}

void error(const AbstractNode& node, const Tree& tree, const std::string& message)
{
  LOG(message_group::Error, node.modinst->location(), tree.getDocumentPath(), "%1$s", message);
}

}  // namespace OcctBridge
