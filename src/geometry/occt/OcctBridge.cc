#include "geometry/occt/OcctBridge.h"

#include <cstdlib>
#include <fstream>
#include <memory>
#include <string>
#include <vector>

#include "core/ModuleInstantiation.h"
#include "core/Tree.h"
#include "core/node.h"
#include "geometry/Geometry.h"
#include "geometry/GeometryEvaluator.h"
#include "geometry/PolySet.h"
#include "geometry/PolySetUtils.h"
#include "geometry/Polygon2d.h"
#include "utils/printutils.h"

namespace OcctBridge {

std::string idString(const Tree& tree, const AbstractNode& node)
{
  return tree.getIdString(node);
}

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
  if (auto ps = std::dynamic_pointer_cast<const PolySet>(geometry)) {
    if (ps->getDimension() != 3) return;
    // Triangles sew reliably; a polygon merged from coplanar triangles can
    // carry a vertex in the middle of a neighbor's edge, which leaves free
    // edges after sewing. Coplanar triangles are merged back afterwards.
    if (!ps->isTriangular()) ps = PolySetUtils::tessellate_faces(*ps);
    MeshData mesh;
    mesh.vertices = ps->vertices;
    for (const auto& polygon : ps->indices) {
      std::vector<size_t> face;
      for (const auto index : polygon) face.push_back(static_cast<size_t>(index));
      mesh.faces.push_back(std::move(face));
    }
    if (const char *dump = std::getenv("OPENSCAD_OCCT_DUMP_MESH")) {
      std::ofstream off(dump);
      off << "OFF\n" << mesh.vertices.size() << " " << mesh.faces.size() << " 0\n";
      off.precision(17);
      for (const auto& v : mesh.vertices) off << v[0] << " " << v[1] << " " << v[2] << "\n";
      for (const auto& f : mesh.faces) {
        off << f.size();
        for (const auto i : f) off << " " << i;
        off << "\n";
      }
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
