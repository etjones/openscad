#pragma once

#include <TopoDS_Shape.hxx>
#include <string>
#include <vector>

#include "geometry/linalg.h"
#include "geometry/occt/OcctBridge.h"

class Polygon2d;

// Conversions between OpenSCAD's own geometry and OCCT B-rep.
namespace OcctMesh {

// Planar faces from the polygon's outlines: one face per positive outline,
// with the negative outlines it contains as holes. The polygon is
// sanitized first when it is not already. Returns a compound of faces, or
// a null shape for an empty polygon.
TopoDS_Shape facesFromPolygon2d(const Polygon2d& polygon);

// A polygon face in the XY plane from a list of vertices. Null on failure
// (fewer than three distinct points, or a non-planar ring).
TopoDS_Shape faceFromRing(const std::vector<Vector3d>& ring);

// Solids from a closed polygon mesh: faces are sewn into shells, each shell
// becomes a solid oriented outward, shells inside another become its
// cavities, and coplanar facets are merged. Returns a compound of solids,
// or a null shape if the mesh does not close.
TopoDS_Shape solidsFromMesh(const OcctBridge::MeshData& mesh, std::string& why);

// Solids from explicit points and faces (the polyhedron() primitive).
// Faces follow OpenSCAD's convention: vertices ordered clockwise seen from
// outside. Sets nonPlanar when a face could not be built as a plane.
TopoDS_Shape solidsFromFaces(const std::vector<Vector3d>& points,
                             const std::vector<std::vector<size_t>>& faces, bool& nonPlanar,
                             int& freeEdges);

}  // namespace OcctMesh
