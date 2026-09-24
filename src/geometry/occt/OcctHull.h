#pragma once

#include <TopoDS_Shape.hxx>
#include <string>
#include <vector>

// Closed-form hull() and minkowski(): every case with an exact B-rep
// answer, ported from solid123d. Each rung classifies the *built*
// children (module-heavy code wraps a primitive in layers of group and
// multmatrix bookkeeping that the builder has already resolved), builds
// the result without booleans where it can, and self-checks against a
// closed form. A null shape means no rung applies and the caller should
// fall back to OpenSCAD's mesh.
namespace OcctHull {

// Each child exploded into its independent solids (3D) or faces (2D):
// hull(A union B) == hull(A, B), and a module body arrives as one
// pre-fused group.
std::vector<TopoDS_Shape> components(const std::vector<TopoDS_Shape>& children, unsigned int dim);

// hull():
// - N equal-radius spheres: offset(convex hull of centers, r); collinear
//   centers become a capsule.
// - Exactly two spheres of any radii: two spherical caps sewn to the
//   external tangent cone.
// - N equal-radius parallel cylinders sharing one axial span: the 2D case
//   of the above, extruded.
// - All-polyhedral children: the convex hull of their vertices.
// - Two discs in the plane: arcs joined by the external tangent lines.
// - Straight-edged planar children: the 2D convex hull of their vertices.
TopoDS_Shape hull(const std::vector<TopoDS_Shape>& components, unsigned int dim, std::string& rung);

// The last resort for a 3D hull no closed form covers: tessellate each
// component at `segments` per full turn, and hull the vertices. This is
// what OpenSCAD itself does, at a fineness the exporter chooses rather than
// the model's $fn: above the facet threshold that $fn buys nothing the
// exact tier honours, while a fallback rendered at $fn=180 carried 17k
// triangles that 2k describe to 0.1% of the volume, and OpenCASCADE's
// booleans are badly superlinear in facet count. Polyhedral, so callers
// must report it as a fallback. Null if nothing could be tessellated.
TopoDS_Shape hullOfTessellation(const std::vector<TopoDS_Shape>& components, int segments);

// minkowski() of exactly two operands where one is a ball at the origin
// (a sphere or circle primitive, or a many-vertex polyhedron/polygon whose
// vertices are equidistant from the origin): an offset of the other.
TopoDS_Shape minkowski(const std::vector<TopoDS_Shape>& children, unsigned int dim, std::string& rung);

}  // namespace OcctHull
