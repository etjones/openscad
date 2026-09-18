#pragma once

#include <TopoDS_Shape.hxx>
#include <string>
#include <vector>

// Boolean operations on OCCT shapes, with the checks OCCT itself does not
// make. Every function returns a null TopoDS_Shape for an empty result.
//
// The invariants and retry ladders here are ported from solid123d, where
// they were arrived at by converting a corpus of ~76k OpenSCAD files:
// OCCT's booleans occasionally return a *valid* shape that is not the
// answer (a fuse missing most of its material, a cut that kept the
// material it was cutting away, a union in more pieces than it was
// given). Each is caught by a statement about the answer -- its volume
// against the bounds the inputs imply, its piece count, whether material
// remains inside a cutting tool -- and retried with a fuzzy tolerance
// sized as a fraction of the model.
namespace OcctBoolean {

// The pieces a boolean is defined on: solids for 3D geometry, faces for 2D.
// OCCT mishandles a compound of touching bodies as one operand, so every
// operand is decomposed into these before an operation.
std::vector<TopoDS_Shape> piecesOf(const TopoDS_Shape& shape, unsigned int dim);

// Volume for 3D geometry, area for 2D: the quantity a boolean must not
// silently leave unchanged.
double extent(const TopoDS_Shape& shape, unsigned int dim);

double diagonal(const std::vector<TopoDS_Shape>& shapes);

TopoDS_Shape makeCompound(const std::vector<TopoDS_Shape>& shapes);

// Does `result` still hold every one of `operands`? A union contains
// everything it was given, so a body of the input that no body of the
// result contains means the fuse dropped it. Sampled at up to eight
// interior points per operand, so it is a detector, not a proof: an
// operand partly swallowed can still pass. 2D geometry always passes.
bool keepsOperands(const std::vector<TopoDS_Shape>& operands, const TopoDS_Shape& result,
                   unsigned int dim);

// Does `result` still hold the material of `args` that no tool covers? A
// cut removes only what its tools cover, so an uncovered interior point
// missing from the result means the cut took too much. Same sampling
// caveat as above.
bool keepsUncoveredMaterial(const std::vector<TopoDS_Shape>& args,
                            const std::vector<TopoDS_Shape>& tools, const TopoDS_Shape& result,
                            unsigned int dim);

// N-ary fuse, checked: the result must hold at least the largest input
// and at most their sum, must not come back in more pieces than it was
// given, must have no open shells, and must still contain every operand
// (sampled at several interior points each). On failure, retried with
// fuzzy tolerances and then one operand at a time; if none helps the
// operands are returned unjoined, `verified` (when given) is cleared and
// a warning is logged. A caller that can render the node another way
// should do so when `verified` comes back false.
TopoDS_Shape fuse(const std::vector<TopoDS_Shape>& operands, unsigned int dim, bool *verified = nullptr);

// N-ary cut, checked: no result body may sit inside a tool, no open shells.
TopoDS_Shape cut(const std::vector<TopoDS_Shape>& args, const std::vector<TopoDS_Shape>& tools,
                 unsigned int dim, bool *verified = nullptr);

TopoDS_Shape common(const std::vector<TopoDS_Shape>& args, const std::vector<TopoDS_Shape>& tools,
                    unsigned int dim, bool *verified = nullptr);

// ShapeUpgrade_UnifySameDomain, adopted only when it conserves the
// shape's extent (it can delete faces crossed by a surface's seam).
TopoDS_Shape unify(const TopoDS_Shape& shape, unsigned int dim);

// Splits faces that wrap all the way around a periodic surface, so no
// face in the result carries a seam. Written STEP survives the trip
// through other CAD systems that way; see the comment on the definition.
// Adopted only when it conserves the shape's extent.
TopoDS_Shape splitClosedFaces(const TopoDS_Shape& shape, unsigned int dim);

}  // namespace OcctBoolean
