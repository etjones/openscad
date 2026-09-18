#pragma once

#include <TopoDS_Shape.hxx>
#include <vector>

// hull() of one solid of revolution repeated by translation perpendicular
// to a shared axis: the ``hull() cornercopy(...)`` idiom (tapered pads,
// stacked bevels, filleted cavity posts, turned legs).
//
// hull(identical translates of X) == conv(centers) (+) conv(X). For X a
// solid of revolution, conv(X) is the upper convex envelope of its (z, r)
// profile, computed as a support-function sweep, and the Minkowski sum's
// boundary decomposes over the normal fan of the centers polygon: each
// envelope piece extrudes along every polygon edge and revolves around
// every polygon vertex. Lines become planes and cones, arcs become
// horizontal cylinders and sphere or torus bands; all native surfaces,
// sewn rather than fused (every junction is tangent contact). The axis
// may point anywhere as long as every child shares it. Every result is
// checked against the closed-form Steiner volume integral. Null when the
// children do not fit or the check fails.
namespace OcctRevolvedHull {

TopoDS_Shape hull(const std::vector<TopoDS_Shape>& components);

}  // namespace OcctRevolvedHull
