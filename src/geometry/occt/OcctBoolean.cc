#include "geometry/occt/OcctBoolean.h"

#include <BRepAlgoAPI_BooleanOperation.hxx>
#include <BRepAlgoAPI_Common.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepAlgoAPI_Fuse.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepBndLib.hxx>
#include <BRepBuilderAPI_Copy.hxx>
#include <BRepTools.hxx>
#include <ShapeAnalysis_ShapeTolerance.hxx>
#include <BRepBuilderAPI_MakeSolid.hxx>
#include <BRepClass3d_SolidClassifier.hxx>
#include <BRepGProp.hxx>
#include <BRep_Builder.hxx>
#include <Bnd_Box.hxx>
#include <GProp_GProps.hxx>
#include <ShapeUpgrade_ShapeDivideClosed.hxx>
#include <ShapeUpgrade_UnifySameDomain.hxx>
#include <Standard_Failure.hxx>
#include <TopAbs_ShapeEnum.hxx>
#include <TopExp_Explorer.hxx>
#include <TopTools_ListOfShape.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Compound.hxx>
#include <gp_Vec.hxx>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <functional>
#include <memory>
#include <vector>

#include "geometry/occt/OcctBridge.h"

namespace OcctBoolean {

namespace {

// Relative slack on the input-implied bounds before a result is declared
// implausible; integration noise is orders of magnitude below this.
constexpr double VOLUME_SLACK = 1e-4;

// Relative tolerance for "unify conserved the extent".
constexpr double UNIFY_RTOL = 1e-5;

// A unify may tighten topology, not loosen it: a result whose largest
// tolerance grew past this (absolute, or relative to the input's) is
// rejected. Merging two sphere faces was seen to inflate one vertex to
// half the model's size, after which every boolean downstream was fuzzy at
// model scale and a plane cut kept both halves.
constexpr double UNIFY_TOLERANCE_CAP = 1e-5;
constexpr double UNIFY_TOLERANCE_GROWTH = 100.0;

// How far splitting seams may move a shape's extent. Splitting a periodic
// face re-approximates it near the seam, which moved a scaled sphere's
// volume by 7e-4 relative (towards its exact value, as it happens), so
// this only guards against a split that mangles the shape.
constexpr double SPLIT_RTOL = 5e-3;

// Tolerance for classifying a point against a solid.
constexpr double POINT_TOL = 1e-7;

// Interior points sampled per body when asking whether a boolean kept it.
// One is enough for a body dropped whole; a body half swallowed needs
// several, and eight caught every case seen in the corpus.
constexpr size_t SAMPLES_PER_BODY = 8;

// Fuzzy values tried, as fractions of the operands' bounding diagonal.
constexpr std::array<double, 4> FUZZ_FRACTIONS{1e-7, 1e-6, 1e-5, 1e-4};

TopAbs_ShapeEnum pieceType(unsigned int dim)
{
  return dim == 2 ? TopAbs_FACE : TopAbs_SOLID;
}

double signedExtent(const TopoDS_Shape& shape, unsigned int dim)
{
  if (shape.IsNull()) return 0.0;
  GProp_GProps props;
  if (dim == 2) {
    BRepGProp::SurfaceProperties(shape, props);
  } else {
    BRepGProp::VolumeProperties(shape, props);
  }
  return props.Mass();
}

bool sane(const TopoDS_Shape& piece)
{
  // OCCT signals an outright failure by including a body whose bounding
  // box spans 1e100. Such a body is contagious in later operations.
  Bnd_Box box;
  BRepBndLib::Add(piece, box);
  if (box.IsVoid()) return false;
  double xmin, ymin, zmin, xmax, ymax, zmax;
  box.Get(xmin, ymin, zmin, xmax, ymax, zmax);
  for (const double v : {xmin, ymin, zmin, xmax, ymax, zmax}) {
    if (!std::isfinite(v) || std::fabs(v) > 1e50) return false;
  }
  return true;
}

bool damaged(const TopoDS_Shape& result, unsigned int dim)
{
  for (const auto& piece : piecesOf(result, dim)) {
    if (!sane(piece)) return true;
    if (signedExtent(piece, dim) < 0) return true;
  }
  return false;
}

using Operation = std::function<std::unique_ptr<BRepAlgoAPI_BooleanOperation>()>;

TopoDS_Shape run(const Operation& make, const std::vector<TopoDS_Shape>& args,
                 const std::vector<TopoDS_Shape>& tools, double fuzzy)
{
  try {
    auto op = make();
    TopTools_ListOfShape argList, toolList;
    for (const auto& s : args) argList.Append(s);
    for (const auto& s : tools) toolList.Append(s);
    op->SetArguments(argList);
    op->SetTools(toolList);
    if (fuzzy > 0) op->SetFuzzyValue(fuzzy);
    // Deterministic results: parallel booleans were seen to return
    // different volumes for the same model on successive runs.
    op->SetRunParallel(false);
    op->Build();
    if (!op->IsDone() || op->HasErrors()) return {};
    return op->Shape();
  } catch (const Standard_Failure&) {
    return {};
  }
}

struct Bounds {
  double lo, hi;
  [[nodiscard]] bool plausible(double value) const
  {
    const double slack = VOLUME_SLACK * std::max(hi, 1.0);
    return lo - slack <= value && value <= hi + slack;
  }
};

double sumExtent(const std::vector<TopoDS_Shape>& shapes, unsigned int dim)
{
  double total = 0;
  for (const auto& s : shapes) total += extent(s, dim);
  return total;
}

double maxExtent(const std::vector<TopoDS_Shape>& shapes, unsigned int dim)
{
  double best = 0;
  for (const auto& s : shapes) best = std::max(best, extent(s, dim));
  return best;
}

bool inside(const TopoDS_Shape& solid, const gp_Pnt& point);

// Classifying a point against a solid builds an explorer over the whole
// solid, which costs far more than the classification itself, so the
// classifier is built once per solid and reused for every point. A corpus
// model spent 86% of its export rebuilding these before this existed.
class PointInSolid
{
public:
  explicit PointInSolid(const TopoDS_Shape& solid) : classifier_(solid) {}
  [[nodiscard]] TopAbs_State state(const gp_Pnt& point) const
  {
    classifier_.Perform(point, POINT_TOL);
    return classifier_.State();
  }
  [[nodiscard]] bool contains(const gp_Pnt& point) const { return state(point) == TopAbs_IN; }

private:
  mutable BRepClass3d_SolidClassifier classifier_;
};

// BRepClass3d_SolidClassifier can be neither copied nor moved, so the
// per-solid classifiers are held by pointer.
using Classifiers = std::vector<std::unique_ptr<PointInSolid>>;

Classifiers classifiers(const std::vector<TopoDS_Shape>& solids)
{
  Classifiers out;
  out.reserve(solids.size());
  for (const auto& solid : solids) out.push_back(std::make_unique<PointInSolid>(solid));
  return out;
}

bool anyContains(const Classifiers& solids, const gp_Pnt& point)
{
  return std::any_of(solids.begin(), solids.end(),
                     [&](const auto& solid) { return solid->contains(point); });
}

}  // namespace

double diagonal(const std::vector<TopoDS_Shape>& shapes);

namespace {

// Points strictly inside the piece: its center of mass when that is
// inside (almost always), then points just inside faces spread evenly over
// the piece, for ring-shaped bodies whose center falls in the hole and for
// asking whether a body survived a boolean in more than one place.
std::vector<gp_Pnt> interiorPoints(const TopoDS_Shape& piece, size_t limit)
{
  std::vector<gp_Pnt> found;
  PointInSolid held(piece);
  GProp_GProps props;
  BRepGProp::VolumeProperties(piece, props);
  const gp_Pnt center = props.CentreOfMass();
  if (held.contains(center)) found.push_back(center);
  if (found.size() >= limit) return found;
  std::vector<TopoDS_Face> faces;
  for (TopExp_Explorer it(piece, TopAbs_FACE); it.More(); it.Next())
    faces.push_back(TopoDS::Face(it.Current()));
  if (faces.empty()) return found;
  const double step = diagonal({piece}) * 1e-3;
  const size_t stride = std::max<size_t>(1, faces.size() / (2 * limit));
  size_t tried = 0;
  for (size_t i = 0; i < faces.size() && found.size() < limit && tried < 4 * limit;
       i += stride, ++tried) {
    const auto& face = faces[i];
    BRepAdaptor_Surface surface(face);
    const double u = 0.5 * (surface.FirstUParameter() + surface.LastUParameter());
    const double v = 0.5 * (surface.FirstVParameter() + surface.LastVParameter());
    gp_Pnt middle;
    gp_Vec du, dv;
    surface.D1(u, v, middle, du, dv);
    gp_Vec normal = du.Crossed(dv);
    if (normal.Magnitude() < 1e-12) continue;
    normal.Normalize();
    if (face.Orientation() == TopAbs_REVERSED) normal.Reverse();
    const gp_Pnt candidate = middle.Translated(normal * -step);
    if (held.contains(candidate)) found.push_back(candidate);
  }
  return found;
}

bool interiorPoint(const TopoDS_Shape& piece, gp_Pnt& point)
{
  const auto points = interiorPoints(piece, 1);
  if (points.empty()) return false;
  point = points.front();
  return true;
}

bool inside(const TopoDS_Shape& solid, const gp_Pnt& point)
{
  BRepClass3d_SolidClassifier classifier(solid, point, 1e-7);
  return classifier.State() == TopAbs_IN;
}

// A cut's invariant: no body of the result may sit inside a tool. A
// detector, not a proof (one sample per body), but every observed OCCT cut
// failure violated it.
bool materialLeftInTools(const TopoDS_Shape& result, const std::vector<TopoDS_Shape>& tools,
                         unsigned int dim)
{
  if (dim != 3) return false;
  const auto cutters = classifiers(tools);
  for (const auto& body : piecesOf(result, dim)) {
    gp_Pnt point;
    if (!interiorPoint(body, point)) continue;
    if (anyContains(cutters, point)) return true;
  }
  return false;
}

// A union joins material; it can never break it apart, so it cannot come
// back in more pieces than it was given.
bool gainedPieces(const std::vector<TopoDS_Shape>& operands, const TopoDS_Shape& result,
                  unsigned int dim)
{
  const auto pieces = piecesOf(result, dim).size();
  return pieces > operands.size();
}

void dumpOperands(const char *what, const std::vector<TopoDS_Shape>& args,
                  const std::vector<TopoDS_Shape>& tools)
{
  const char *dir = std::getenv("OPENSCAD_OCCT_DUMP_DIR");
  if (!dir) return;
  static int index = 0;
  const std::string base = std::string(dir) + "/" + what + std::to_string(index++);
  for (size_t i = 0; i < args.size(); ++i) {
    BRepTools::Write(args[i], (base + "-arg" + std::to_string(i) + ".brep").c_str());
  }
  for (size_t i = 0; i < tools.size(); ++i) {
    BRepTools::Write(tools[i], (base + "-tool" + std::to_string(i) + ".brep").c_str());
  }
}

void trace(const std::string& message)
{
  if (std::getenv("OPENSCAD_OCCT_DEBUG")) OcctBridge::warn("occt trace: " + message);
}

TopoDS_Shape emptyIfNoPieces(const TopoDS_Shape& result, unsigned int dim)
{
  if (result.IsNull()) return {};
  const auto pieces = piecesOf(result, dim);
  if (pieces.empty()) return {};
  if (pieces.size() == 1) return pieces.front();
  return makeCompound(pieces);
}

}  // namespace

// A cut keeps everything of its argument that no tool covers: a point
// inside the argument and outside every tool must still be inside the
// result. OCCT's cut was seen to hand back half of a helix, valid-looking
// and within the volume bounds, after a plane cut through its base.
bool keepsUncoveredMaterial(const std::vector<TopoDS_Shape>& args,
                            const std::vector<TopoDS_Shape>& tools, const TopoDS_Shape& result,
                            unsigned int dim)
{
  if (dim != 3) return true;
  const auto bodies = classifiers(piecesOf(result, dim));
  const auto cutters = classifiers(tools);
  for (const auto& arg : args) {
    for (const auto& point : interiorPoints(arg, SAMPLES_PER_BODY)) {
      const bool covered = std::any_of(cutters.begin(), cutters.end(), [&](const auto& tool) {
        return tool->state(point) != TopAbs_OUT;
      });
      if (covered) continue;
      if (!anyContains(bodies, point)) return false;
    }
  }
  return true;
}

// A union contains everything it was given: an operand whose interior
// point is no longer inside the result was dropped. OCCT's fuse was seen
// to lose a plate wedged between a torus and a prism that shared its
// face, with a plausible volume and no extra pieces.
bool keepsOperands(const std::vector<TopoDS_Shape>& operands, const TopoDS_Shape& result,
                   unsigned int dim)
{
  if (dim != 3) return true;
  const auto bodies = classifiers(piecesOf(result, dim));
  for (const auto& operand : operands) {
    for (const auto& point : interiorPoints(operand, SAMPLES_PER_BODY)) {
      if (!anyContains(bodies, point)) return false;
    }
  }
  return true;
}

std::vector<TopoDS_Shape> piecesOf(const TopoDS_Shape& shape, unsigned int dim)
{
  std::vector<TopoDS_Shape> pieces;
  if (shape.IsNull()) return pieces;
  for (TopExp_Explorer it(shape, pieceType(dim)); it.More(); it.Next()) {
    pieces.push_back(it.Current());
  }
  return pieces;
}

double extent(const TopoDS_Shape& shape, unsigned int dim)
{
  double total = 0;
  for (const auto& piece : piecesOf(shape, dim)) {
    total += std::fabs(signedExtent(piece, dim));
  }
  return total;
}

double diagonal(const std::vector<TopoDS_Shape>& shapes)
{
  Bnd_Box box;
  for (const auto& s : shapes) {
    if (!s.IsNull()) BRepBndLib::Add(s, box);
  }
  if (box.IsVoid()) return 1.0;
  double xmin, ymin, zmin, xmax, ymax, zmax;
  box.Get(xmin, ymin, zmin, xmax, ymax, zmax);
  const double d = std::sqrt((xmax - xmin) * (xmax - xmin) + (ymax - ymin) * (ymax - ymin) +
                             (zmax - zmin) * (zmax - zmin));
  return d > 0 ? d : 1.0;
}

TopoDS_Shape makeCompound(const std::vector<TopoDS_Shape>& shapes)
{
  TopoDS_Compound compound;
  BRep_Builder builder;
  builder.MakeCompound(compound);
  for (const auto& s : shapes) {
    if (!s.IsNull()) builder.Add(compound, s);
  }
  return compound;
}

double maxTolerance(const TopoDS_Shape& shape)
{
  ShapeAnalysis_ShapeTolerance analysis;
  return analysis.Tolerance(shape, 1);
}

TopoDS_Shape unify(const TopoDS_Shape& shape, unsigned int dim)
{
  if (shape.IsNull()) return shape;
  try {
    const double before = extent(shape, dim);
    const double toleranceBefore = maxTolerance(shape);
    // UnifySameDomain updates vertex and edge tolerances on the shapes it
    // was given, so a rejected trial would still have damaged the input:
    // it runs on a copy.
    const TopoDS_Shape trial = BRepBuilderAPI_Copy(shape).Shape();
    ShapeUpgrade_UnifySameDomain unifier(trial, true, true, true);
    unifier.Build();
    const auto after = unifier.Shape();
    if (after.IsNull()) return shape;
    const double afterExtent = extent(after, dim);
    if (std::fabs(afterExtent - before) > UNIFY_RTOL * std::max(before, 1e-9) + 1e-9) return shape;
    const double toleranceAfter = maxTolerance(after);
    if (toleranceAfter > std::max(UNIFY_TOLERANCE_CAP, UNIFY_TOLERANCE_GROWTH * toleranceBefore)) {
      return shape;
    }
    return after;
  } catch (const Standard_Failure&) {
    return shape;
  }
}

// A face that wraps all the way around a periodic surface (a whole
// sphere, cylinder or surface of revolution) closes on a seam edge. STEP
// can express such a face, but reading one back is fragile: an ellipsoid
// fused to a cylinder wrote a file whose ellipsoid came back with no
// volume at all (30.74 before the trip, 14.01 after) in OCCT's own reader
// and in the viewer the file was opened in. Splitting every closed face
// in two removes the seam and the file survives. The cost is a couple of
// faces per curved body.
TopoDS_Shape splitClosedFaces(const TopoDS_Shape& shape, unsigned int dim)
{
  if (shape.IsNull() || dim != 3) return shape;
  try {
    const double before = extent(shape, dim);
    ShapeUpgrade_ShapeDivideClosed divider(shape);
    divider.SetNbSplitPoints(1);
    if (!divider.Perform()) return shape;
    const auto after = divider.Result();
    if (after.IsNull()) return shape;
    if (std::fabs(extent(after, dim) - before) > SPLIT_RTOL * std::max(before, 1e-9) + 1e-9) {
      return shape;
    }
    return after;
  } catch (const Standard_Failure&) {
    return shape;
  }
}

TopoDS_Shape fuse(const std::vector<TopoDS_Shape>& operands, unsigned int dim, bool *verified)
{
  if (verified) *verified = true;
  std::vector<TopoDS_Shape> pieces;
  for (const auto& s : operands) {
    for (const auto& p : piecesOf(s, dim)) pieces.push_back(p);
  }
  if (pieces.empty()) return {};
  if (pieces.size() == 1) return pieces.front();

  const Bounds bounds{maxExtent(pieces, dim), sumExtent(pieces, dim)};
  const std::vector<TopoDS_Shape> args{pieces.front()};
  const std::vector<TopoDS_Shape> tools(pieces.begin() + 1, pieces.end());
  const Operation make = [] { return std::make_unique<BRepAlgoAPI_Fuse>(); };
  dumpOperands("fuse", args, tools);

  auto failed = [&](const TopoDS_Shape& candidate, const char *rung) {
    if (candidate.IsNull()) {
      trace(std::string("fuse ") + rung + ": null");
      return true;
    }
    const bool bad = damaged(candidate, dim) || gainedPieces(pieces, candidate, dim) ||
                     !bounds.plausible(extent(candidate, dim)) || !keepsOperands(pieces, candidate, dim);
    trace(std::string("fuse ") + rung + ": extent " + std::to_string(extent(candidate, dim)) +
          " pieces " + std::to_string(piecesOf(candidate, dim).size()) +
          (bad ? " rejected" : " accepted"));
    return bad;
  };

  auto result = run(make, args, tools, 0.0);
  if (!failed(result, "plain")) return unify(emptyIfNoPieces(result, dim), dim);

  const double diag = diagonal(pieces);
  for (const double fraction : FUZZ_FRACTIONS) {
    auto candidate = run(make, args, tools, diag * fraction);
    if (!failed(candidate, ("fuzzy " + std::to_string(fraction)).c_str())) {
      return unify(emptyIfNoPieces(candidate, dim), dim);
    }
  }
  {
    // One operand at a time.
    TopoDS_Shape folded = pieces.front();
    bool ok = true;
    for (size_t i = 1; i < pieces.size(); ++i) {
      auto step = run(make, {folded}, {pieces[i]}, 0.0);
      if (step.IsNull() || damaged(step, dim)) {
        ok = false;
        break;
      }
      folded = emptyIfNoPieces(step, dim);
    }
    if (ok && !failed(folded, "fold")) return unify(folded, dim);
  }
  if (verified) *verified = false;
  OcctBridge::warn(
    "STEP export: a union came back damaged, in more pieces than it was given, missing one of its "
    "operands, or with an implausible volume; its operands are kept unjoined instead");
  return makeCompound(pieces);
}

TopoDS_Shape cut(const std::vector<TopoDS_Shape>& args, const std::vector<TopoDS_Shape>& tools,
                 unsigned int dim, bool *verified)
{
  if (verified) *verified = true;
  std::vector<TopoDS_Shape> argPieces, toolPieces;
  for (const auto& s : args) {
    for (const auto& p : piecesOf(s, dim)) argPieces.push_back(p);
  }
  for (const auto& s : tools) {
    for (const auto& p : piecesOf(s, dim)) toolPieces.push_back(p);
  }
  if (argPieces.empty()) return {};
  if (toolPieces.empty()) return argPieces.size() == 1 ? argPieces.front() : makeCompound(argPieces);

  const double argExtent = sumExtent(argPieces, dim);
  const Bounds bounds{std::max(argExtent - sumExtent(toolPieces, dim), 0.0), argExtent};
  const Operation make = [] { return std::make_unique<BRepAlgoAPI_Cut>(); };
  dumpOperands("cut", argPieces, toolPieces);

  auto failed = [&](const TopoDS_Shape& candidate) {
    if (candidate.IsNull()) return true;
    if (damaged(candidate, dim)) return true;
    if (!bounds.plausible(extent(candidate, dim))) return true;
    return materialLeftInTools(candidate, toolPieces, dim) ||
           !keepsUncoveredMaterial(argPieces, toolPieces, candidate, dim);
  };

  auto result = run(make, argPieces, toolPieces, 0.0);
  if (!failed(result)) return unify(emptyIfNoPieces(result, dim), dim);

  const double diag = diagonal(argPieces);
  for (const double fraction : FUZZ_FRACTIONS) {
    auto candidate = run(make, argPieces, toolPieces, diag * fraction);
    if (!failed(candidate)) return unify(emptyIfNoPieces(candidate, dim), dim);
  }
  if (toolPieces.size() > 1) {
    // One tool at a time.
    std::vector<TopoDS_Shape> current = argPieces;
    bool ok = true;
    for (const auto& tool : toolPieces) {
      auto step = run(make, current, {tool}, 0.0);
      if (step.IsNull() || damaged(step, dim)) {
        ok = false;
        break;
      }
      current = piecesOf(step, dim);
      if (current.empty()) return {};
    }
    if (ok) {
      auto folded = current.size() == 1 ? current.front() : makeCompound(current);
      if (!failed(folded)) return unify(folded, dim);
    }
  }
  if (verified) *verified = false;
  OcctBridge::warn(
    "STEP export: a difference kept material inside the shapes it was cutting with, came back "
    "with an open shell, or returned an implausible volume, and no retry fixed it");
  return unify(emptyIfNoPieces(result, dim), dim);
}

TopoDS_Shape common(const std::vector<TopoDS_Shape>& args, const std::vector<TopoDS_Shape>& tools,
                    unsigned int dim, bool *verified)
{
  if (verified) *verified = true;
  std::vector<TopoDS_Shape> argPieces, toolPieces;
  for (const auto& s : args) {
    for (const auto& p : piecesOf(s, dim)) argPieces.push_back(p);
  }
  for (const auto& s : tools) {
    for (const auto& p : piecesOf(s, dim)) toolPieces.push_back(p);
  }
  if (argPieces.empty() || toolPieces.empty()) return {};

  const Bounds bounds{0.0, std::min(sumExtent(argPieces, dim), sumExtent(toolPieces, dim))};
  const Operation make = [] { return std::make_unique<BRepAlgoAPI_Common>(); };

  auto failed = [&](const TopoDS_Shape& candidate) {
    if (candidate.IsNull()) return true;
    if (damaged(candidate, dim)) return true;
    return !bounds.plausible(extent(candidate, dim));
  };

  auto result = run(make, argPieces, toolPieces, 0.0);
  if (!failed(result)) return unify(emptyIfNoPieces(result, dim), dim);

  const double diag = diagonal(argPieces);
  for (const double fraction : FUZZ_FRACTIONS) {
    auto candidate = run(make, argPieces, toolPieces, diag * fraction);
    if (!failed(candidate)) return unify(emptyIfNoPieces(candidate, dim), dim);
  }
  if (verified) *verified = false;
  OcctBridge::warn("STEP export: an intersection returned an implausible volume and no retry fixed it");
  return unify(emptyIfNoPieces(result, dim), dim);
}

}  // namespace OcctBoolean
