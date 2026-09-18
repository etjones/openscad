#include "geometry/occt/OcctBoolean.h"

#include <BRepAlgoAPI_BooleanOperation.hxx>
#include <BRepAlgoAPI_Common.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepAlgoAPI_Fuse.hxx>
#include <BRepBuilderAPI_Copy.hxx>
#include <ShapeAnalysis_ShapeTolerance.hxx>
#include <BRepBndLib.hxx>
#include <BRepBuilderAPI_MakeSolid.hxx>
#include <BRepClass3d_SolidClassifier.hxx>
#include <BRepGProp.hxx>
#include <BRep_Builder.hxx>
#include <Bnd_Box.hxx>
#include <GProp_GProps.hxx>
#include <ShapeUpgrade_UnifySameDomain.hxx>
#include <Standard_Failure.hxx>
#include <TopAbs_ShapeEnum.hxx>
#include <TopExp_Explorer.hxx>
#include <TopTools_ListOfShape.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Compound.hxx>
#include <array>
#include <cmath>
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

// Center of mass, if it lies strictly inside the piece (it need not, for
// a non-convex body); used as a cheap interior sample.
bool interiorPoint(const TopoDS_Shape& piece, gp_Pnt& point)
{
  GProp_GProps props;
  BRepGProp::VolumeProperties(piece, props);
  point = props.CentreOfMass();
  BRepClass3d_SolidClassifier classifier(piece, point, 1e-7);
  return classifier.State() == TopAbs_IN;
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
  for (const auto& body : piecesOf(result, dim)) {
    gp_Pnt point;
    if (!interiorPoint(body, point)) continue;
    for (const auto& tool : tools) {
      if (inside(tool, point)) return true;
    }
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

TopoDS_Shape emptyIfNoPieces(const TopoDS_Shape& result, unsigned int dim)
{
  if (result.IsNull()) return {};
  const auto pieces = piecesOf(result, dim);
  if (pieces.empty()) return {};
  if (pieces.size() == 1) return pieces.front();
  return makeCompound(pieces);
}

}  // namespace

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

TopoDS_Shape fuse(const std::vector<TopoDS_Shape>& operands, unsigned int dim)
{
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

  auto failed = [&](const TopoDS_Shape& candidate) {
    if (candidate.IsNull()) return true;
    if (damaged(candidate, dim) || gainedPieces(pieces, candidate, dim)) return true;
    return !bounds.plausible(extent(candidate, dim));
  };

  auto result = run(make, args, tools, 0.0);
  if (!failed(result)) return unify(emptyIfNoPieces(result, dim), dim);

  const double diag = diagonal(pieces);
  for (const double fraction : FUZZ_FRACTIONS) {
    auto candidate = run(make, args, tools, diag * fraction);
    if (!failed(candidate)) return unify(emptyIfNoPieces(candidate, dim), dim);
  }
  OcctBridge::warn(
    "STEP export: a union came back damaged, in more pieces than it was given, or with an "
    "implausible volume; its operands are kept unjoined instead");
  return makeCompound(pieces);
}

TopoDS_Shape cut(const std::vector<TopoDS_Shape>& args, const std::vector<TopoDS_Shape>& tools,
                 unsigned int dim)
{
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

  auto failed = [&](const TopoDS_Shape& candidate) {
    if (candidate.IsNull()) return true;
    if (damaged(candidate, dim)) return true;
    if (!bounds.plausible(extent(candidate, dim))) return true;
    return materialLeftInTools(candidate, toolPieces, dim);
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
  OcctBridge::warn(
    "STEP export: a difference kept material inside the shapes it was cutting with, or returned "
    "an implausible volume, and no retry fixed it; the first result is used");
  return unify(emptyIfNoPieces(result, dim), dim);
}

TopoDS_Shape common(const std::vector<TopoDS_Shape>& args, const std::vector<TopoDS_Shape>& tools,
                    unsigned int dim)
{
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
  OcctBridge::warn(
    "STEP export: an intersection returned an implausible volume and no retry fixed it; the first "
    "result is used");
  return unify(emptyIfNoPieces(result, dim), dim);
}

}  // namespace OcctBoolean
