#include "geometry/occt/OcctRevolvedHull.h"

#include <BRepAdaptor_Surface.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakePolygon.hxx>
#include <BRepBuilderAPI_MakeSolid.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepBuilderAPI_Sewing.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepOffsetAPI_ThruSections.hxx>
#include <BRepPrimAPI_MakePrism.hxx>
#include <BRepPrimAPI_MakeRevol.hxx>
#include <BRep_Tool.hxx>
#include <GC_MakeArcOfCircle.hxx>
#include <GeomAbs_SurfaceType.hxx>
#include <Geom_TrimmedCurve.hxx>
#include <ShapeFix_Solid.hxx>
#include <Standard_Failure.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Shell.hxx>
#include <TopoDS_Vertex.hxx>
#include <TopoDS_Wire.hxx>
#include <gp_Ax1.hxx>
#include <gp_Circ.hxx>
#include <gp_Pnt.hxx>
#include <gp_Trsf.hxx>
#include <gp_Vec.hxx>
#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <optional>
#include <set>
#include <utility>
#include <vector>

#include "geometry/occt/OcctBoolean.h"

namespace OcctRevolvedHull {

namespace {

constexpr double PARALLEL_TOL = 1e-6;
constexpr double ANG_EPS = 1e-12;
constexpr double TWO_PI = 2 * M_PI;

using P2 = std::array<double, 2>;
using P3 = std::array<double, 3>;
using Span = std::pair<double, double>;
using Trims = std::map<int, std::vector<Span>>;

double dist2(const P2& a, const P2& b)
{
  return std::hypot(a[0] - b[0], a[1] - b[1]);
}

// An element of the (z, r) profile: a point, or a circle (center, rho)
// whose trimmed angular spans are recorded under `key`.
struct Elem {
  bool isPoint;
  P2 p;  // the point, or the circle's center (z, r)
  double rho = 0;
  int key = -1;
};

// A piece of the envelope chain: a line from a to b, or an arc of the
// circle (c, rho) with the direction angle t running from hi down to lo.
struct Piece {
  bool isLine;
  P2 a, b;
  P2 c;
  double rho = 0, hi = 0, lo = 0;
  int key = -1;
};

P2 arcPoint(const Piece& piece, double t)
{
  return {piece.c[0] + piece.rho * std::cos(t), piece.c[1] + piece.rho * std::sin(t)};
}

P2 pieceStart(const Piece& piece)
{
  return piece.isLine ? piece.a : arcPoint(piece, piece.hi);
}
P2 pieceEnd(const Piece& piece)
{
  return piece.isLine ? piece.b : arcPoint(piece, piece.lo);
}

double mod2pi(double a)
{
  a = std::fmod(a, TWO_PI);
  return a < 0 ? a + TWO_PI : a;
}

bool psiInSpans(double psi, const std::vector<Span>& spans)
{
  for (const auto& [a, b] : spans) {
    if (b - a >= TWO_PI - 1e-6) return true;
    if (mod2pi(psi - a) <= mod2pi(b - a) + 1e-9) return true;
  }
  return false;
}

double elemSupport(const Elem& e, double psi, const Trims& trims)
{
  const double ux = std::cos(psi), uy = std::sin(psi);
  if (e.isPoint) return e.p[0] * ux + e.p[1] * uy;
  auto it = trims.find(e.key);
  if (it != trims.end() && !psiInSpans(psi, it->second)) return -1e300;
  return e.p[0] * ux + e.p[1] * uy + e.rho;
}

// psi in (0, pi) where two elements' support values tie.
std::vector<double> switchAngles(const Elem& e1, const Elem& e2)
{
  const double r1 = e1.isPoint ? 0.0 : e1.rho, r2 = e2.isPoint ? 0.0 : e2.rho;
  const double a = e1.p[0] - e2.p[0], b = e1.p[1] - e2.p[1], k = r2 - r1;
  const double m = std::hypot(a, b);
  std::vector<double> out;
  if (m < 1e-12 || std::fabs(k) > m) return out;
  const double base = std::atan2(b, a);
  const double off = std::acos(std::clamp(k / m, -1.0, 1.0));
  for (const double cand : {mod2pi(base + off), mod2pi(base - off)}) {
    if (cand > ANG_EPS && cand < M_PI - ANG_EPS) out.push_back(cand);
  }
  return out;
}

// Upper convex envelope of points and circles in the (z, r) plane, as an
// ordered chain from min z to max z. Between consecutive switch angles
// the supporting element is constant and contributes its boundary.
std::vector<Piece> upperEnvelope(const std::vector<Elem>& elements, const Trims& trims)
{
  std::set<double> angles{ANG_EPS, M_PI - ANG_EPS};
  for (size_t i = 0; i < elements.size(); ++i) {
    for (size_t j = i + 1; j < elements.size(); ++j) {
      for (const double a : switchAngles(elements[i], elements[j])) angles.insert(a);
    }
  }
  for (const auto& [key, spans] : trims) {
    for (const auto& span : spans) {
      for (const double raw : {span.first, span.second}) {
        const double a = mod2pi(raw);
        if (a > ANG_EPS && a < M_PI - ANG_EPS) angles.insert(a);
      }
    }
  }
  std::vector<double> ordered;
  for (auto it = angles.rbegin(); it != angles.rend(); ++it) {
    if (ordered.empty() || ordered.back() - *it > 1e-9) ordered.push_back(*it);
  }

  struct Interval {
    size_t best;
    double hi, lo;
  };
  std::vector<Interval> intervals;
  for (size_t i = 0; i + 1 < ordered.size(); ++i) {
    const double hi = ordered[i], lo = ordered[i + 1], mid = (hi + lo) / 2;
    size_t best = 0;
    double bestValue = -1e301;
    for (size_t e = 0; e < elements.size(); ++e) {
      const double v = elemSupport(elements[e], mid, trims);
      if (v > bestValue) {
        bestValue = v;
        best = e;
      }
    }
    if (!intervals.empty() && intervals.back().best == best) {
      intervals.back().lo = lo;
    } else {
      intervals.push_back({best, hi, lo});
    }
  }

  std::vector<Piece> chain;
  std::optional<P2> prevEnd;
  for (const auto& iv : intervals) {
    const Elem& e = elements[iv.best];
    if (e.isPoint) {
      if (prevEnd && dist2(*prevEnd, e.p) > 1e-6) {
        chain.push_back({true, *prevEnd, e.p, {}, 0, 0, 0, -1});
      }
      prevEnd = e.p;
    } else {
      const P2 start{e.p[0] + e.rho * std::cos(iv.hi), e.p[1] + e.rho * std::sin(iv.hi)};
      const P2 end{e.p[0] + e.rho * std::cos(iv.lo), e.p[1] + e.rho * std::sin(iv.lo)};
      if (prevEnd && dist2(*prevEnd, start) > 1e-6) {
        chain.push_back({true, *prevEnd, start, {}, 0, 0, 0, -1});
      }
      if (!chain.empty() && !chain.back().isLine && chain.back().key == e.key &&
          std::fabs(chain.back().lo - iv.hi) < 1e-9) {
        chain.back().lo = iv.lo;
      } else {
        chain.push_back({false, {}, {}, e.p, e.rho, iv.hi, iv.lo, e.key});
      }
      prevEnd = end;
    }
  }
  return chain;
}

// Is [lo, hi] covered by the union of (possibly wrapping) spans?
bool spanCovered(double lo, double hi, const std::vector<Span>& spans)
{
  std::vector<Span> segs;
  for (auto [a, b] : spans) {
    if (b - a >= TWO_PI - 1e-6) return true;
    a = mod2pi(a);
    b = mod2pi(b);
    if (b < a - 1e-12) {
      segs.emplace_back(a, TWO_PI);
      segs.emplace_back(0.0, b);
    } else {
      segs.emplace_back(a, b);
    }
  }
  std::vector<Span> need;
  if (hi - lo > 2e-6) need.emplace_back(lo + 1e-6, hi - 1e-6);
  std::sort(segs.begin(), segs.end());
  for (const auto& [a, b] : segs) {
    std::vector<Span> next;
    for (const auto& [nlo, nhi] : need) {
      for (const Span piece : {Span{nlo, std::min(nhi, a)}, Span{std::max(nlo, b), nhi}}) {
        if (piece.second - piece.first > 1e-9) next.push_back(piece);
      }
    }
    need = std::move(next);
  }
  return need.empty();
}

double arcFR(double rc, double rho, double t)
{
  return -rc * rho * std::cos(t) + rho * rho * (t - std::sin(t) * std::cos(t)) / 2;
}

double arcFR2(double rc, double rho, double t)
{
  const double c = std::cos(t);
  return -rc * rc * rho * c + rc * rho * rho * (t - std::sin(t) * c) +
         rho * rho * rho * (-c + c * c * c / 3);
}

// (dz, integral r dz, integral r^2 dz) over the chain: everything the
// Steiner cross-section volume needs, exactly.
std::array<double, 3> chainIntegrals(const std::vector<Piece>& chain)
{
  double dzTotal = 0, intR = 0, intR2 = 0;
  for (const auto& piece : chain) {
    if (piece.isLine) {
      const double z1 = piece.a[0], r1 = piece.a[1], z2 = piece.b[0], r2 = piece.b[1];
      const double dz = z2 - z1;
      dzTotal += dz;
      intR += dz * (r1 + r2) / 2;
      intR2 += dz * (r1 * r1 + r1 * r2 + r2 * r2) / 3;
    } else {
      const double rc = piece.c[1], rho = piece.rho;
      dzTotal += rho * (std::cos(piece.lo) - std::cos(piece.hi));
      intR += arcFR(rc, rho, piece.hi) - arcFR(rc, rho, piece.lo);
      intR2 += arcFR2(rc, rho, piece.hi) - arcFR2(rc, rho, piece.lo);
    }
  }
  return {dzTotal, intR, intR2};
}

bool chainsMatch(const std::vector<Piece>& a, const std::vector<Piece>& b)
{
  if (a.size() != b.size()) return false;
  auto close = [](double x, double y) { return std::fabs(x - y) <= 1e-7 + 1e-6 * std::fabs(y); };
  for (size_t i = 0; i < a.size(); ++i) {
    if (a[i].isLine != b[i].isLine) return false;
    if (a[i].isLine) {
      if (!close(a[i].a[0], b[i].a[0]) || !close(a[i].a[1], b[i].a[1]) || !close(a[i].b[0], b[i].b[0]) ||
          !close(a[i].b[1], b[i].b[1])) {
        return false;
      }
    } else if (!close(a[i].c[0], b[i].c[0]) || !close(a[i].c[1], b[i].c[1]) ||
               !close(a[i].rho, b[i].rho) || !close(a[i].hi, b[i].hi) || !close(a[i].lo, b[i].lo)) {
      return false;
    }
  }
  return true;
}

// --- geometry helpers ------------------------------------------------------

TopoDS_Shape faceFromWire(const TopoDS_Wire& wire)
{
  BRepBuilderAPI_MakeFace maker(wire, true);
  return maker.IsDone() ? TopoDS_Shape(maker.Face()) : TopoDS_Shape();
}

TopoDS_Edge lineEdge(const gp_Pnt& a, const gp_Pnt& b)
{
  return BRepBuilderAPI_MakeEdge(a, b).Edge();
}

TopoDS_Edge tangentArc(const gp_Pnt& start, const gp_Vec& tangent, const gp_Pnt& end)
{
  return BRepBuilderAPI_MakeEdge(GC_MakeArcOfCircle(start, tangent, end).Value()).Edge();
}

TopoDS_Edge threePointArc(const gp_Pnt& a, const gp_Pnt& m, const gp_Pnt& b)
{
  return BRepBuilderAPI_MakeEdge(GC_MakeArcOfCircle(a, m, b).Value()).Edge();
}

TopoDS_Shape circleFace(const P2& c, double r, double z)
{
  const gp_Circ circ(gp_Ax2(gp_Pnt(c[0], c[1], z), gp_Dir(0, 0, 1)), r);
  return faceFromWire(BRepBuilderAPI_MakeWire(BRepBuilderAPI_MakeEdge(circ).Edge()).Wire());
}

TopoDS_Shape stadiumFace(const P2& a, const P2& b, double r, double z)
{
  const double length = dist2(a, b);
  const double ux = (b[0] - a[0]) / length, uy = (b[1] - a[1]) / length;
  const double nx = -uy, ny = ux;
  const gp_Pnt c1(a[0] + nx * r, a[1] + ny * r, z), c2(b[0] + nx * r, b[1] + ny * r, z);
  const gp_Pnt c3(b[0] - nx * r, b[1] - ny * r, z), c4(a[0] - nx * r, a[1] - ny * r, z);
  BRepBuilderAPI_MakeWire wire;
  wire.Add(lineEdge(c1, c2));
  wire.Add(tangentArc(c2, gp_Vec(ux, uy, 0), c3));
  wire.Add(lineEdge(c3, c4));
  wire.Add(tangentArc(c4, gp_Vec(-ux, -uy, 0), c1));
  return wire.IsDone() ? faceFromWire(wire.Wire()) : TopoDS_Shape();
}

// The rounded polygon (2D hull of the centers offset by r) as a face at
// height z, as one deterministic wire so every section of a loft has
// identical edge structure and orientation.
TopoDS_Shape roundedSection(const std::vector<P2>& hullPts, double r, double z)
{
  const size_t k = hullPts.size();
  if (k == 1) return circleFace(hullPts[0], r, z);
  if (k == 2) return stadiumFace(hullPts[0], hullPts[1], r, z);
  std::vector<P2> normals;
  for (size_t i = 0; i < k; ++i) {
    const auto& a = hullPts[i];
    const auto& b = hullPts[(i + 1) % k];
    const double d = dist2(a, b);
    normals.push_back({(b[1] - a[1]) / d, -(b[0] - a[0]) / d});
  }
  BRepBuilderAPI_MakeWire wire;
  for (size_t i = 0; i < k; ++i) {
    const auto& a = hullPts[i];
    const auto& b = hullPts[(i + 1) % k];
    const auto& n = normals[i];
    const auto& n2 = normals[(i + 1) % k];
    wire.Add(lineEdge(gp_Pnt(a[0] + r * n[0], a[1] + r * n[1], z),
                      gp_Pnt(b[0] + r * n[0], b[1] + r * n[1], z)));
    const double d = dist2(a, b);
    wire.Add(tangentArc(gp_Pnt(b[0] + r * n[0], b[1] + r * n[1], z),
                        gp_Vec((b[0] - a[0]) / d, (b[1] - a[1]) / d, 0),
                        gp_Pnt(b[0] + r * n2[0], b[1] + r * n2[1], z)));
  }
  return wire.IsDone() ? faceFromWire(wire.Wire()) : TopoDS_Shape();
}

TopoDS_Wire outerWire(const TopoDS_Shape& face)
{
  TopoDS_Wire wire;
  for (TopExp_Explorer it(face, TopAbs_WIRE); it.More(); it.Next()) wire = TopoDS::Wire(it.Current());
  return wire;
}

struct Frame {
  std::vector<std::pair<size_t, size_t>> edges;
  std::vector<double> normals;  // outward-normal angle of each directed edge
};

Frame polygonFrame(const std::vector<P2>& hullPts)
{
  Frame frame;
  const size_t k = hullPts.size();
  if (k >= 3) {
    for (size_t i = 0; i < k; ++i) frame.edges.emplace_back(i, (i + 1) % k);
  } else if (k == 2) {
    frame.edges = {{0, 1}, {1, 0}};
  }
  for (const auto& [i, j] : frame.edges) {
    const double dx = hullPts[j][0] - hullPts[i][0], dy = hullPts[j][1] - hullPts[i][1];
    frame.normals.push_back(std::atan2(-dx, dy));
  }
  return frame;
}

TopoDS_Edge curveEdge(const Piece& piece, double cx, double cy, double theta)
{
  const double ct = std::cos(theta), st = std::sin(theta);
  auto p3 = [&](const P2& zr) { return gp_Pnt(cx + zr[1] * ct, cy + zr[1] * st, zr[0]); };
  if (piece.isLine) return lineEdge(p3(piece.a), p3(piece.b));
  return threePointArc(p3(arcPoint(piece, piece.hi)), p3(arcPoint(piece, (piece.hi + piece.lo) / 2)),
                       p3(arcPoint(piece, piece.lo)));
}

TopoDS_Shape solidFromFaces(const std::vector<TopoDS_Shape>& faces, double tolerance)
{
  BRepBuilderAPI_Sewing sewing(tolerance);
  for (const auto& f : faces) sewing.Add(f);
  sewing.Perform();
  TopoDS_Shell shell;
  for (TopExp_Explorer it(sewing.SewedShape(), TopAbs_SHELL); it.More(); it.Next()) {
    shell = TopoDS::Shell(it.Current());
    break;
  }
  if (shell.IsNull()) return {};
  BRepBuilderAPI_MakeSolid maker(shell);
  if (!maker.IsDone()) return {};
  ShapeFix_Solid fix(maker.Solid());
  fix.Perform();
  return fix.Solid();
}

// Realize conv(centers) (+) conv(profile) face by face over the normal
// fan: every envelope piece extrudes along every polygon edge and
// revolves through every polygon vertex wedge, plus flat caps.
TopoDS_Shape fanSolid(const std::vector<P2>& hullPts, const std::vector<Piece>& chain)
{
  const size_t k = hullPts.size();
  const auto frame = polygonFrame(hullPts);
  std::vector<TopoDS_Shape> faces;
  for (const auto& piece : chain) {
    for (size_t e = 0; e < frame.edges.size(); ++e) {
      const auto [i, j] = frame.edges[e];
      const auto edge = curveEdge(piece, hullPts[i][0], hullPts[i][1], frame.normals[e]);
      const gp_Vec vec(hullPts[j][0] - hullPts[i][0], hullPts[j][1] - hullPts[i][1], 0);
      faces.push_back(BRepPrimAPI_MakePrism(edge, vec).Shape());
    }
    for (size_t i = 0; i < k; ++i) {
      double aIn = 0, sweep = TWO_PI;
      if (k >= 2) {
        const size_t n = frame.edges.size();
        aIn = frame.normals[(i + n - 1) % n];
        const double aOut = k >= 3 ? frame.normals[i % n] : aIn + M_PI;
        sweep = mod2pi(aOut - aIn);
      }
      const auto edge = curveEdge(piece, hullPts[i][0], hullPts[i][1], aIn);
      const gp_Ax1 axis(gp_Pnt(hullPts[i][0], hullPts[i][1], 0), gp_Dir(0, 0, 1));
      faces.push_back(BRepPrimAPI_MakeRevol(edge, axis, sweep).Shape());
    }
  }
  const P2 first = pieceStart(chain.front()), last = pieceEnd(chain.back());
  for (const auto& [z, r] : {first, last}) {
    if (r > 1e-7) {
      faces.push_back(roundedSection(hullPts, r, z));
    } else if (k >= 3) {
      BRepBuilderAPI_MakePolygon polygon;
      for (const auto& p : hullPts) polygon.Add(gp_Pnt(p[0], p[1], z));
      polygon.Close();
      faces.push_back(faceFromWire(polygon.Wire()));
    }
  }
  for (const auto& f : faces) {
    if (f.IsNull()) return {};
  }
  return solidFromFaces(faces, 1e-5);
}

// All-line profiles: a ruled loft of rounded-polygon sections at the
// envelope breakpoints.
TopoDS_Shape loftSolid(const std::vector<P2>& hullPts, const std::vector<Piece>& chain)
{
  std::vector<P2> breakpoints{pieceStart(chain.front())};
  for (const auto& piece : chain) breakpoints.push_back(pieceEnd(piece));
  for (const auto& [z, r] : breakpoints) {
    if (r < 1e-9) return {};
  }
  BRepOffsetAPI_ThruSections loft(true, true);
  for (const auto& [z, r] : breakpoints) {
    auto section = roundedSection(hullPts, r, z);
    if (section.IsNull()) return {};
    loft.AddWire(outerWire(section));
  }
  loft.Build();
  return loft.IsDone() ? loft.Shape() : TopoDS_Shape();
}

// --- classification --------------------------------------------------------

struct Profile {
  P2 center;  // axis position in the XY plane
  std::vector<Elem> elements;
  Trims trims;
};

// Circle keys: (zc, rc, rho) matched to 1e-6, shared across children so
// a torus band split by a boolean still lands on one element.
struct KeyRegistry {
  std::vector<std::array<double, 3>> keys;
  int lookup(double zc, double rc, double rho)
  {
    for (size_t i = 0; i < keys.size(); ++i) {
      if (std::fabs(keys[i][0] - zc) < 1e-6 && std::fabs(keys[i][1] - rc) < 1e-6 &&
          std::fabs(keys[i][2] - rho) < 1e-6) {
        return static_cast<int>(i);
      }
    }
    keys.push_back({zc, rc, rho});
    return static_cast<int>(keys.size() - 1);
  }
};

bool alongZ(const gp_Dir& d)
{
  return std::fabs(std::fabs(d.Z()) - 1.0) <= PARALLEL_TOL;
}

std::optional<Profile> profileElements(const TopoDS_Shape& shape, KeyRegistry& registry)
{
  bool hasSolid = false;
  for (TopExp_Explorer it(shape, TopAbs_SOLID); it.More(); it.Next()) hasSolid = true;
  if (!hasSolid) return std::nullopt;

  struct Circle {
    double zc, rc, rho;
    Span span;
  };
  std::optional<P2> axisXY;
  std::vector<Circle> circles;
  auto sameAxis = [&](const P2& xy) {
    if (!axisXY) {
      axisXY = xy;
      return true;
    }
    return dist2(xy, *axisXY) <= 1e-6;
  };
  for (TopExp_Explorer it(shape, TopAbs_FACE); it.More(); it.Next()) {
    BRepAdaptor_Surface surface(TopoDS::Face(it.Current()));
    const auto kind = surface.GetType();
    if (kind == GeomAbs_Plane) {
      if (!alongZ(surface.Plane().Axis().Direction())) return std::nullopt;
      continue;
    }
    P2 xy;
    if (kind == GeomAbs_Sphere) {
      const auto sph = surface.Sphere();
      const auto loc = sph.Location();
      xy = {loc.X(), loc.Y()};
      const double s = sph.Position().Direction().Z();
      if (std::fabs(std::fabs(s) - 1.0) > PARALLEL_TOL) return std::nullopt;
      const double v1 = surface.FirstVParameter(), v2 = surface.LastVParameter();
      const Span span = s > 0 ? Span{M_PI / 2 - v2, M_PI / 2 - v1} : Span{M_PI / 2 + v1, M_PI / 2 + v2};
      circles.push_back({loc.Z(), 0.0, sph.Radius(), span});
    } else if (kind == GeomAbs_Torus) {
      const auto tor = surface.Torus();
      const auto ax = tor.Axis();
      if (!alongZ(ax.Direction())) return std::nullopt;
      if (tor.MinorRadius() > tor.MajorRadius() - 1e-9) return std::nullopt;
      const auto loc = ax.Location();
      xy = {loc.X(), loc.Y()};
      const double s = ax.Direction().Z();
      const double v1 = surface.FirstVParameter(), v2 = surface.LastVParameter();
      const Span span = s > 0 ? Span{M_PI / 2 - v2, M_PI / 2 - v1} : Span{M_PI / 2 + v1, M_PI / 2 + v2};
      circles.push_back({loc.Z(), tor.MajorRadius(), tor.MinorRadius(), span});
    } else if (kind == GeomAbs_Cylinder) {
      const auto ax = surface.Cylinder().Axis();
      if (!alongZ(ax.Direction())) return std::nullopt;
      xy = {ax.Location().X(), ax.Location().Y()};
    } else if (kind == GeomAbs_Cone) {
      const auto ax = surface.Cone().Axis();
      if (!alongZ(ax.Direction())) return std::nullopt;
      xy = {ax.Location().X(), ax.Location().Y()};
    } else {
      return std::nullopt;
    }
    if (!sameAxis(xy)) return std::nullopt;
  }
  if (!axisXY) return std::nullopt;

  Profile profile;
  profile.center = *axisXY;
  for (TopExp_Explorer it(shape, TopAbs_VERTEX); it.More(); it.Next()) {
    const auto p = BRep_Tool::Pnt(TopoDS::Vertex(it.Current()));
    profile.elements.push_back(
      {true, {p.Z(), std::hypot(p.X() - (*axisXY)[0], p.Y() - (*axisXY)[1])}, 0, -1});
  }
  std::set<int> seen;
  for (const auto& c : circles) {
    const int key = registry.lookup(c.zc, c.rc, c.rho);
    if (seen.insert(key).second) profile.elements.push_back({false, {c.zc, c.rc}, c.rho, key});
    profile.trims[key].push_back(c.span);
    for (const double a : {c.span.first, c.span.second}) {
      profile.elements.push_back(
        {true, {c.zc + c.rho * std::cos(a), c.rc + c.rho * std::sin(a)}, 0, -1});
    }
  }
  return profile;
}

std::vector<P2> centersHull2d(const std::vector<P2>& centers)
{
  // Extremes, collinearity, then Andrew's monotone chain.
  double span = 0;
  P2 a = centers[0], b = centers[0];
  for (size_t i = 0; i < centers.size(); ++i) {
    for (size_t j = i + 1; j < centers.size(); ++j) {
      const double d = dist2(centers[i], centers[j]);
      if (d > span) {
        span = d;
        a = centers[i];
        b = centers[j];
      }
    }
  }
  if (span < 1e-9) return {centers[0]};
  bool collinear = true;
  const double ux = (b[0] - a[0]) / span, uy = (b[1] - a[1]) / span;
  for (const auto& p : centers) {
    const double t = (p[0] - a[0]) * ux + (p[1] - a[1]) * uy;
    if (dist2(p, {a[0] + t * ux, a[1] + t * uy}) > 1e-6 * span) collinear = false;
  }
  if (collinear) return {a, b};
  std::vector<P2> pts = centers;
  std::sort(pts.begin(), pts.end());
  pts.erase(
    std::unique(pts.begin(), pts.end(), [](const P2& x, const P2& y) { return dist2(x, y) < 1e-9; }),
    pts.end());
  auto cross = [](const P2& o, const P2& p, const P2& q) {
    return (p[0] - o[0]) * (q[1] - o[1]) - (p[1] - o[1]) * (q[0] - o[0]);
  };
  std::vector<P2> hull(2 * pts.size());
  size_t k = 0;
  for (const auto& p : pts) {
    while (k >= 2 && cross(hull[k - 2], hull[k - 1], p) <= 1e-12) --k;
    hull[k++] = p;
  }
  for (size_t i = pts.size() - 1, t = k + 1; i-- > 0;) {
    while (k >= t && cross(hull[k - 2], hull[k - 1], pts[i]) <= 1e-12) --k;
    hull[k++] = pts[i];
  }
  hull.resize(k - 1);
  return hull;
}

TopoDS_Shape hullVerticalTranslates(const std::vector<TopoDS_Shape>& shapes)
{
  KeyRegistry registry;
  std::vector<Profile> groups;
  for (const auto& s : shapes) {
    auto profile = profileElements(s, registry);
    if (!profile) return {};
    bool merged = false;
    for (auto& g : groups) {
      if (dist2(g.center, profile->center) <= 1e-6) {
        g.elements.insert(g.elements.end(), profile->elements.begin(), profile->elements.end());
        for (const auto& [key, spans] : profile->trims) {
          g.trims[key].insert(g.trims[key].end(), spans.begin(), spans.end());
        }
        merged = true;
        break;
      }
    }
    if (!merged) groups.push_back(std::move(*profile));
  }

  std::vector<std::vector<Piece>> chains;
  for (const auto& g : groups) {
    auto chain = upperEnvelope(g.elements, g.trims);
    for (const auto& piece : chain) {
      if (piece.isLine) continue;
      auto it = g.trims.find(piece.key);
      const std::vector<Span> spans = it == g.trims.end() ? std::vector<Span>{} : it->second;
      if (!spanCovered(piece.lo, piece.hi, spans)) return {};
    }
    chains.push_back(std::move(chain));
  }
  const auto& chain = chains.front();
  for (size_t i = 1; i < chains.size(); ++i) {
    if (!chainsMatch(chain, chains[i])) return {};
  }
  if (chain.empty()) return {};
  const P2 first = pieceStart(chain.front()), last = pieceEnd(chain.back());
  if (last[0] - first[0] < 1e-9) return {};

  std::vector<P2> centers;
  for (const auto& g : groups) centers.push_back(g.center);
  const auto hullPts = centersHull2d(centers);
  bool allLines = true;
  for (const auto& piece : chain) allLines = allLines && piece.isLine;
  TopoDS_Shape solid;
  try {
    solid = allLines ? loftSolid(hullPts, chain) : fanSolid(hullPts, chain);
  } catch (const Standard_Failure&) {
    return {};
  }
  if (solid.IsNull() || !BRepCheck_Analyzer(solid).IsValid()) return {};

  double area = 0, perimeter = 0;
  const size_t k = hullPts.size();
  if (k < 3) {
    perimeter = 2 * dist2(hullPts.front(), hullPts.back());
  } else {
    for (size_t i = 0; i < k; ++i) {
      const auto& p = hullPts[i];
      const auto& q = hullPts[(i + 1) % k];
      area += p[0] * q[1] - q[0] * p[1];
      perimeter += dist2(p, q);
    }
    area = std::fabs(area) / 2;
  }
  const auto [dz, intR, intR2] = chainIntegrals(chain);
  const double exact = area * dz + perimeter * intR + M_PI * intR2;
  const double got = OcctBoolean::extent(solid, 3);
  if (std::fabs(got - exact) > 1e-6 * std::fabs(exact)) return {};
  return solid;
}

// The unit axis shared by a child's cylinder/cone/torus faces. `free` when
// nothing constrains it (spheres and planes only); nullopt when faces
// disagree or an unsupported surface appears.
std::optional<P3> childAxisDirection(const TopoDS_Shape& shape, bool& free, bool& ok)
{
  ok = true;
  free = true;
  bool hasSolid = false;
  for (TopExp_Explorer it(shape, TopAbs_SOLID); it.More(); it.Next()) hasSolid = true;
  if (!hasSolid) {
    ok = false;
    return std::nullopt;
  }
  std::optional<P3> direction;
  for (TopExp_Explorer it(shape, TopAbs_FACE); it.More(); it.Next()) {
    BRepAdaptor_Surface surface(TopoDS::Face(it.Current()));
    const auto kind = surface.GetType();
    if (kind == GeomAbs_Plane || kind == GeomAbs_Sphere) continue;
    gp_Dir d;
    if (kind == GeomAbs_Cylinder) d = surface.Cylinder().Axis().Direction();
    else if (kind == GeomAbs_Cone) d = surface.Cone().Axis().Direction();
    else if (kind == GeomAbs_Torus) d = surface.Torus().Axis().Direction();
    else {
      ok = false;
      return std::nullopt;
    }
    const P3 v{d.X(), d.Y(), d.Z()};
    if (!direction) {
      direction = v;
      free = false;
    } else if (std::fabs(
                 std::fabs((*direction)[0] * v[0] + (*direction)[1] * v[1] + (*direction)[2] * v[2]) -
                 1.0) > 1e-9) {
      ok = false;
      return std::nullopt;
    }
  }
  return direction;
}

}  // namespace

TopoDS_Shape hull(const std::vector<TopoDS_Shape>& components)
{
  if (components.empty()) return {};
  std::optional<P3> d0;
  for (const auto& s : components) {
    bool free = true, ok = true;
    auto d = childAxisDirection(s, free, ok);
    if (!ok) return {};
    if (free) continue;
    if (!d0) {
      d0 = d;
    } else if (std::fabs(std::fabs((*d0)[0] * (*d)[0] + (*d0)[1] * (*d)[1] + (*d0)[2] * (*d)[2]) - 1.0) >
               1e-9) {
      return {};
    }
  }
  if (!d0) d0 = P3{0, 0, 1};
  try {
    if (std::fabs(std::fabs((*d0)[2]) - 1.0) < 1e-12) return hullVerticalTranslates(components);

    // Conjugate: rotate the shared axis onto +Z, solve, rotate back.
    const double nx = (*d0)[1], ny = -(*d0)[0];
    const double norm = std::hypot(nx, ny);
    const gp_Ax1 axis(gp_Pnt(0, 0, 0), gp_Dir(nx / norm, ny / norm, 0));
    const double angle = std::acos(std::clamp((*d0)[2], -1.0, 1.0));
    gp_Trsf rotate, back;
    rotate.SetRotation(axis, angle);
    back.SetRotation(axis, -angle);
    std::vector<TopoDS_Shape> rotated;
    for (const auto& s : components)
      rotated.push_back(BRepBuilderAPI_Transform(s, rotate, true).Shape());
    auto result = hullVerticalTranslates(rotated);
    if (result.IsNull()) return {};
    return BRepBuilderAPI_Transform(result, back, true).Shape();
  } catch (const Standard_Failure&) {
    return {};
  }
}

}  // namespace OcctRevolvedHull
