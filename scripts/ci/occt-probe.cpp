// Minimal reproduction for the msys2 OpenCASCADE segfault, following the
// exact sequence of OcctBoolean::cut and its invariant: an N-ary cut with
// parallelism off, a unify pass on a copy, then a point classified
// against the *result* -- whose new intersection edges are what the
// classifier's ray-versus-edge extrema chokes on in OpenSCAD's process.
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepAlgoAPI_Fuse.hxx>
#include <BRepBuilderAPI_Copy.hxx>
#include <BRepClass3d_SolidClassifier.hxx>
#include <BRepGProp.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepPrimAPI_MakeSphere.hxx>
#include <GProp_GProps.hxx>
#include <ShapeUpgrade_UnifySameDomain.hxx>
#include <TopTools_ListOfShape.hxx>
#include <cfenv>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#ifdef _WIN32
#include <windows.h>
#include <float.h>
#endif

static TopoDS_Shape boolean(bool cut, const TopoDS_Shape& a, const TopoDS_Shape& b)
{
  TopTools_ListOfShape args, tools;
  args.Append(a);
  tools.Append(b);
  if (cut) {
    BRepAlgoAPI_Cut op;
    op.SetArguments(args);
    op.SetTools(tools);
    op.SetRunParallel(false);
    op.Build();
    return op.Shape();
  }
  BRepAlgoAPI_Fuse op;
  op.SetArguments(args);
  op.SetTools(tools);
  op.SetRunParallel(false);
  op.Build();
  return op.Shape();
}

static double volume(const TopoDS_Shape& s)
{
  GProp_GProps p;
  BRepGProp::VolumeProperties(s, p);
  return p.Mass();
}

static void classify(const char *what, const TopoDS_Shape& s)
{
  GProp_GProps p;
  BRepGProp::VolumeProperties(s, p);
  const gp_Pnt c = p.CentreOfMass();
  std::cout << what << ": volume " << volume(s) << ", classifying centroid..." << std::flush;
  BRepClass3d_SolidClassifier classifier(s);
  classifier.Perform(c, 1e-7);
  std::cout << " state " << classifier.State() << std::endl;
}

// Process-level conditions OpenSCAD imposes that a bare probe does not,
// selected by environment so one binary covers them all:
//   PROBE_ROUND=up|down|zero   rounding mode, as CGAL's interval arithmetic sets it
//   PROBE_X87=up               the x87 control word too, which MinGW keeps separately
//   PROBE_LOAD=a.dll;b.dll     load these DLLs first, so their DllMain runs
static void imposeConditions()
{
  if (const char *r = std::getenv("PROBE_ROUND")) {
    const int mode = !std::strcmp(r, "up") ? FE_UPWARD : !std::strcmp(r, "down") ? FE_DOWNWARD : FE_TOWARDZERO;
    std::cout << "fesetround(" << r << ") -> " << std::fesetround(mode) << std::endl;
  }
#ifdef _WIN32
  if (std::getenv("PROBE_X87")) {
    unsigned int cw = 0;
    _controlfp_s(&cw, _RC_UP, _MCW_RC);
    std::cout << "_controlfp x87 rounding up, cw=" << std::hex << cw << std::dec << std::endl;
  }
  if (const char *l = std::getenv("PROBE_LOAD")) {
    std::string list(l);
    size_t start = 0;
    while (start < list.size()) {
      const size_t end = list.find(';', start);
      const std::string name = list.substr(start, end == std::string::npos ? std::string::npos : end - start);
      HMODULE h = LoadLibraryA(name.c_str());
      std::cout << "LoadLibrary(" << name << ") -> " << (h ? "ok" : "FAILED") << std::endl;
      if (end == std::string::npos) break;
      start = end + 1;
    }
  }
#endif
}

int main()
{
  imposeConditions();
  const TopoDS_Shape cube = BRepPrimAPI_MakeBox(gp_Pnt(-5, -5, -5), 10, 10, 10).Shape();
  const TopoDS_Shape ball = BRepPrimAPI_MakeSphere(gp_Pnt(0, 0, 0), 6.0).Shape();
  const TopoDS_Shape cyl = BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(0, 0, -20), gp_Dir(0, 0, 1)), 2.0, 40.0).Shape();
  std::cout << "primitives built" << std::endl;

  const TopoDS_Shape cut = boolean(true, cube, ball);
  classify("cube - ball", cut);
  ShapeUpgrade_UnifySameDomain unifier(BRepBuilderAPI_Copy(cut).Shape(), true, true, true);
  unifier.Build();
  classify("cube - ball, unified", unifier.Shape());

  const TopoDS_Shape fused = boolean(false, ball, cyl);
  classify("ball + cylinder", fused);
  const TopoDS_Shape cut2 = boolean(true, fused, cube);
  classify("(ball + cylinder) - cube", cut2);
  std::cout << "PROBE COMPLETE" << std::endl;
  return 0;
}
