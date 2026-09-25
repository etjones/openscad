// Minimal reproduction for the msys2 OpenCASCADE segfault: a point
// classified against a plain sphere, then a sphere fused with a box.
// Both crash inside Extrema_ExtCC::Perform in the 7.9.3-3 package.
#include <BRepAlgoAPI_Fuse.hxx>
#include <BRepClass3d_SolidClassifier.hxx>
#include <BRepGProp.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeSphere.hxx>
#include <GProp_GProps.hxx>
#include <iostream>

int main()
{
  const TopoDS_Shape sphere = BRepPrimAPI_MakeSphere(gp_Pnt(0, 0, 0), 5.0).Shape();
  std::cout << "sphere built" << std::endl;
  BRepClass3d_SolidClassifier classifier(sphere, gp_Pnt(1, 1, 1), 1e-7);
  std::cout << "classified: state " << classifier.State() << std::endl;
  const TopoDS_Shape box = BRepPrimAPI_MakeBox(gp_Pnt(3, -5, -5), 10, 10, 10).Shape();
  BRepAlgoAPI_Fuse fuse(sphere, box);
  fuse.Build();
  GProp_GProps props;
  BRepGProp::VolumeProperties(fuse.Shape(), props);
  std::cout << "fused volume " << props.Mass() << " (expect ~1339)" << std::endl;
  return 0;
}
