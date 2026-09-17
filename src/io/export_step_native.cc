#include <APIHeaderSection_MakeHeader.hxx>
#include <BRepBndLib.hxx>
#include <BRepGProp.hxx>
#include <Bnd_Box.hxx>
#include <GProp_GProps.hxx>
#include <STEPCAFControl_Reader.hxx>
#include <TDF_LabelSequence.hxx>
#include <TopExp_Explorer.hxx>
#include <Quantity_Color.hxx>
#include <cmath>
#include <map>
#include <random>
#include <IFSelect_ReturnStatus.hxx>
#include <Interface_Static.hxx>
#include <Message.hxx>
#include <Message_Messenger.hxx>
#include <Message_PrinterOStream.hxx>
#include <Quantity_ColorRGBA.hxx>
#include <STEPCAFControl_Writer.hxx>
#include <STEPControl_StepModelType.hxx>
#include <Standard_Failure.hxx>
#include <TCollection_ExtendedString.hxx>
#include <TCollection_HAsciiString.hxx>
#include <TDF_Label.hxx>
#include <TDataStd_Name.hxx>
#include <TDocStd_Document.hxx>
#include <TopLoc_Location.hxx>
#include <XCAFApp_Application.hxx>
#include <XCAFDoc_ColorTool.hxx>
#include <XCAFDoc_ColorType.hxx>
#include <XCAFDoc_DocumentTool.hxx>
#include <XCAFDoc_ShapeTool.hxx>
#include <cstdio>
#include <filesystem>
#include <string>
#include <utility>
#include <vector>

#include "geometry/occt/OcctBoolean.h"
#include "geometry/occt/OcctBuilder.h"
#include "geometry/occt/OcctGeometry.h"
#include "io/export_step.h"
#include "json/json.hpp"
#include "geometry/occt/OcctBridge.h"

namespace fs = std::filesystem;

namespace {

void quietOcct()
{
  const auto messenger = Message::DefaultMessenger();
  for (Message_SequenceOfPrinters::Iterator it(messenger->Printers()); it.More(); it.Next()) {
    it.Value()->SetTraceLevel(Message_Fail);
  }
}

std::string colorLabel(const Color4f& color)
{
  int r, g, b, a;
  color.getRgba(r, g, b, a);
  char buffer[16];
  std::snprintf(buffer, sizeof(buffer), "#%02x%02x%02x", r, g, b);
  return buffer;
}

void setName(const TDF_Label& label, const std::string& text)
{
  TDataStd_Name::Set(label, TCollection_ExtendedString(text.c_str()));
}

TDF_Label addBody(const Handle(XCAFDoc_ShapeTool) & shapes, const Handle(XCAFDoc_ColorTool) & colors,
                  const TDF_Label *parent, const TopoDS_Shape& piece, const Color4f& color)
{
  auto part = shapes->AddShape(piece, false);
  setName(part, color.isValid() ? colorLabel(color) : "solid");
  if (color.isValid()) {
    float r, g, b, a;
    color.getRgba(r, g, b, a);
    const Quantity_ColorRGBA rgba(r, g, b, a);
    // A solid carries a generic color; a sheet has no volume to color, so
    // its color is written as a face style as well.
    colors->SetColor(part, rgba, XCAFDoc_ColorGen);
    colors->SetColor(part, rgba, XCAFDoc_ColorSurf);
  }
  if (parent) shapes->AddComponent(*parent, part, TopLoc_Location());
  return part;
}

struct Piece {
  TopoDS_Shape shape;
  Color4f color;
};

}  // namespace

namespace {

bool writeStep(const OcctGeometry& geometry, const fs::path& outputPath, const std::string& title)
{
  std::vector<Piece> pieces;
  for (const auto& body : geometry.bodies) {
    for (const auto& piece : OcctBoolean::piecesOf(body.shape, geometry.dim)) {
      pieces.push_back({piece, body.color});
    }
  }
  if (pieces.empty()) {
    OcctBridge::error("STEP export: the model produced no geometry");
    return false;
  }

  Handle(TDocStd_Document) doc;
  XCAFApp_Application::GetApplication()->NewDocument("MDTV-XCAF", doc);
  auto shapes = XCAFDoc_DocumentTool::ShapeTool(doc->Main());
  auto colors = XCAFDoc_DocumentTool::ColorTool(doc->Main());

  if (pieces.size() == 1) {
    addBody(shapes, colors, nullptr, pieces.front().shape, pieces.front().color);
  } else {
    auto rootLabel = shapes->NewShape();
    setName(rootLabel, title.empty() ? "model" : title);
    if (geometry.hasColor()) {
      // One group per color, in order of first appearance, so every region
      // of a color can be selected or hidden together.
      std::vector<std::pair<Color4f, TDF_Label>> groups;
      for (const auto& piece : pieces) {
        TDF_Label *group = nullptr;
        for (auto& [color, label] : groups) {
          if ((!color.isValid() && !piece.color.isValid()) ||
              (color.isValid() && piece.color.isValid() && color == piece.color)) {
            group = &label;
            break;
          }
        }
        if (!group) {
          auto label = shapes->NewShape();
          setName(label, piece.color.isValid() ? colorLabel(piece.color) : "uncolored");
          shapes->AddComponent(rootLabel, label, TopLoc_Location());
          groups.emplace_back(piece.color, label);
          group = &groups.back().second;
        }
        addBody(shapes, colors, group, piece.shape, piece.color);
      }
    } else {
      for (const auto& piece : pieces) addBody(shapes, colors, &rootLabel, piece.shape, piece.color);
    }
  }
  shapes->UpdateAssemblies();

  try {
    Interface_Static::SetCVal("write.step.unit", "MM");
    STEPCAFControl_Writer writer;
    writer.SetColorMode(true);
    writer.SetNameMode(true);
    writer.SetLayerMode(true);
    APIHeaderSection_MakeHeader header(writer.ChangeWriter().Model());
    if (header.IsDone()) {
      header.SetOriginatingSystem(new TCollection_HAsciiString("OpenSCAD"));
      if (!title.empty()) header.SetName(new TCollection_HAsciiString(title.c_str()));
    }
    if (!writer.Transfer(doc, STEPControl_AsIs)) {
      OcctBridge::error("STEP export: transfer to the STEP model failed");
      return false;
    }
    if (writer.Write(outputPath.string().c_str()) != IFSelect_RetDone) {
      OcctBridge::error("Can't write STEP file \"" + outputPath.string() + "\"");
      return false;
    }
  } catch (const Standard_Failure& e) {
    OcctBridge::error(std::string("STEP export failed in OpenCASCADE: ") + e.GetMessageString());
    return false;
  }
  return true;
}

double rounded(double v)
{
  const double r = std::round(v * 1e4) / 1e4;
  return r == 0 ? 0.0 : r;
}

nlohmann::json measure(const TopoDS_Shape& shape, unsigned int dim)
{
  nlohmann::json out;
  out["extent"] = rounded(OcctBoolean::extent(shape, dim));
  Bnd_Box box;
  if (!shape.IsNull()) BRepBndLib::AddOptimal(shape, box, false, false);
  if (box.IsVoid()) {
    out["bbox"] = nlohmann::json::array({0.0, 0.0, 0.0});
    out["centroid"] = nlohmann::json::array({0.0, 0.0, 0.0});
  } else {
    double xmin, ymin, zmin, xmax, ymax, zmax;
    box.Get(xmin, ymin, zmin, xmax, ymax, zmax);
    out["bbox"] =
      nlohmann::json::array({rounded(xmax - xmin), rounded(ymax - ymin), rounded(zmax - zmin)});
    GProp_GProps props;
    if (dim == 2) BRepGProp::SurfaceProperties(shape, props);
    else BRepGProp::VolumeProperties(shape, props);
    const auto c = props.CentreOfMass();
    out["centroid"] = nlohmann::json::array({rounded(c.X()), rounded(c.Y()), rounded(c.Z())});
  }
  size_t faces = 0, solids = 0;
  for (TopExp_Explorer it(shape, TopAbs_FACE); it.More(); it.Next()) ++faces;
  for (TopExp_Explorer it(shape, TopAbs_SOLID); it.More(); it.Next()) ++solids;
  out["faces"] = faces;
  out["solids"] = solids;
  return out;
}

std::string rgbLabel(double r, double g, double b)
{
  char buffer[16];
  std::snprintf(buffer, sizeof(buffer), "#%02x%02x%02x", static_cast<int>(std::lround(r * 255)),
                static_cast<int>(std::lround(g * 255)), static_cast<int>(std::lround(b * 255)));
  return buffer;
}

// Every simple shape reachable from the document's free shapes, with the
// generic color on its own label or the referred one.
void collectReadBack(const Handle(XCAFDoc_ShapeTool) & shapes, const Handle(XCAFDoc_ColorTool) & colors,
                     const TDF_Label& label, std::vector<TopoDS_Shape>& all,
                     std::map<std::string, double>& byColor, unsigned int dim)
{
  TDF_Label target = label;
  if (XCAFDoc_ShapeTool::IsReference(label)) XCAFDoc_ShapeTool::GetReferredShape(label, target);
  if (XCAFDoc_ShapeTool::IsAssembly(target)) {
    TDF_LabelSequence components;
    XCAFDoc_ShapeTool::GetComponents(target, components);
    for (int i = 1; i <= components.Length(); ++i) {
      collectReadBack(shapes, colors, components.Value(i), all, byColor, dim);
    }
    return;
  }
  const auto shape = XCAFDoc_ShapeTool::GetShape(label);
  if (shape.IsNull()) return;
  all.push_back(shape);
  // The STEP reader attaches what was written as a surface style as a
  // surface color, so look for that before the generic one.
  Quantity_ColorRGBA rgba;
  std::string key = "uncolored";
  if (colors->GetColor(target, XCAFDoc_ColorSurf, rgba) ||
      colors->GetColor(label, XCAFDoc_ColorSurf, rgba) ||
      colors->GetColor(target, XCAFDoc_ColorGen, rgba) ||
      colors->GetColor(label, XCAFDoc_ColorGen, rgba)) {
    const auto rgb = rgba.GetRGB();
    key = rgbLabel(rgb.Red(), rgb.Green(), rgb.Blue());
  }
  byColor[key] += OcctBoolean::extent(shape, dim);
}

}  // namespace

std::string step_metrics_json(const Tree& tree, const AbstractNode& root, int facetThreshold)
{
  quietOcct();
  OcctBuilder builder(tree, facetThreshold);
  OcctGeometry geometry;
  try {
    geometry = builder.build(root);
  } catch (const Standard_Failure& e) {
    OcctBridge::error(std::string("STEP metrics failed in OpenCASCADE: ") + e.GetMessageString());
    return {};
  }
  nlohmann::json out;
  out["dim"] = geometry.dim;
  out["fallbacks"] = builder.fallbacks().size();
  std::vector<TopoDS_Shape> shapes;
  std::map<std::string, double> byColor;
  for (const auto& body : geometry.bodies) {
    shapes.push_back(body.shape);
    byColor[body.color.isValid() ? colorLabel(body.color) : "uncolored"] +=
      OcctBoolean::extent(body.shape, geometry.dim);
  }
  out["built"] = measure(OcctBoolean::makeCompound(shapes), geometry.dim);
  for (auto& [key, value] : byColor) out["built"]["colors"][key] = rounded(value);

  if (!geometry.isEmpty()) {
    std::random_device rd;
    const auto path =
      fs::temp_directory_path() / ("openscad-step-metrics-" + std::to_string(rd()) + ".step");
    nlohmann::json readBack;
    try {
      if (writeStep(geometry, path, "metrics")) {
        Handle(TDocStd_Document) doc;
        XCAFApp_Application::GetApplication()->NewDocument("MDTV-XCAF", doc);
        STEPCAFControl_Reader reader;
        reader.SetColorMode(true);
        reader.SetNameMode(true);
        if (reader.ReadFile(path.string().c_str()) == IFSelect_RetDone && reader.Transfer(doc)) {
          auto shapeTool = XCAFDoc_DocumentTool::ShapeTool(doc->Main());
          auto colorTool = XCAFDoc_DocumentTool::ColorTool(doc->Main());
          TDF_LabelSequence free;
          shapeTool->GetFreeShapes(free);
          std::vector<TopoDS_Shape> all;
          std::map<std::string, double> colorsBack;
          for (int i = 1; i <= free.Length(); ++i) {
            collectReadBack(shapeTool, colorTool, free.Value(i), all, colorsBack, geometry.dim);
          }
          readBack = measure(OcctBoolean::makeCompound(all), geometry.dim);
          for (auto& [key, value] : colorsBack) readBack["colors"][key] = rounded(value);
        }
      }
    } catch (const Standard_Failure& e) {
      OcctBridge::error(std::string("STEP metrics round trip failed: ") + e.GetMessageString());
    }
    std::error_code ec;
    fs::remove(path, ec);
    out["roundtrip"] = readBack;
  }
  return out.dump(1) + "\n";
}

bool export_step_native(const Tree& tree, const AbstractNode& root, const fs::path& outputPath,
                        int facetThreshold, const std::string& title)
{
  quietOcct();
  OcctBuilder builder(tree, facetThreshold);
  OcctGeometry geometry;
  try {
    geometry = builder.build(root);
  } catch (const Standard_Failure& e) {
    OcctBridge::error(std::string("STEP export failed in OpenCASCADE: ") + e.GetMessageString());
    return false;
  }
  if (geometry.isEmpty()) {
    OcctBridge::error("STEP export: the model produced no geometry");
    return false;
  }

  return writeStep(geometry, outputPath, title);
}
