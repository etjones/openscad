#include <APIHeaderSection_MakeHeader.hxx>
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
