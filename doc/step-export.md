# Native STEP export

This describes the design of the built-in STEP exporter, the choices it
makes and why, and what it does not do. It is written for review: each
section says what was decided and what the evidence was, so a maintainer
can disagree with the decision rather than guess at it.

The feature is a compile-time option, `ENABLE_OCCT`, off by default. With
it off nothing here is built and the binary is unchanged.

## Why a second evaluator

STEP is a boundary-representation format: a cylinder is a cylinder, not a
fan of triangles. That is the whole point of exporting it. Every other
exporter OpenSCAD has receives a `Geometry` that has already been meshed,
and by that time the curvature is gone. Only CSG export sees the node tree.

So the STEP exporter takes the tree path. It walks the evaluated
`AbstractNode` tree and rebuilds it in OpenCASCADE, which means it is a
second geometry evaluator living beside CGAL and Manifold, not a format
writer bolted onto the existing one.

A model exported this way opens in FreeCAD, Fusion or SolidWorks as
solids whose faces can be selected, filleted and dimensioned.

## What it covers, and what happens when it cannot

Leaves: `cube`, `sphere`, `cylinder`, `polyhedron`, `square`, `circle`,
`polygon`, `text`. Operations: `union`, `difference`, `intersection`,
transforms, `linear_extrude`, `rotate_extrude`, and the closed forms of
`hull` and `minkowski` (equal spheres, two unequal spheres, parallel
cylinders, polyhedra, discs, polygons, Minkowski with a ball, and
revolved translates).

Anything else -- `offset`, `projection`, `resize`, extrusion with twist,
`roof`, and hulls or Minkowski sums outside those forms -- is rendered
through OpenSCAD's own evaluator and sewn into a solid. Each such region
is warned about and counted, and the warning names the node and says that
surfaces there are not analytic.

The intent is a guarantee rather than a best effort: **STEP export should
never be worse than STL export of the same model**. Where an exact form
exists it is used; where it does not, the output is what OpenSCAD itself
would have rendered.

## Faceting

`circle($fn=6)` is a hexagon, and a user who wrote that meant a hexagon.
A cylinder at default fineness is a circle that happens to be drawn with
a few segments, and a user who wrote that meant a cylinder.

The exporter facets a curve when its segment count is both deliberate and
below a threshold, default 20 and configurable. Deliberate means an
explicit `$fn`, or `$fa`/`$fs`/`$fe` set away from their defaults. Default
fineness always produces exact curves.

A consequence worth stating plainly: the volume of an exported STEP is
legitimately a few percent above the volume of a default-fineness mesh of
the same model, because an exact cylinder is larger than the prism
inscribed in it. On a corpus sample this accounted for every difference
between 2% and 6% that was not a defect.

## Colour

The requirement was that a STEP file look like the preview. `color()`
follows the preview's rule, where the outermost `color()` wins, and the
bodies are grouped by colour in the file so a region can be selected or
hidden as a unit.

Geometry never depends on colour. Where two differently coloured bodies
overlap, the later operand claims the contested material; a body used as
a cutter never paints what it cuts. Where a colour-preserving boolean
cannot be done, the affected bodies merge as uncoloured and the loss is
local to that collision rather than to the model.

## Trusting OpenCASCADE

OpenCASCADE occasionally returns a valid, self-consistent solid that is
not the answer. This is the part of the design that took the most work,
because such a result passes every obvious check.

Each operation is therefore checked against statements about the answer,
which need no reference render:

- A boolean's extent must lie inside the bounds its inputs imply.
- A union may not come back in more pieces than it was given.
- A union must still contain every operand. *(A five-way union of a hand's
  digits kept half of one finger; a union dropped a plate wedged between a
  ring and a wedge.)*
- A cut must not leave a body inside one of its tools, and must keep the
  argument material that no tool covers. *(A cut through a coil's base
  returned half a revolution out of four.)*
- A `UnifySameDomain` pass is adopted only if it conserves the extent and
  does not loosen tolerances. It also runs on a copy, because it updates
  the tolerances of the shape it is given even when its result is
  discarded.

The containment checks sample up to eight interior points per body, so
they are detectors, not proofs: a body partly swallowed can still pass.

A rejected result is retried with fuzzy tolerances sized as fractions of
the model, and then by folding one body at a time.

When nothing works, the policy differs by operation, and the difference
is empirical rather than aesthetic:

- A **union** hands back its operands unjoined. That is the right material
  in the right places, merely unmerged.
- A **cut or intersection** falls back to a mesh, because its failures are
  wrong rather than unmerged.

Sending unions down the mesh path as well was tried and measured on 500
corpus models: three models were repaired and five went from agreeing
with the mesh to between 50% and 76% off, with two new timeouts and
double the export time. Hence the split.

## Seams

A face that wraps all the way around a periodic surface closes on a seam
edge. Such a face is legal in STEP but fragile to read back: an ellipsoid
fused to a cylinder wrote a file whose ellipsoid came back with no volume
at all, in OpenCASCADE's own reader and in an unrelated viewer.

Closed faces are therefore split in two before writing, which removes the
seam at a cost of a couple of faces per curved body. The split is
geometrically lossless: sampled against the exact ellipsoid the deviation
is 9e-16 before and after. It is rejected if it moves the extent.

## Checking the written file

Every export reads its own file back and compares it with the geometry it
was given, warning if they differ by more than half a percent. This is the
check that would have caught the seam defect above without a user
reporting it.

It is not free. Over a 500-model corpus it added 14% to total export time,
and on the largest files it costs tens of seconds, because reading a STEP
back is comparable work to writing it. Whether that is the right default,
or whether it belongs behind a setting, is a fair question for review.

## Time budget

There is an optional wall-clock budget for the B-rep work, off by default.

It cannot be enforced by checking the clock between operations, because a
single boolean is what runs for minutes. It is expressed as an
OpenCASCADE progress indicator whose `UserBreak` turns true past the
deadline, which the boolean algorithms poll. Abandoned work takes the
same mesh fallback as any other failure, so a model that otherwise never
finishes still produces a file: a corpus gear that ran past three minutes
exports in 44 seconds under a ten second budget.

The mesh fallback is exempt from the budget, since it is what the budget
degrades to.

It is off by default because a budget makes the output depend on the
machine, and a file people diff and archive should not change with the
speed of the computer that wrote it.

## The dependency

OpenCASCADE is a large dependency and the honest numbers are these, from a
macOS universal build:

| | Upstream snapshot | With OpenCASCADE |
|---|---|---|
| Disk image | 70.1 MB | 107.1 MB |
| OpenCASCADE libraries, compressed | - | 33.4 MB |
| Launch to `--version` | 36 ms | 73 ms |
| Resident memory at launch | 23.4 MB | 31.3 MB |

Twenty-six libraries are linked, of which about 11 MB is visualisation
toolkits pulled in transitively by `TKXCAF`. With `ENABLE_OCCT=OFF` the
binary contains no OpenCASCADE symbols and links none of them.

If the launch cost or the download size is the sticking point, the
evaluator is already isolated behind a small interface and could be built
as a module loaded on first use, at the cost of packaging work on three
platforms and a Windows export problem. That has not been done, because
it should be a maintainer's call.

Platforms: enabled on macOS and on Linux where OpenCASCADE 7.6 or newer
ships a CMake config. Disabled on Windows, where the msys2 package
segfaults inside its own `Extrema_ExtCC` under GCC 16 before any of our
code runs.

## Testing

The exporter is covered inside the existing harness, not beside it:

- 41 regression models under `tests/data/scad/step/`, compared through
  `--export-format step-metrics`, a JSON summary of extent, bounding box,
  centroid, face and solid counts, volume per colour, and the same again
  after a write and read-back. Values are rounded so they are stable
  across platforms, and compared exactly.
- 12 Catch2 unit tests tagged `[occt]` over the boolean invariants, the
  mesh conversion and the closed-form hulls.

Both run under `ctest` with everything else.

Beyond the tree, the evaluator has been run over a 500-model sample of a
public corpus of OpenSCAD files, comparing every STEP against the mesh
OpenSCAD renders from the same source.

## Known limitations

- Three corpus models out of 500 do not finish within three minutes. Two
  are slow inside OpenCASCADE's sewing in the mesh fallback, which the
  time budget deliberately does not interrupt.
- A thread built from tilted half-tori exports as a mesh: the cut that
  flattens its base fails and the region degrades. The geometry has an
  exact form; OpenCASCADE will not produce it.
- The exporter runs on the GUI thread with no progress or cancel.
- Layers are written but few consumers read them.

## Questions for maintainers

1. Is an optional dependency of this size acceptable at all, and if so
   should it be a loadable module rather than a link-time option?
2. Is the faceting rule the right one, in particular that default
   `$fa`/`$fs` produce exact curves?
3. Should an unverifiable boolean fall back to a mesh, as a cut does now,
   or should it fail the export outright and say so?
4. Does the colour model match what you would want a STEP consumer to see?
