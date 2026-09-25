/*
 *  OpenSCAD (www.openscad.org)
 *  Copyright (C) 2009-2019 Clifford Wolf <clifford@clifford.at> and
 *                          Marius Kintel <marius@kintel.net>
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  As a special exception, you have permission to link this program
 *  with the CGAL library and distribute executables, as long as you
 *  follow the requirements of the GNU GPL in regard to all of the
 *  software in the executable aside from CGAL.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

#pragma once

#include <filesystem>
#include <string>

// STEP export via an external converter (scad2step, from the scad123d
// project). The evaluated CSG tree is written to a temporary .csg file
// and handed to the configured command, which rebuilds the analytic
// primitives in a B-rep kernel and writes the STEP file.
//
//   <command> <tmp>.csg -o <outputPath>
//
// Returns true when the command exited successfully and the output file
// exists. Everything the command prints is forwarded to the log.
bool export_step_external(const std::string& csgText, const std::filesystem::path& outputPath,
                          const std::filesystem::path& workDir, const std::string& command);

#ifdef ENABLE_OCCT
class Tree;
class AbstractNode;

// STEP export through the built-in OpenCASCADE evaluator: the node tree
// is rebuilt as B-rep geometry and written as an XDE document, one
// product per body, colored and grouped by color.
// `groupByColor` puts every body of one colour under its own node of the
// file's assembly tree, so a region can be selected or hidden as a unit in
// consumers that ignore layers; off, bodies sit directly under the root
// and carry their colour only as a style, which is what select-by-colour
// uses anyway. The tree is the one grouping every consumer honours, and
// spending it on colour is a default, not a design: see the discussion in
// doc/step-export.md.
//
// `timeBudget` is a wall-clock limit in seconds for the B-rep work, or 0
// for none. Past it OpenCASCADE abandons whatever it is doing and the
// region falls back to a mesh, so a pathological model still produces a
// file. Off by default, so that the same model always exports the same
// geometry whatever the machine.
bool export_step_native(const Tree& tree, const AbstractNode& root,
                        const std::filesystem::path& outputPath, int facetThreshold, int timeBudget,
                        bool groupByColor, const std::string& title);

// Measurements of what the built-in evaluator produces, as JSON, for
// regression tests: volume (area for 2D), bounding box, centroid, face
// and solid counts, volume per color, and the same again after writing
// a STEP file and reading it back. Values are rounded so the output is
// stable across platforms.
std::string step_metrics_json(const Tree& tree, const AbstractNode& root, int facetThreshold,
                              bool groupByColor);
#endif
