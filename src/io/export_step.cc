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

#include "io/export_step.h"

#include <array>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <random>
#include <string>

#ifndef _WIN32
#include <sys/wait.h>
#endif

#include "utils/printutils.h"

namespace fs = std::filesystem;

namespace {

std::string shellQuote(const std::string& s)
{
#ifdef _WIN32
  return "\"" + s + "\"";
#else
  std::string out = "'";
  for (const char c : s) {
    if (c == '\'') {
      out += "'\\''";
    } else {
      out += c;
    }
  }
  return out + "'";
#endif
}

FILE *openPipe(const std::string& command)
{
#ifdef _WIN32
  return _popen(command.c_str(), "r");
#else
  return popen(command.c_str(), "r");
#endif
}

int closePipe(FILE *pipe)
{
#ifdef _WIN32
  return _pclose(pipe);
#else
  const int status = pclose(pipe);
  if (status == -1) return -1;
  return WIFEXITED(status) ? WEXITSTATUS(status) : -1;
#endif
}

fs::path temporaryCsgPath(const fs::path& outputPath)
{
  std::random_device rd;
  std::uniform_int_distribution<unsigned> dist(0, 0xffffff);
  char suffix[16];
  std::snprintf(suffix, sizeof(suffix), "-%06x.csg", dist(rd));
  return fs::temp_directory_path() / ("openscad-" + outputPath.stem().string() + suffix);
}

}  // namespace

bool export_step_external(const std::string& csgText, const fs::path& outputPath,
                          const fs::path& workDir, const std::string& command)
{
  if (command.empty()) {
    LOG(message_group::Error,
        "No STEP export command configured, check Preferences -> Advanced -> STEP Export");
    return false;
  }

  const auto csgPath = temporaryCsgPath(outputPath);
  {
    std::ofstream csgFile(csgPath);
    if (!csgFile.is_open()) {
      LOG(message_group::Error, "Can't write temporary CSG file \"%1$s\"", csgPath.string());
      return false;
    }
    csgFile << csgText << "\n";
  }

#ifdef _WIN32
  const std::string commandLine = "cd /d " + shellQuote(workDir.string()) + " && " + command + " " +
                                  shellQuote(csgPath.string()) + " -o " +
                                  shellQuote(fs::absolute(outputPath).string()) + " 2>&1";
#else
  const std::string commandLine = "cd " + shellQuote(workDir.string()) + " && " + command + " " +
                                  shellQuote(csgPath.string()) + " -o " +
                                  shellQuote(fs::absolute(outputPath).string()) + " 2>&1";
#endif
  PRINTD("Executing: " + commandLine);

  FILE *pipe = openPipe(commandLine);
  if (!pipe) {
    LOG(message_group::Error, "Could not start STEP export command '%1$s'", command);
    fs::remove(csgPath);
    return false;
  }

  std::array<char, 4096> buffer{};
  std::string line;
  while (std::fgets(buffer.data(), buffer.size(), pipe)) {
    line += buffer.data();
    if (!line.empty() && line.back() == '\n') {
      line.pop_back();
      LOG("%1$s", line);
      line.clear();
    }
  }
  if (!line.empty()) LOG("%1$s", line);

  const int exitCode = closePipe(pipe);
  fs::remove(csgPath);

  if (exitCode != 0) {
    LOG(message_group::Error, "STEP export command '%1$s' failed with exit code %2$d", command,
        exitCode);
    return false;
  }
  if (!fs::exists(outputPath)) {
    LOG(message_group::Error, "STEP export command '%1$s' produced no file \"%2$s\"", command,
        outputPath.string());
    return false;
  }
  return true;
}
