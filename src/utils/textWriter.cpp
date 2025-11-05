#include "textWriter.h"

#include <fstream>
#include <iomanip>
#include <sstream>

#include <pde/system.h>

namespace Writer {
void writeTextVelocity(const PDESystem& system, double currentTime)
{
  static uint32_t fileNo = 0;
  // Assemble the filename
  std::stringstream fileName;
  fileName << "out/output_" << std::setw(4) << std::setfill('0') << fileNo << ".txt";

  // increment file no.
  fileNo++;

  // open file
  std::ofstream file(fileName.str().c_str());

  if (!file.is_open())
  {
    std::cout << "Could not write to file \"" << fileName.str() << "\".";
    return;
  }

  // write time
  file << "t: " << currentTime << std::endl;

  // write mesh width
  file << "nCells: " << system.settings.nCells[0] << "x" << system.settings.nCells[1]
       << ", dx: " << system.h.y << ", dy: " << system.h.y << std::endl
       << std::endl;

  const int fieldWidth = 9; // number of characters to use for a single value

  auto writeGrid = [&](const Grid2D& grid, const std::string& name) {
    uint16_t sizeX = grid.end.x - grid.begin.x + 2;
    uint16_t sizeY = grid.end.y - grid.begin.y + 2;

    file << name << " (" << sizeX << "x" << sizeY << "): " << std::endl
         << std::string(fieldWidth, ' ') << "|";
    for (int i = -1; i < sizeX - 1; i++)
    {
      file << std::setw(fieldWidth) << i;
    }
    file << std::endl
         << std::string(fieldWidth * (sizeX + 2) + 1, '-') << std::endl;

    // write values
    for (int j = grid.end.y + 1; j >= grid.begin.y - 1; j--)
    {
      file << std::setw(fieldWidth) << j - grid.begin.y - 1 << "|";
      for (int i = grid.begin.x - 1; i <= grid.end.x + 1; i++)
      {
        file << std::setw(fieldWidth) << std::setprecision(fieldWidth - 6) << grid[i, j];
      }
      file << std::endl;
    }
    file << std::endl;
  };

  // write u
  // ---------
  // write header lines
  writeGrid(system.u, "u");

  // write v
  // ---------
  // write header lines
  writeGrid(system.v, "v");

  // write p
  // ---------
  // write header lines
  writeGrid(system.p, "p");

  // write f
  // ---------
  // write header lines
  writeGrid(system.F, "F");

  // write g
  // ---------
  // write header lines
  writeGrid(system.G, "G");

  // write rhs
  // ---------
  // write header lines
  writeGrid(system.rhs, "rhs");
}

void writeTextPressure(const PDESystem& system)
{
  // counter for files, counter value is part of the file name
  static int pressurefileNo = 0;

  // Assemble the filename
  std::stringstream fileName;
  fileName << "out/pressure_" << std::setw(4) << std::setfill('0') << pressurefileNo++ << ".txt";

  // open file
  std::ofstream file(fileName.str().c_str());

  if (!file.is_open())
  {
    std::cout << "Could not write to file \"" << fileName.str() << "\".";
    return;
  }

  // write mesh width
  file << "nCells: " << system.settings.nCells[0] << "x" << system.settings.nCells[1]
       << ", dx: " << system.h.x << ", dy: " << system.h.y << std::endl
       << std::endl;

  const int fieldWidth = 9; // number of characters to use for a single value

  // write p
  // ---------
  // write header lines

  const Grid2D& p = system.p;
  uint16_t sizeX = p.end.x - p.begin.x + 2;
  uint16_t sizeY = p.end.y - p.begin.y + 2;
  file << "p (" << sizeX << "x" << sizeY << "): " << std::endl
       << std::string(fieldWidth, ' ') << "|";
  for (int i = -1; i < sizeX - 1; i++)
  {
    file << std::setw(fieldWidth) << i;
  }
  file << std::endl
       << std::string(fieldWidth * (sizeX + 2) + 1, '-') << std::endl;

  // write p values
  for (int j = p.end.y + 1; j >= p.begin.y - 1; j--)
  {
    file << std::setw(fieldWidth) << j - p.begin.y - 1 << "|";
    for (int i = p.begin.x - 1; i <= p.end.y + 1; i++)
    {
      file << std::setw(fieldWidth) << std::setprecision(fieldWidth - 6) << p[i, j];
    }
    file << std::endl;
  }
  file << std::endl;
}
} // namespace Writer
