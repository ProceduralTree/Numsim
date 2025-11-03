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

  auto writeGrid = [&](const Grid2D& grid, const std::string& name, uint16_t offsetX, uint16_t offsetY) {
    uint16_t sizeX = grid.size_x + offsetX;
    uint16_t sizeY = grid.size_y + offsetY;

    file << name << " (" << sizeX << "x" << sizeY << "): " << std::endl
         << std::string(fieldWidth, ' ') << "|";
    for (int i = -1; i < sizeX - 1; i++)
    {
      file << std::setw(fieldWidth) << i;
    }
    file << std::endl
         << std::string(fieldWidth * (sizeX + 2) + 1, '-') << std::endl;

    // write u values
    for (int j = sizeY - 1; j >= -1; j--)
    {
      file << std::setw(fieldWidth) << j << "|";
      for (int i = 0; i < sizeX; i++)
      {
        file << std::setw(fieldWidth) << std::setprecision(fieldWidth - 6) << grid[i, j + 1];
      }
      file << std::endl;
    }
    file << std::endl;
  };

  // write u
  // ---------
  // write header lines
  writeGrid(system.u, "u", -1, 0);

  // write v
  // ---------
  // write header lines
  writeGrid(system.v, "v", 0, -1);

  // write p
  // ---------
  // write header lines
  writeGrid(system.p, "p", 0, 0);

  // write f
  // ---------
  // write header lines
  writeGrid(system.F, "F", -1, 0);

  // write g
  // ---------
  // write header lines
  writeGrid(system.G, "G", 0, -1);

  // write rhs
  // ---------
  // write header lines
  writeGrid(system.b, "rhs", 0, 0);
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
  file << "p (" << system.p.size_x << "x" << system.p.size_y << "): " << std::endl
       << std::string(fieldWidth, ' ') << "|";
  for (int i = -1; i < system.p.size_x - 1; i++)
  {
    file << std::setw(fieldWidth) << i;
  }
  file << std::endl
       << std::string(fieldWidth * (system.p.size_x + 2) + 1, '-') << std::endl;

  // write p values
  for (int j = system.p.size_y - 2; j >= -1; j--)
  {
    file << std::setw(fieldWidth) << j << "|";
    for (int i = 0; i < system.p.size_y; i++)
    {
      file << std::setw(fieldWidth) << std::setprecision(fieldWidth - 6) << system.p[i, j];
    }
    file << std::endl;
  }
  file << std::endl;
}
} // namespace Writer
