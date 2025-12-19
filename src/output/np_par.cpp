
#include "np_par.h"
#include "../pde/system.h"
#include "utils/partitioning.h"
#include "utils/settings.h"
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <grid/grid.h>
#include <mpi.h>

namespace np_par {
static std::ofstream outFileP;
static std::ofstream outFileU;
static std::ofstream outFileV;
struct grid
{
  std::vector<double> _data;
  Index beginI;
  Index endI;
  Index sizeI;
  inline void setSize(size_t x, size_t y)
  {
    _data.resize(x * y);
    sizeI = { static_cast<uint16_t>(x), static_cast<uint16_t>(y) };
  }
  inline void copyFromTo(const Grid2D& src, Range srcR, Range dstR)
  {
    assert(dstR.end <= sizeI);
    assert(srcR.count() == dstR.count());
// faster with memcpy but doesnt allow different scaled ranges which we dont need anyways
#if 1
    for (size_t srcY = srcR.begin.y, dstY = dstR.begin.y; srcY <= srcR.end.y; srcY++, dstY++)
    {
      for (size_t srcX = srcR.begin.x, dstX = dstR.begin.x; srcX <= srcR.end.x; srcX++, dstX++)
      {
        getAt(dstX, dstY) = src[{ static_cast<uint16_t>(srcX), static_cast<uint16_t>(srcY) }];
      }
    }
#else
    size_t length = srcR.end.x - srcR.begin.x + 1;
    for (size_t srcY = srcR.begin.y, dstY = dstR.begin.y; srcY <= srcR.end.y; srcY++, dstY++)
    {
      // DebugF("copy for rank {} from row: {}", Settings::get().mpi.rank, srcY);
      memcpy(&getAt(dstR.begin.x, dstY), &src[{ srcR.begin.x, static_cast<uint16_t>(srcY) }], length * sizeof(double));
    }
#endif
  }
  inline double& getAt(size_t x, size_t y)
  {
    return _data[x + sizeI.x * y];
  }
  inline double& operator[](size_t index)
  {
    return _data[index];
  }
  inline double* data() { return _data.data(); }
  inline size_t size() { return _data.size(); }
  inline double interpolate(Index at, Offset offset)
  {
    return (getAt(at.x, at.y) + getAt(at.x + offset.x, at.y + offset.y)) / 2.0;
  }
  inline double interpolate4(Index at)
  {
    return (interpolate(at, Ix) + interpolate(at + Iy, Ix)) / 2.0;
    // return (getAt(at.x, at.y) + getAt(at.x + 1, at.y) + getAt(at.x, at.y + 1) + getAt(at.x + 1, at.y + 1)) / 4.0;
  }
};
constexpr int root_rank = 0;
static grid InterpolatedpressureGrid;
static grid InterpolateduGrid;
static grid InterpolatedvGrid;
static grid GlobalpressureGrid;
static grid GlobaluGrid;
static grid GlobalvGrid;
static grid pressureGrid;
static grid uGrid;
static grid vGrid;

void checkForDir(const std::string& name)
{
  if (std::filesystem::is_directory(name))
    return;
  std::filesystem::create_directory(name);
}
void set_filename(int& fileNumber)
{

  // Assemble the filename
  std::stringstream fileNameP;
  std::stringstream fileNameU;
  std::stringstream fileNameV;
  fileNameP << "out/output_" << std::setw(4) << std::setfill('0') << fileNumber << ".p"
            << ".npy";
  fileNameU << "out/output_" << std::setw(4) << std::setfill('0') << fileNumber << ".u"
            << ".npy";
  fileNameV << "out/output_" << std::setw(4) << std::setfill('0') << fileNumber << ".v"
            << ".npy";
  checkForDir("out");
  // increment file no.
  // assign the new file name to the output vtkWriter_
  outFileP.open(fileNameP.str(), std::ios::binary);
  outFileU.open(fileNameU.str(), std::ios::binary);
  outFileV.open(fileNameV.str(), std::ios::binary);
  fileNumber++;
}
void initializeHeader(const PDESystem& system)
{
  struct __attribute__((packed)) Header
  {
    const unsigned char id = 0x93;
    const char idStr[5] = { 'N', 'U', 'M', 'P', 'Y' };
    const unsigned char majorVersion = 0x01;
    const unsigned char minorVersion = 0x00;
    uint16_t headerLength = 0;
  };
  std::stringstream headerStr;
  headerStr << "descr=numpy.double" << ' ' << "fortran_order=false" << ' ' << "shape=(" << system.settings.nCells[0] << ',' << system.settings.nCells[1] << ')' << "\n";
  size_t headerSize = sizeof(Header) + headerStr.str().size() + 1;
  if (headerSize % 64 != 0)
  {
    headerSize = ((headerSize / 64) + 1) * 64;
    // add padding to headerStr;
    headerStr << "                                                                 ";
  }
  headerSize -= sizeof(Header);
  Header header;
  header.headerLength = headerSize;

  outFileP.write((const char*)&header, sizeof(Header));
  outFileU.write((const char*)&header, sizeof(Header));
  outFileV.write((const char*)&header, sizeof(Header));
  outFileP.write(headerStr.str().c_str(), headerSize);
  outFileU.write(headerStr.str().c_str(), headerSize);
  outFileV.write(headerStr.str().c_str(), headerSize);
}
void init(const PDESystem& system)
{
  const size_t offsetCount = 2;
  pressureGrid.setSize(system.settings.nCells[0] + offsetCount, system.settings.nCells[1] + offsetCount);
  uGrid.setSize(system.settings.nCells[0] + offsetCount, system.settings.nCells[1] + offsetCount);
  vGrid.setSize(system.settings.nCells[0] + offsetCount, system.settings.nCells[1] + offsetCount);
  if (system.settings.mpi.rank != root_rank)
    return;
  GlobalpressureGrid.setSize(system.settings.nCells[0] + offsetCount, system.settings.nCells[1] + offsetCount);
  GlobaluGrid.setSize(system.settings.nCells[0] + offsetCount, system.settings.nCells[1] + offsetCount);
  GlobalvGrid.setSize(system.settings.nCells[0] + offsetCount, system.settings.nCells[1] + offsetCount);
  InterpolatedpressureGrid.setSize(system.settings.nCells[0], system.settings.nCells[1]);
  InterpolateduGrid.setSize(system.settings.nCells[0], system.settings.nCells[1]);
  InterpolatedvGrid.setSize(system.settings.nCells[0], system.settings.nCells[1]);
}

std::pair<Range, Range> calcCopyRanges(const Grid2D& grid, const Partitioning::MPIInfo& mpi, Offset o)
{
  Range srcR = grid.range;

  Range dstR = { { 1, 1 }, {} };
  for (size_t rankX = 0; rankX < mpi.getGridPos().x; rankX++)
  {
    dstR.begin.x += Partitioning::getInfo(rankX, 0).nCells[0];
  }
  for (size_t rankY = mpi.Partitions[1] - 1; rankY > mpi.getGridPos().y; rankY--)
  {
    dstR.begin.y += Partitioning::getInfo(0, rankY).nCells[1];
  }
  dstR.end = dstR.begin + srcR.size() - II;
  if (mpi.top_neighbor < 0)
  {
    srcR.end.y++;
    dstR.end.y++;
  } else
  {
    srcR.end.y -= o.y;
    dstR.end.y -= o.y;
  }
  if (mpi.bottom_neighbor < 0)
  {
    srcR.begin.y--;
    dstR.begin.y--;
  } else
  {
    dstR.begin.y -= o.y;
    dstR.end.y -= o.y;
  }
  if (mpi.left_neighbor < 0)
  {
    srcR.begin.x--;
    dstR.begin.x--;
  } else
  {
    dstR.begin.x -= o.x;
    dstR.end.x -= o.x;
  }
  if (mpi.right_neighbor < 0)
  {
    srcR.end.x++;
    dstR.end.x++;
  } else
  {
    srcR.end.x -= o.x;
    dstR.end.x -= o.x;
  }

  return std::make_pair(srcR, dstR);
}
void reduceAll(const PDESystem& system, const Partitioning::MPIInfo& mpi)
{
#define debugRanges(name, src, dst) DebugF("rank: {}, grid " #name " for src:({},{})({},{}) and dst:({},{})({},{})", mpi.rank, src.begin.x, src.begin.y, src.end.x, src.end.y, dst.begin.x, dst.begin.y, dst.end.x, dst.end.y)
  auto [srcRP, dstRP] = calcCopyRanges(system.p, mpi, Offset { 0, 0 });
  // debugRanges(p, srcRP, dstRP);
  pressureGrid.copyFromTo(system.p, srcRP, dstRP);
  auto [srcRU, dstRU] = calcCopyRanges(system.u, mpi, Ix);
  // debugRanges(u, srcRU, dstRU);
  uGrid.copyFromTo(system.u, srcRU, dstRU);
  auto [srcRV, dstRV] = calcCopyRanges(system.v, mpi, Iy);
  // debugRanges(v, srcRV, dstRV);
  vGrid.copyFromTo(system.v, srcRV, dstRV);

  MPI_Reduce(pressureGrid.data(), GlobalpressureGrid.data(), pressureGrid.size(), MPI_DOUBLE, MPI_SUM, root_rank, MPI_COMM_WORLD);
  MPI_Reduce(uGrid.data(), GlobaluGrid.data(), uGrid.size(), MPI_DOUBLE, MPI_SUM, root_rank, MPI_COMM_WORLD);
  MPI_Reduce(vGrid.data(), GlobalvGrid.data(), vGrid.size(), MPI_DOUBLE, MPI_SUM, root_rank, MPI_COMM_WORLD);
}
void writeNP(const PDESystem& system, double dt)
{
  Partitioning::MPIInfo mpi = system.settings.mpi;
  reduceAll(system, mpi);

  if (mpi.rank != root_rank)
    return;

  static int fileNumber = 0;
  set_filename(fileNumber);
  initializeHeader(system);

  size_t index = 0;
  for (size_t j = 0; j <= system.settings.nCells[1]; j++)
  {
    for (size_t i = 0; i <= system.settings.nCells[0]; i++, index++)
    {
      InterpolatedpressureGrid[index] = GlobalpressureGrid.interpolate4({ static_cast<uint16_t>(i), static_cast<uint16_t>(j) });
    }
  }
  assert(index == dataSet->GetNumberOfPoints());
  index = 0; // index for the vtk data structure
  Index I;
  for (size_t j = 0; j <= system.settings.nCells[1]; j++)
  {
    for (size_t i = 0; i <= system.settings.nCells[0]; i++, index++)
    {
      I = { static_cast<uint16_t>(i), static_cast<uint16_t>(j) };

      InterpolateduGrid[index] = GlobaluGrid.interpolate(I, Iy);
      InterpolatedvGrid[index] = GlobalvGrid.interpolate(I, Ix);
    }
  }
  assert(index == dataSet->GetNumberOfPoints());

  outFileP.write((const char*)InterpolatedpressureGrid.data(), InterpolatedpressureGrid.size());
  outFileU.write((const char*)InterpolateduGrid.data(), InterpolateduGrid.size());
  outFileV.write((const char*)InterpolatedvGrid.data(), InterpolatedvGrid.size());

  outFileP.close();
  outFileU.close();
  outFileV.close();
}
}
