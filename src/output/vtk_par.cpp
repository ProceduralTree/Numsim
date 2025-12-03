#include "vtk_par.h"
#include "../pde/system.h"
#include "utils/partitioning.h"
#include "utils/settings.h"
#include <cstddef>
#include <grid/grid.h>
#include <mpi.h>
#include <vtkDoubleArray.h>
#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkXMLImageDataWriter.h>

namespace vtk_par {
static auto vtkWriter_ = vtkSmartPointer<vtkXMLImageDataWriter>::New();
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
#if 0
    for (size_t srcY = srcR.begin.y, dstY = dstR.begin.y; srcY <= srcR.end.y; srcY++, dstY++)
    {
      for (size_t srcX = srcR.begin.x, dstX = dstR.begin.x; srcX <= srcR.end.x; srcX++, dstX++)
      {
        getAt(dstX, dstY) = getAt(srcX, srcY);
      }
    }
#endif
    size_t length = srcR.end.x - srcR.begin.x + 1;
    for (size_t srcY = srcR.begin.y, dstY = dstR.begin.y; srcY <= srcR.end.y; srcY++, dstY++)
    {
      // DebugF("copy for rank {} from row: {}", Settings::get().mpi.rank, srcY);
      memcpy(&getAt(dstR.begin.x, dstY), &src[{ srcR.begin.x, static_cast<uint16_t>(srcY) }], length * sizeof(double));
    }
  }
  inline double& getAt(size_t x, size_t y)
  {
    return _data[x + sizeI.x * y];
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
void set_filename(const vtkSmartPointer<vtkXMLImageDataWriter> writer, int& fileNumber)
{

  // Assemble the filename
  std::stringstream fileName;
  fileName << "out/output_" << std::setw(4) << setfill('0') << fileNumber << "."
           << writer->GetDefaultFileExtension();
  checkForDir("out");
  // increment file no.
  // assign the new file name to the output vtkWriter_
  writer->SetFileName(fileName.str().c_str());
  fileNumber++;
}
vtkSmartPointer<vtkImageData> initialize_dataset(const PDESystem& system)
{
  vtkSmartPointer<vtkImageData> dataSet = vtkSmartPointer<vtkImageData>::New();
  dataSet->SetOrigin(0, 0, 0);

  // set spacing of mesh
  const double dx = system.h.x;
  const double dy = system.h.y;
  const double dz = 1;
  dataSet->SetSpacing(dx, dy, dz);

  // set number of points in each dimension, 1 cell in z direction
  dataSet->SetDimensions(
    system.settings.nCells[0] + 1, system.settings.nCells[1] + 1, 1); // we want to have points at each corner of each cell

  return dataSet;
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
}

std::pair<Range, Range> calcCopyRanges(const Grid2D& grid, const Partitioning::MPIInfo& mpi)
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
  if (mpi.left_neighbor < 0)
  {
    srcR.begin.x--;
    dstR.begin.x--;
  }
  if (mpi.bottom_neighbor < 0)
  {
    srcR.begin.y--;
    dstR.begin.y--;
  }
  return std::make_pair(srcR, dstR);
}
void reduceAll(const PDESystem& system, const Partitioning::MPIInfo& mpi)
{
  auto [srcRP, dstRP] = calcCopyRanges(system.p, mpi);
  pressureGrid.copyFromTo(system.p, srcRP, dstRP);
  auto [srcRU, dstRU] = calcCopyRanges(system.u, mpi);
  uGrid.copyFromTo(system.u, srcRU, dstRU);
  auto [srcRV, dstRV] = calcCopyRanges(system.u, mpi);
  vGrid.copyFromTo(system.v, srcRV, dstRV);

  MPI_Reduce(pressureGrid.data(), GlobalpressureGrid.data(), pressureGrid.size(), MPI_DOUBLE, MPI_SUM, root_rank, MPI_COMM_WORLD);
  MPI_Reduce(uGrid.data(), GlobaluGrid.data(), uGrid.size(), MPI_DOUBLE, MPI_SUM, root_rank, MPI_COMM_WORLD);
  MPI_Reduce(vGrid.data(), GlobalvGrid.data(), vGrid.size(), MPI_DOUBLE, MPI_SUM, root_rank, MPI_COMM_WORLD);
}
void writeVTK(const PDESystem& system, double dt)
{
  Partitioning::MPIInfo mpi = system.settings.mpi;
  reduceAll(system, mpi);

  if (mpi.rank != root_rank)
    return;

  static int fileNumber = 0;
  set_filename(vtkWriter_, fileNumber);
  auto dataSet = initialize_dataset(system);

  vtkSmartPointer<vtkDoubleArray> arrayPressure = vtkDoubleArray::New();
  // the pressure is a scalar which means the number of components is 1
  arrayPressure->SetNumberOfComponents(1);
  // Set the number of pressure values and allocate memory for it. We
  // already know the number, it has to be the same as there are nodes in
  // the mesh.
  arrayPressure->SetNumberOfTuples(dataSet->GetNumberOfPoints());
  arrayPressure->SetName("pressure");
  size_t index = 0;
  for (size_t j = 0; j <= system.settings.nCells[1]; j++)
  {
    for (size_t i = 0; i <= system.settings.nCells[0]; i++, index++)
    {
      arrayPressure->SetValue(index, GlobalpressureGrid.interpolate4({ static_cast<uint16_t>(i), static_cast<uint16_t>(j) }));
    }
  }
  assert(index == dataSet->GetNumberOfPoints());
  dataSet->GetPointData()->AddArray(arrayPressure);

  // loop over the nodes of the mesh and assign the interpolated p values
  // in the vtk data structure we only consider the cells that are the
  // actual computational domain, not the helper values in the "halo"

  // now, we should have added as many values as there are points in the
  // vtk data structure

  // add the field variable to the data set

  // add velocity field variable
  // ---------------------------
  vtkSmartPointer<vtkDoubleArray> arrayVelocity = vtkDoubleArray::New();

  // here we have two components (u,v), but ParaView will only allow
  // vector glyphs if we have an ℝ^3 vector, therefore we use a
  // 3-dimensional vector and set the 3rd component to zero
  arrayVelocity->SetNumberOfComponents(3);

  // set the number of values
  arrayVelocity->SetNumberOfTuples(dataSet->GetNumberOfPoints());

  arrayVelocity->SetName("velocity");

  // loop over the mesh where p is defined and assign the values in the
  // vtk data structure
  index = 0; // index for the vtk data structure
  Index I;
  for (size_t j = 0; j <= system.settings.nCells[1]; j++)
  {
    for (size_t i = 0; i <= system.settings.nCells[0]; i++, index++)
    {
      I = { static_cast<uint16_t>(i), static_cast<uint16_t>(j) };
      std::array<double, 3> velocityVector;
      velocityVector[0] = GlobaluGrid.interpolate(I, Iy);
      velocityVector[1] = GlobalvGrid.interpolate(I, Ix);
      velocityVector[2] = 0.0;

      arrayVelocity->SetTuple(index, velocityVector.data());
    }
  }
  // now, we should have added as many values as there are points in the
  // vtk data structure
  assert(index == dataSet->GetNumberOfPoints());

  // add the field variable to the data set
  dataSet->GetPointData()->AddArray(arrayVelocity);

  // add current time
  vtkSmartPointer<vtkDoubleArray> arrayTime = vtkDoubleArray::New();
  arrayTime->SetName("TIME");
  arrayTime->SetNumberOfTuples(1);
  arrayTime->SetTuple1(0, dt);
  dataSet->GetFieldData()->AddArray(arrayTime);

  // Remove unused memory
  dataSet->Squeeze();

  // Write the data
  vtkWriter_->SetInputData(dataSet);

  // vtkWriter_->SetDataModeToAscii();     // comment this in to get ascii
  // text files: those can be checked in an editor
  vtkWriter_->SetDataModeToBinary(); // set file mode to binary files:
                                     // smaller file sizes

  // finally write out the data
  vtkWriter_->Write();
}
/*
// gather all local Grid values into a buffer, holding all values.
// layout is: and should be:
// +---+---+      +---+---+
// |000|000|      |000|111|
// |111|111|      |000|111|
// +---+---+      +---+---+
// |222|222|      |222|333|
// |333|333|      |222|333|
// +---+---+      +---+---+
// so first gatherAll and then unscrambleGrid
void gatherAllNonRoot(const Grid2D& grid, const Range& range)
{
  // TODO: account for ghost border either here or later (which would be faster)
  MPI_Gather(grid.data(), grid.size(), MPI_DOUBLE, nullptr, 0, MPI_DOUBLE, root_rank, MPI_COMM_WORLD);
}
[[nodiscard]] double* gatherAllRoot(const Grid2D& grid, const Range& range)
{
  const auto& settings = Settings::get();

  double* buffer = (double*)malloc(grid.size() * sizeof(double) + ghostBarrierCount * sizeof(double));
  MPI_Gather(grid.data(), grid.size(), MPI_DOUBLE, buffer, range.count(), MPI_DOUBLE, root_rank, MPI_COMM_WORLD);
  return buffer;
}
// unscramble by using temp buffer of size:
// +---+---+
// |000|000|
// |111|111|
// +---+---+
// or precompute swapindices for faster parallel copy but use 2nd grid
void unscrambleGrid(const double* __restrict__ src, double* __restrict__ dst)
{
  const auto& settings = Settings::get();
  const Partitioning::MPIInfo& mpi = *settings.mpi;

  // precompute indices
  std::vector<size_t> blockSizes(mpi.xParts * mpi.yParts);
  for (size_t i = 0; i < mpi.yParts; i++)
    for (size_t j = 0; i < mpi.xParts; j++)
      blockSizes[i * mpi.xParts + j] = Partitioning::getInfos()[i, j].nCells[0];
  std::vector<size_t> prefixes(mpi.xParts * mpi.yParts);
  prefixes[0] = 0;
  for (size_t i = 0; i < mpi.xParts * mpi.yParts; i++)
    prefixes[i + 1] = prefixes[i] + blockSizes[i];

#pragma omp parallel for
  for (size_t j = 0; j < mpi.yParts; j++)
  {
    for (size_t i = 0; i < mpi.xParts; i++)
    {
      size_t oldIndex = i * mpi.xParts + j;
      size_t newIndex = j * mpi.yParts + i;
      memcpy(dst + prefixes[newIndex], src + prefixes[oldIndex], blockSizes[oldIndex] * sizeof(double));
    }
  }
}
[[nodiscard]] double* gatherPartition(const Grid2D& grid, const Range& range)
{
  double* buffer = gatherAllRoot(grid, range);
  double* out = (double*)malloc(grid.range.count() * sizeof(double));
  unscrambleGrid(buffer, out);
  free(buffer);
  return out;
}
void writeVTK(const PDESystem& system, double dt)
{
  double* pressureBuffer = nullptr;
  double* uBuffer = nullptr;
  double* vBuffer = nullptr;
  if (Settings::get().mpi.rank != root_rank)
  {
    gatherAllNonRoot(system.p, system.p.range);
  } else
  {
    pressureBuffer = gatherPartition(system.p, system.p.range);
  }
  if (Settings::get().mpi.rank != root_rank)
  {
    gatherAllNonRoot(system.u, system.p.range);
  } else
  {
    uBuffer = gatherPartition(system.p, system.p.range);
  }
  if (Settings::get().mpi.rank != root_rank)
  {
    gatherAllNonRoot(system.u, system.p.range);
  } else
  {
    vBuffer = gatherPartition(system.p, system.p.range);
  }
  static int fileNumber = 0;
  set_filename(vtkWriter_, fileNumber);
  auto dataSet = initialize_dataset(system);

  // initialize data set that will be output to the file
  dataSet->SetOrigin(0, 0, 0);

  // set spacing of mesh
  //
  vtkSmartPointer<vtkDoubleArray> arrayPressure = vtkDoubleArray::New();
  arrayPressure->SetArray(pressureBuffer, 5, 1);
  // the pressure is a scalar which means the number of components is 1
  arrayPressure->SetNumberOfComponents(1);
  // Set the number of pressure values and allocate memory for it. We
  // already know the number, it has to be the same as there are nodes in
  // the mesh.
  arrayPressure->SetNumberOfTuples(dataSet->GetNumberOfPoints());
  arrayPressure->SetName("pressure");
  double index = write_pressure(arrayPressure, system);
  dataSet->GetPointData()->AddArray(arrayPressure);
  dataSet->GetPointData()->AddArray(pressureBuffer);
  assert(index == dataSet->GetNumberOfPoints());

  // loop over the nodes of the mesh and assign the interpolated p values
  // in the vtk data structure we only consider the cells that are the
  // actual computational domain, not the helper values in the "halo"

  // now, we should have added as many values as there are points in the
  // vtk data structure

  // add the field variable to the data set

  // add velocity field variable
  // ---------------------------
  vtkSmartPointer<vtkDoubleArray> arrayVelocity = vtkDoubleArray::New();

  // here we have two components (u,v), but ParaView will only allow
  // vector glyphs if we have an ℝ^3 vector, therefore we use a
  // 3-dimensional vector and set the 3rd component to zero
  arrayVelocity->SetNumberOfComponents(3);

  // set the number of values
  arrayVelocity->SetNumberOfTuples(dataSet->GetNumberOfPoints());

  arrayVelocity->SetName("velocity");

  // loop over the mesh where p is defined and assign the values in the
  // vtk data structure
  index = 0; // index for the vtk data structure
  Index I;
  arrayVelocity.for (int j = system.begin.y - 1; j <= system.end.y; j++)
  {
    for (int i = system.begin.x - 1; i <= system.end.x; i++, index++)
    {
      I = { static_cast<uint16_t>(i), static_cast<uint16_t>(j) };
      std::array<double, 3> velocityVector;
      velocityVector[0] = interpolate_u(system, system.u, I);
      velocityVector[1] = interpolate_v(system, system.v, I);
      velocityVector[2] = 0.0;

      arrayVelocity->SetTuple(index, velocityVector.data());
    }
  }
  // now, we should have added as many values as there are points in the
  // vtk data structure
  assert(index == dataSet->GetNumberOfPoints());

  // add the field variable to the data set
  dataSet->GetPointData()->AddArray(arrayVelocity);

  // add current time
  vtkSmartPointer<vtkDoubleArray> arrayTime = vtkDoubleArray::New();
  arrayTime->SetName("TIME");
  arrayTime->SetNumberOfTuples(1);
  arrayTime->SetTuple1(0, time);
  dataSet->GetFieldData()->AddArray(arrayTime);

  // Remove unused memory
  dataSet->Squeeze();

  // Write the data
  vtkWriter_->SetInputData(dataSet);

  // vtkWriter_->SetDataModeToAscii();     // comment this in to get ascii
  // text files: those can be checked in an editor
  vtkWriter_->SetDataModeToBinary(); // set file mode to binary files:
                                     // smaller file sizes

  // finally write out the data
  vtkWriter_->Write();
}
*/
}
