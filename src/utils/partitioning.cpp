#include "utils/partitioning.h"
#include <cmath>
#include <mpi.h>
#include <utils/settings.h>

namespace Partitioning {
void setMPIInfo(MPIInfo& mpiInfo, const Settings& settings, int rank, int size)
{
  mpiInfo.rank = rank;
  mpiInfo.size = size;

  int nCellsX = settings.nCells[0];
  int nCellsY = settings.nCells[1];

  int px = 0;
  int py = 0;
  double diff = INFINITY;
  for (int i = 1; i <= size; i++)
  {
    double new_diff = std::abs((i / (size / i)) - (nCellsX / nCellsY));
    if (new_diff < diff && size % i == 0)
    {
      diff = new_diff;
      px = i;
      py = size / i;
    }
  };

  int mpiCellsx;
  int mpiCellsy;

  if (nCellsX % px == 0 && nCellsY % py == 0)
  {
    mpiCellsx = nCellsX / px;
    mpiCellsy = nCellsY / py;
  } else if (nCellsY % py != 0 && rank >= (py - 1) * px && (rank + 1) % px != 0)
  {
    mpiCellsx = std::floor(nCellsX / px);
    mpiCellsy = (nCellsY / py) + nCellsY % py;
  } else if (nCellsX % px != 0 && (rank + 1) % px == 0 && rank <= (py - 1) * px)
  {
    mpiCellsx = (nCellsX / px) + nCellsX % px;
    mpiCellsy = std::floor(nCellsY / py);
  } else if ((rank + 1) % px == 0 && rank >= (py - 1) * px)
  {
    mpiCellsx = (nCellsX / px) + nCellsX % px;
    mpiCellsy = (nCellsY / py) + nCellsY % py;
  } else
  {
    mpiCellsx = std::floor(nCellsX / px);
    mpiCellsy = std::floor(nCellsY / py);
  }

  mpiInfo.nCells[0] = mpiCellsx;
  mpiInfo.nCells[1] = mpiCellsy;

  mpiInfo.nCellsWithGhostcells[0] = mpiCellsx + 2;
  mpiInfo.nCellsWithGhostcells[1] = mpiCellsy + 2;

  int lefneighbor;
  int rightneighbor;
  int topneighbor;
  int bottomneighbor;

  if (py == 1)
  {
    topneighbor = -2;
    bottomneighbor = -4;
    if (rank == 0)
    {
      lefneighbor = -1;
      if (size == 1)
      {
        rightneighbor = -3;
      } else
      {
        rightneighbor = rank + 1;
      }
    } else if (rank < px - 1)
    {
      lefneighbor = rank - 1;
      rightneighbor = rank + 1;
    } else
    {
      lefneighbor = rank - 1;
      rightneighbor = -3;
    }
  } else if (rank < px - 1)
  {
    topneighbor = -2;
    bottomneighbor = rank + px;
    if (rank == 0)
    {
      lefneighbor = -1;
      if (px > 1)
      {
        rightneighbor = rank + 1;
      } else
      {
        rightneighbor = -3;
      }
    } else if (rank < px - 1)
    {
      lefneighbor = rank - 1;
      rightneighbor = rank + 1;
    } else
    {
      lefneighbor = rank - 1;
      rightneighbor = -3;
    }
  } else if (rank < px * (py - 1))
  {
    topneighbor = rank - px;
    bottomneighbor = rank + px;
    if (px > 1)
    {
      if (rank % px == 0)
      {
        lefneighbor = -1;
        rightneighbor = rank + 1;
      } else if ((rank + 1) % px == 0)
      {
        lefneighbor = rank - 1;
        rightneighbor = -3;
      } else
      {
        lefneighbor = rank - 1;
        rightneighbor = rank + 1;
      }
    } else
    {
      lefneighbor = -1;
      rightneighbor = -3;
    }
  } else
  {
    topneighbor = rank - px;
    bottomneighbor = -4;
    if (px > 1)
    {
      if (rank % px == 0)
      {
        lefneighbor = -1;
        rightneighbor = rank + 1;
      } else if ((rank + 1) % px == 0)
      {
        lefneighbor = rank - 1;
        rightneighbor = -3;
      } else
      {
        lefneighbor = rank - 1;
        rightneighbor = rank + 1;
      }
    } else
    {
      lefneighbor = -1;
      rightneighbor = -3;
    }
  }

  // boundary conditions for -1 to -4 clockwise: left, top, right, bottom
  mpiInfo.left_neighbor = lefneighbor;
  mpiInfo.right_neighbor = rightneighbor;
  mpiInfo.top_neighbor = topneighbor;
  mpiInfo.bottom_neighbor = bottomneighbor;
  mpiInfo.Partitions[0] = px;
  mpiInfo.Partitions[1] = py;
}
const std::vector<MPIInfo>& getInfos()
{
  static std::vector<MPIInfo> infos;

  // Determine current MPI world size (fall back to 1 if MPI isn't initialized)
  int initialized = 0;
  MPI_Initialized(&initialized);
  int currentSize = 1;
  if (initialized)
  {
    MPI_Comm_size(MPI_COMM_WORLD, &currentSize);
  }

  // If the cached vector doesn't match the current number of ranks, rebuild it
  if ((int)infos.size() != currentSize)
  {
    infos.clear();
    infos.reserve(currentSize);

    const Settings& settings = Settings::get();
    for (int rank = 0; rank < currentSize; ++rank)
    {
      MPIInfo info; // zero-initialized
      setMPIInfo(info, settings, rank, currentSize);
      infos.push_back(std::move(info));
    }
  }

  return infos;
}
// Get MPIInfo for a specific (x,y) partition starting with (0,0) at top-left
const MPIInfo& getInfo(size_t x, size_t y)
{
  const auto& infos = getInfos();
  size_t index = y * infos[0].nCells[0] + x;
  return infos[index];
}
} // namespace Partitioning
