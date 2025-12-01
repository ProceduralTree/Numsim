#ifndef PARTITIONING_H_
#define PARTITIONING_H_
#include <array>
#include <mpi.h>
#include <utils/settings.h>

struct MPIInfo
{
  int rank;
  int size;
  int Top_neighbor;
  int bottom_neighbor;
  int left_neighbor;
  int right_neighbor;
  int nCells[2];
  int nCellsWithGhostcells[2];

  inline std::array<std::array<int, 2>, 4> neighbours()
  {
    return { std::array<int, 2> { Top_neighbor, 1 }, { bottom_neighbor, 0 }, { left_neighbor, 3 }, { right_neighbor, 2 } };
  }
};
void setMPIInfo(MPIInfo& mpiInfo, const Settings& settings, int rank, int size);

#endif // PARTITIONING_H_
