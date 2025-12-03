#ifndef PARTITIONING_H_
#define PARTITIONING_H_
#include <array>
#include <mpi.h>
#include <vector>

struct Settings;

namespace Partitioning {
struct MPIInfo
{
  int rank;
  int size;
  int top_neighbor;
  int bottom_neighbor;
  int left_neighbor;
  int right_neighbor;
  int nCells[2];
  int nCellsWithGhostcells[2];

  inline std::array<std::array<int, 2>, 4> neighbours() const
  {
    return { std::array<int, 2> { top_neighbor, 1 }, { bottom_neighbor, 0 }, { left_neighbor, 3 }, { right_neighbor, 2 } };
  }
  int Partitions[2];
};
void setMPIInfo(MPIInfo& mpiInfo, const Settings& settings, int rank, int size);
const std::vector<MPIInfo>& getInfos();
const MPIInfo& getInfo(size_t x, size_t y);

}

#endif // PARTITIONING_H_
