#ifndef PARTITIONING_H_
#define PARTITIONING_H_
#include <mpi.h>

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
};

#endif // PARTITIONING_H_
