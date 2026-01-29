#ifndef DERIVATIVES_H_
#define DERIVATIVES_H_

#include "grid/sparsegrid.h"
#include <cassert>
#include <grid/grid.h>
#include <pde/system.h>
#include <utils/index.h>
template <typename T>
concept Grid = std::same_as<T, Grid2D> || std::same_as<T, SparseGrid2D<double>>;

#define ASSERT(condition, message)                               \
  do                                                             \
  {                                                              \
    if (!(condition))                                            \
    {                                                            \
      std::cerr << "Assertion failed: " << message << std::endl; \
      assert(condition);                                         \
    }                                                            \
  } while (0)

template <Grid G>
inline double d(Offset Direction, const G& field, Index I, double h)
{
  assert(Direction.x <= I.x + 1);
  assert(Direction.y <= I.y + 1);
  return 1 / h * (field[I + Direction] - field[I]);
}
template <Grid G>
inline double dd(Offset Direction, const G& field, Index I, double h_squared)
{
  assert(Direction.x <= I.x);
  assert(Direction.y <= I.y);
  return 1 / h_squared * (field[I + Direction] + field[I - Direction] - 2 * field[I]);
}
template <Grid G>
inline double duv(Offset Direction, const G& field1, const G& field2, Index I, double h, double alpha)
{
  assert(Direction.x <= I.x);
  assert(Direction.y <= I.y);
  if (Direction == Ix)
  {
    double donor_cell_correction = alpha * (1 / h) * ((std::abs(field1[I + Iy] + field1[I]) * (field2[I] - field2[I + Ix])) / 4 - (std::abs(field1[I - Ix] + field1[I - Ix + Iy]) * (field2[I - Ix] - field2[I])) / 4);
    return (1 / h) * (((field1[I + Iy] + field1[I]) * (field2[I + Ix] + field2[I])) / 4 - ((field1[I - Ix] + field1[I - Ix + Iy]) * (field2[I] + field2[I - Ix])) / 4) + donor_cell_correction;
  } else if (Direction == Iy)
  {
    double donor_cell_correction = alpha * (1 / h) * ((std::abs(field2[I + Ix] + field2[I]) * (field1[I] - field1[I + Iy])) / 4 - (std::abs(field2[I - Iy] + field2[I - Iy + Ix]) * (field1[I - Iy] - field1[I])) / 4);
    return (1 / h) * (((field1[I + Iy] + field1[I]) * (field2[I + Ix] + field2[I])) / 4 - ((field1[I] + field1[I - Iy]) * (field2[I - Iy] + field2[I - Iy + Ix])) / 4) + donor_cell_correction;
  } else
  {
    assert(false && "Invalid Direction for duv");
    return 0.0;
  }
}

template <Grid G>
inline double dxx(Offset Direction, const G& field1, const G& field2, Index I, double h, double alpha)
{
  assert(Direction.x <= I.x);
  assert(Direction.y <= I.y);
  double donor_cell_correction = alpha * (1 / h) * ((std::abs(field1[I + Direction] + field1[I]) * (field2[I] - field2[I + Direction])) / 4 - (std::abs(field1[I - Direction] + field1[I]) * (field2[I - Direction] - field2[I])) / 4);
  return (1 / h) * (((field1[I + Direction] + field1[I]) * (field2[I + Direction] + field2[I])) / 4 - ((field1[I - Direction] + field1[I]) * (field2[I] + field2[I - Direction])) / 4) + donor_cell_correction;
}

#endif // DERIVATIVES_H_
