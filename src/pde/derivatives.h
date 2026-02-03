#ifndef DERIVATIVES_H_
#define DERIVATIVES_H_

#include "grid/adjacencymap.h"
#include "grid/sparsegrid.h"
#include <cassert>
#include <cstddef>
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

template <Grid G, Offset Direction>
inline double d(const G& field, size_t index, double h, const AdjMap& map)
{
  return 1 / h * (field[get_index<Direction, Sign::Plus>(index, map)] - field[index]);
}
template <Grid G, Offset Direction>
inline double dd(const G& field, size_t index, double h_squared, const AdjMap& map)
{
  return 1 / h_squared * (field[get_index<Direction, Sign::Plus>(index, map)] + field[get_index<Direction, Sign::Minus>(index, map)] - 2 * field[index]);
}
template <Grid G, Offset Direction>
inline double duv(const G& field1, const G& field2, size_t index, double h, double alpha, const AdjMap& map)
{
  size_t top = get_index<Iy, Sign::Plus>(index, map);
  size_t bottom = get_index<Iy, Sign::Minus>(index, map);
  size_t left = get_index<Ix, Sign::Minus>(index, map);
  size_t right = get_index<Ix, Sign::Plus>(index, map);
  size_t top_left = get_index<Iy, Sign::Plus>(left, map);
  size_t bottom_right = get_index<Ix, Sign::Plus>(bottom, map);
  if constexpr (Direction == Ix)
  {
    double donor_cell_correction = alpha * (1 / h) * ((std::abs(field1[top] + field1[index]) * (field2[index] - field2[right])) / 4 - (std::abs(field1[left] + field1[top_left]) * (field2[left] - field2[index])) / 4);
    return (1 / h) * (((field1[top] + field1[index]) * (field2[right] + field2[index])) / 4 - ((field1[left] + field1[top_left]) * (field2[index] + field2[left])) / 4) + donor_cell_correction;
  } else if constexpr (Direction == Iy)
  {
    double donor_cell_correction = alpha * (1 / h) * ((std::abs(field2[right] + field2[index]) * (field1[index] - field1[top])) / 4 - (std::abs(field2[bottom] + field2[bottom_right]) * (field1[bottom] - field1[index])) / 4);
    return (1 / h) * (((field1[top] + field1[index]) * (field2[right] + field2[index])) / 4 - ((field1[index] + field1[bottom]) * (field2[bottom] + field2[bottom_right])) / 4) + donor_cell_correction;
  }
}

template <Grid G, Offset Direction>
inline double dxx(const G& field1, const G& field2, size_t I, double h, double alpha, const AdjMap& map)
{
  size_t Iplus = get_index<Direction, Sign::Plus>(I, map);
  size_t Iminus = get_index<Direction, Sign::Minus>(I, map);
  double donor_cell_correction = alpha * (1 / h) * ((std::abs(field1[Iplus] + field1[I]) * (field2[I] - field2[Iplus])) / 4 - (std::abs(field1[Iminus] + field1[I]) * (field2[Iminus] - field2[I])) / 4);
  return (1 / h) * (((field1[Iplus] + field1[I]) * (field2[Iplus] + field2[I])) / 4 - ((field1[Iminus] + field1[I]) * (field2[I] + field2[Iminus])) / 4) + donor_cell_correction;
}

#endif // DERIVATIVES_H_
