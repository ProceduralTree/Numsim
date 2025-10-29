#ifndef DERIVATIVES_H_
#define DERIVATIVES_H_

#include <grid/grid.h>
#include <grid/indexing.h>
#include <pde/system.h>

constexpr double dx(const Grid2D& field, uint32_t i, uint32_t j,
  const Gridsize& h)
{
  return 1 / h.x * (field[i + 1, j] - field[i, j]);
};
constexpr double dy(const Grid2D& field, uint32_t i, uint32_t j,
  const Gridsize& h)
{
  return 1 / h.y * (field[i + 1, j] - field[i, j]);
};

constexpr double dx_interpolated(const Grid2D& field1, const Grid2D& field2,
  uint16_t i, uint16_t j, const Gridsize& h)
{
  return 1 / h.x * ((field1[i + 1, j] * field2[i + 1, j] + field1[i, j] * field2[i, j]) / 2 - (field1[i, j] * field2[i, j] + field1[i - 1, j] * field2[i - 1, j]) / 2);
};
constexpr double dy_interpolated(const Grid2D& field1, const Grid2D& field2,
  uint16_t i, uint16_t j, const Gridsize& h)
{
  return 1 / h.y * ((field1[i, j + 1] * field2[i, j + 1] + field1[i, j] * field2[i, j]) / 2 - (field1[i, j] * field2[i, j] + field1[i, j - 1] * field2[i, j - 1]) / 2);
};

constexpr double ddx(const Grid2D& field, uint16_t i, uint16_t j,
  const Gridsize& h)
{
  return 1 / h.x_squared * (field[i + 1, j] + field[i - 1, j] - 2 * field[i, j]);
};

constexpr double ddy(const Grid2D& field, uint16_t i, uint16_t j,
  const Gridsize& h)
{
  return 1 / h.y_squared * (field[i, j + 1] + field[i, j - 1] - 2 * field[i, j]);
};
#endif // DERIVATIVES_H_
