#ifndef MATRIX_H_
#define MATRIX_H_

#include "grid/grid.h"
#include "pde/system.h"
#include "utils/index.h"
struct SparseMatrixOperator
{
};

struct LaplaceMatrixOperator
{
  const Gridsize h;
  const double h_x_squared_inv;
  const double h_y_squared_inv;
  const double a_ij;

  LaplaceMatrixOperator(const Gridsize& grid)
    : h(grid)
    , h_x_squared_inv(1.0 / (grid.x_squared))
    , h_y_squared_inv(1.0 / (grid.y_squared))
    , a_ij(-2.0 * (1.0 / (grid.x_squared) + 1.0 / (grid.y_squared))) {
    };

  constexpr double operator()(const Grid2D& vec, Index I) const
  {

    double res = ((vec[I - Ix] + vec[I + Ix]) * h_x_squared_inv) + ((vec[I - Iy] + vec[I + Iy]) / h_y_squared_inv);
    double a_ij = -2 * (1 / h.y_squared) - 2 * (1 / h.x_squared);
    res += a_ij * vec[I];
    return res;
  }
};

#endif // MATRIX_H_
