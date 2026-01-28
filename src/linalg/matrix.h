#ifndef MATRIX_H_
#define MATRIX_H_

#include "grid/boundary.h"
#include "grid/densetree.h"
#include "grid/grid.h"
#include "grid/sparsegrid.h"
#include "pde/system.h"
#include "utils/broadcast.h"
#include "utils/index.h"
#include <cstddef>
// struct SparseMatrixOperator
//{
//   SparseMatrixOperator(const SparseMatrixOperator&) = default;
//   SparseMatrixOperator(SparseMatrixOperator&&) = default;
//   SparseMatrixOperator& operator=(const SparseMatrixOperator&) = delete;
//   LaplaceMatrixOperator& operator=(SparseMatrixOperator&&) = delete;
//   SparseMatrixOperator(const BoundaryFlags& flags);
//
//   inline double operator()(const SparseGrid2D<double>& vec, size_t local_index) const
//   {
//     double res = ((vec[I - Ix] + vec[I + Ix]) * h_x_squared_inv) + ((vec[I - Iy] + vec[I + Iy]) * h_y_squared_inv);
//     res += a_ij * vec[I];
//     return res;
//   }
// };

struct LaplaceMatrixOperator
{
  const Gridsize h;
  const double h_x_squared_inv;
  const double h_y_squared_inv;
  const double a_ij;

  LaplaceMatrixOperator(const LaplaceMatrixOperator&) = default;
  LaplaceMatrixOperator(LaplaceMatrixOperator&&) = default;
  LaplaceMatrixOperator& operator=(const LaplaceMatrixOperator&) = delete;
  LaplaceMatrixOperator& operator=(LaplaceMatrixOperator&&) = delete;
  LaplaceMatrixOperator(const Gridsize& grid)
    : h(grid)
    , h_x_squared_inv(1.0 / (grid.x_squared))
    , h_y_squared_inv(1.0 / (grid.y_squared))
    , a_ij(-2.0 * (1.0 / (grid.x_squared) + 1.0 / (grid.y_squared))) {
    };

  inline double operator()(const Grid2D& vec, Index I) const
  {
    double res = ((vec[I - Ix] + vec[I + Ix]) * h_x_squared_inv) + ((vec[I - Iy] + vec[I + Iy]) * h_y_squared_inv);
    res += a_ij * vec[I];
    return res;
  }
};

#endif // MATRIX_H_
