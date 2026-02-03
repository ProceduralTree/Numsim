#ifndef MATRIX_H_
#define MATRIX_H_

#include "grid/adjacencymap.h"
#include "grid/boundary.h"
#include "grid/densetree.h"
#include "grid/grid.h"
#include "grid/sparsegrid.h"
#include "grid/zindex.h"
#include "pde/system.h"
#include "utils/index.h"
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <format>
#include <iostream>
struct SparseMatrixOperator
{
  const Gridsize h;
  const AdjMap& map;
  SparseMatrixOperator(const SparseMatrixOperator&) = default;
  SparseMatrixOperator(SparseMatrixOperator&&) = default;
  SparseMatrixOperator& operator=(const SparseMatrixOperator&) = delete;
  SparseMatrixOperator& operator=(SparseMatrixOperator&&) = delete;
  SparseMatrixOperator(const Gridsize& grid, const AdjMap& map)
    : h(grid)
    , map(map) { };

  constexpr double operator()(size_t index, uint16_t depth, const SparseGrid2D<double>& vec) const
  {
    size_t top = get_index<Iy, Sign::Plus>(index, map);
    size_t bottom = get_index<Iy, Sign::Minus>(index, map);
    size_t left = get_index<Ix, Sign::Minus>(index, map);
    size_t right = get_index<Ix, Sign::Plus>(index, map);

    double local_hx_2_inv = 1. / h.x_squared(depth);
    double local_hy_2_inv = 1. / h.y_squared(depth);
    double res = ((vec[left] + vec[right]) * local_hx_2_inv) + ((vec[bottom] + vec[top]) * local_hy_2_inv);
    double a_ij = -2. * (local_hx_2_inv + local_hy_2_inv);
    res += a_ij * vec[index];
    return res;
  }
  constexpr double operator[](uint16_t depth)
  {
    double local_hx_2_inv = 1. / h.x_squared(depth);
    double local_hy_2_inv = 1. / h.y_squared(depth);
    double a_ij = -2. * (local_hx_2_inv + local_hy_2_inv);
    return a_ij;
  };
};

#endif // MATRIX_H_
