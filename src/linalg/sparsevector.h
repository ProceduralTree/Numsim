#ifndef SPARSEVECTOR_H_
#define SPARSEVECTOR_H_

#include "grid/boundary.h"
#include "grid/sparsegrid.h"
#include "linalg/matrix.h"
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>

namespace SparseVector {

inline void axpy(size_t index, uint16_t depth, SparseGrid2D<double>& result, double a, const SparseGrid2D<double>& x, const SparseGrid2D<double>& y)
{
  result[index] = a * x[index] + y[index];
};
inline void aAxpy(size_t index, uint16_t depth, SparseGrid2D<double>& result, double a, SparseMatrixOperator A, const SparseGrid2D<double>& x, const SparseGrid2D<double>& y)
{
  result[index] = a * A(index, depth, x) + y[index];
};
inline void precondition(size_t index, uint16_t depth, SparseGrid2D<double>& result, SparseMatrixOperator A, const SparseGrid2D<double>& x)
{
  result[index] = (1. / A[depth]) * x[index];
};
inline void axpy_with_jacoby_precondition(size_t index, uint16_t depth, SparseGrid2D<double>& result, double a, SparseMatrixOperator A, const SparseGrid2D<double>& x, const SparseGrid2D<double>& y)
{
  result[index] = (1. / A[depth]) * a * x[index] + y[index];
};
inline void times(size_t index, uint16_t depth, const SparseGrid2D<double>& a, const SparseGrid2D<double>& b, double& result)
{
  result += a[index] * b[index];
}
inline void Atimes(size_t index, uint16_t depth, const SparseGrid2D<double>& a, SparseMatrixOperator A, const SparseGrid2D<double>& b, double& result)
{
  result += a[index] * A(index, depth, b);
}

constexpr double dot(SparseGrid2D<double>& a, SparseGrid2D<double>& b, const BoundaryFlags& flags)
{
  ProfileScope("sparse dot Product");
  size_t p_with_boundary = static_cast<uint16_t>(BoundaryType::P_Inside);
  double result = 0;

  broadcast_cell_type(times, 0, flags.tree.maxDepth, flags, p_with_boundary, a, b, result);
  return result;
  // return distributed_sum(times, a.range, a, b);
};
constexpr double dot(SparseGrid2D<double>& a, SparseMatrixOperator A, SparseGrid2D<double>& b, const BoundaryFlags& flags)
{
  ProfileScope("sparse dot Product");
  size_t p_with_boundary = static_cast<uint16_t>(BoundaryType::P_Inside);
  double result = 0;

  broadcast_cell_type(Atimes, 0, flags.tree.maxDepth, flags, p_with_boundary, a, A, b, result);
  return result;
};
constexpr void max_(size_t index, uint16_t depth, SparseGrid2D<double>& v, double& result)
{
  result = std::max(result, v[index]);
};
constexpr void min_(size_t index, uint16_t depth, SparseGrid2D<double>& v, double& result)
{
  result = std::min(result, v[index]);
};

constexpr double max(SparseGrid2D<double>& v, const BoundaryFlags& flags)
{

  size_t p_with_boundary = static_cast<uint16_t>(BoundaryType::P);
  double result = -INFINITY;
  broadcast_cell_type(max_, 0, flags.tree.maxDepth, flags, p_with_boundary, v, result);
  return result;
};
constexpr double min(SparseGrid2D<double>& v, const BoundaryFlags& flags)
{

  size_t p_with_boundary = static_cast<uint16_t>(BoundaryType::P);
  double result = INFINITY;
  broadcast_cell_type(min_, 0, flags.tree.maxDepth, flags, p_with_boundary, v, result);
  return result;
};

};

#endif // SPARSEVECTOR_H_
