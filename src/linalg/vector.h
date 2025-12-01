#ifndef VECTOR_H_
#define VECTOR_H_

#include "linalg/matrix.h"
#include "utils/index.h"
#include "utils/profiler.h"
#include <cstdint>
#include <grid/grid.h>
#include <mpi.h>
#include <utility>

template <typename Operator, typename... Args>
inline double sum(Operator&& O, Range r, Args&&... args)
{
  double result = 0;
  ProfileScope("Reduction");
#pragma omp parallel for simd collapse(2) reduction(+ : result)
  for (uint16_t j = r.begin.y; j <= r.end.y; j++)
  {
    for (uint16_t i = r.begin.x; i <= r.end.x; i++)
    {
      result += std::forward<Operator>(O)(Index { i, j }, std::forward<Args>(args)...);
    }
  }
  return result;
}
template <typename Operator, typename... Args>
inline double distributed_sum(Operator&& O, Range r, Args&&... args)
{
  double local_sum = sum(std::forward<Operator>(O), r, std::forward<Args>(args)...);
  double global_sum = 0.;
  MPI_Allreduce(&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  return global_sum;
}

inline double times(Index I, const Grid2D& a, const Grid2D& b)
{
  return a[I] * b[I];
}

inline double dot(Grid2D& a, Grid2D& b)
{
  return distributed_sum(times, a.range, a, b);
};

inline double Axy(Index I, LaplaceMatrixOperator A, const Grid2D& x, const Grid2D& y)
{
  return A(x, I) * y[I];
}

inline double Adot(LaplaceMatrixOperator A, Grid2D& a, Grid2D& b)
{

  return distributed_sum(Axy, a.range, A, a, b);
}

inline void axpy(Index I, Grid2D& result, double a, const Grid2D& x, const Grid2D& y)
{
  result[I] = a * x[I] + y[I];
};
inline void aAxpy(Index I, Grid2D& result, double a, LaplaceMatrixOperator A, const Grid2D& x, const Grid2D& y)
{
  result[I] = a * A(x, I) + y[I];
};

#endif // VECTOR_H_
