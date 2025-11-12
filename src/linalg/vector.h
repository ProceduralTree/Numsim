#ifndef VECTOR_H_
#define VECTOR_H_

#include "linalg/matrix.h"
#include "utils/index.h"
#include <cstdint>
#include <grid/grid.h>

template <typename Operator, typename... Args>
inline double reduce(Operator&& O, const Grid2D& g, Args&&... args)
{
  double result = 0;

#pragma omp parallel for collapse(2) reduction(+ : result)
  for (uint16_t j = g.begin.y; j <= g.end.y; j++)
  {
    for (uint16_t i = g.begin.x; i <= g.end.x; i++)
    {
      result += std::forward<Operator>(O)(Index { i, j }, g, std::forward<Args>(args)...);
    }
  }
  return result;
}

constexpr double times(Index I, const Grid2D& a, const Grid2D& b)
{
  return a[I] * b[I];
};

constexpr double dot(Grid2D& a, Grid2D& b)
{
  return reduce(times, a, b);
};

double Adot(LaplaceMatrixOperator A, Grid2D& a, Grid2D& b)
{

  return reduce(
    [&](Index I, const Grid2D& a, const Grid2D& b, LaplaceMatrixOperator A) {
      return A(a, I) * b[I];
    },
    a, b, A);
}

void saxpy(Index I, Grid2D& result, Grid2D& y, Grid2D& x, LaplaceMatrixOperator A)
{
  result[I] = A(x, I) + y[I];
};

#endif // VECTOR_H_
