#ifndef VECTOR_H_
#define VECTOR_H_

#include "linalg/matrix.h"
#include "utils/index.h"
#include <cstdint>
#include <grid/grid.h>

template <typename Operator, typename... Args>
double reduce(Operator&& O, Grid2D&& g, Args&&... args)
{
  double result = 0;

#pragma omp parallel for simd collapse(2) reduction(+ : result)
  for (uint16_t j = g.begin.y; j <= g.end.y; j++)
  {
    for (uint16_t i = g.begin.x; i <= g.end.x; i++)
    {
      result += std::forward<Operator>(O)(Index { i, j }, g, std::forward<Args>(args)...);
    }
  }
  return result;
}

double times(Index I, const Grid2D& a, const Grid2D& b)
{
  return a[I] * b[I];
};

double dot(Grid2D& a, Grid2D& b)
{
  return reduce(times, a, b);
};

double Adot(LaplaceMatrixOperator A, Grid2D& a, Grid2D& b)
{

  return reduce(
    [&](Index I, auto a, auto b, auto A) {
      return A(a, I) * b[I];
    },
    a, b, A);
}

#endif // VECTOR_H_
