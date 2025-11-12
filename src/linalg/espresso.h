#ifndef ESPRESSO_H_
#define ESPRESSO_H_

#include "linalg/matrix.h"
#include <utils/index.h>

// ChatGTP supported gcc abuse:

template <typename L, typename R, typename Op>
struct Expr
{
  const L& l;
  const R& r;
  Op op;
  Expr(const L& l_, const R& r_, Op op_)
    : l(l_)
    , r(r_)
    , op(op_)
  {
  }

  double operator[](Index I) const { return op(l[I], r[I]); }
};
template <typename L, typename R>
auto operator+(const L& l, const R& r)
{
  return Expr<L, R, std::plus<>> { l, r, std::plus<> {} };
}

template <typename L, typename R>
auto operator*(const L& l, const R& r)
{
  return Expr<L, R, std::multiplies<>> { l, r, std::multiplies<> {} };
}

// Scalars
template <typename L>
auto operator*(const L& l, double scalar)
{
  return Expr<L, double, std::multiplies<>> { l, scalar, std::multiplies<> {} };
}
template <typename L>
auto operator*(const L& l, LaplaceMatrixOperator A)
{
  return Expr<L, double, std::multiplies<>> { l, A, [](L l, LaplaceMatrixOperator A) { return [&](Index I) { A(l, I); }; } };
}

template <typename R>
auto operator*(double scalar, const R& r)
{
  return r * scalar;
}

#endif // ESPRESSO_H_
