#ifndef SYSTEM_H_
#define SYSTEM_H_
#include <cassert>
#include <cmath>
#include <cstdint>
#include <grid/grid.h>
#include <utils/index.h>
#include <utils/settings.h>

struct Gridsize
{
  const double x;
  const double y;
  const double x_squared;
  const double y_squared;

  Gridsize(double x, double y)
    : x(x)
    , y(y)
    , x_squared(x * x)
    , y_squared(y * y)
  {
  }
};

struct PDESystem
{
  double residual = 0;
  const Index begin;
  const Index end;
  const double Re;
  double dt;
  Grid2D p;
  Grid2D u;
  Grid2D v;
  Grid2D F;
  Grid2D G;
  Grid2D rhs;
  const uint16_t size_x;
  const uint16_t size_y;
  const Gridsize h;
  const std::array<double, 2> boundaryBottom;
  const std::array<double, 2> boundaryTop;
  const std::array<double, 2> boundaryLeft;
  const std::array<double, 2> boundaryRight;

  PDESystem(double Re, double dt, uint16_t size_x, uint16_t size_y, double hx,
    double hy, std::array<double, 2> boundaryBottom,
    std::array<double, 2> boundaryTop, std::array<double, 2> boundaryLeft,
    std::array<double, 2> boundaryRight)
    : begin({ 2, 2 })
    , end(Index(size_x + 2, size_y + 2)) //+ Offset(1, 1))
    , Re(Re)
    , dt(dt)
    , p(Grid2D(begin, end))
    , u(Grid2D(begin, end - Ix))
    , v(Grid2D(begin, end - Iy))
    , F(Grid2D(begin, end - Ix))
    , G(Grid2D(begin, end - Iy))
    , rhs(Grid2D(begin, end))
    , size_x(size_x)
    , size_y(size_y)
    , h(Gridsize(hx, hy))
    , boundaryBottom(boundaryBottom)
    , boundaryTop(boundaryTop)
    , boundaryLeft(boundaryLeft)
    , boundaryRight(boundaryRight)
  {
  }
  PDESystem(const PDESystem&) = delete;
  PDESystem& operator=(const PDESystem&) = delete;
};

void step(PDESystem& system, uint16_t i);
void print_pde_system(const PDESystem& sys);

double interpolate_at(const PDESystem& sys, const Grid2D& field, Index I, Offset o);

constexpr std::array<double, 2> boundary(PDESystem& system, Offset offset)
{
  if (offset == Offset(1, 0))
    return system.boundaryRight;
  else if (offset == Offset(-1, 0))
    return system.boundaryLeft;
  else if (offset == Offset(0, 1))
    return system.boundaryTop;
  else
    return system.boundaryBottom;
};

#endif // SYSTEM_H_
