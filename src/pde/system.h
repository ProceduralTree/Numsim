#ifndef SYSTEM_H_
#define SYSTEM_H_
#include "grid.h"
#include "indexing.h"
#include <array>
#include <cstdint>

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
    : Re(Re)
    , dt(dt)
    , size_x(size_x)
    , size_y(size_y)
    , p(Grid2D(size_x + 2, size_y + 2))
    , u(Grid2D(size_x + 2, size_y + 2))
    , v(Grid2D(size_x + 2, size_y + 2))
    , F(Grid2D(size_x + 2, size_y + 2))
    , G(Grid2D(size_x + 2, size_y + 2))
    , rhs(Grid2D(size_x + 2, size_y + 2))
    , h(Gridsize(hx, hy))
    , boundaryBottom(boundaryBottom)
    , boundaryTop(boundaryTop)
    , boundaryLeft(boundaryLeft)
    , boundaryRight(boundaryRight)
  {
  }
};

void timestep(PDESystem system);
void print_pde_system(const PDESystem& sys);

#endif // SYSTEM_H_
