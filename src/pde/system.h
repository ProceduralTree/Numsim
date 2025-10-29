#ifndef SYSTEM_H_
#define SYSTEM_H_
#include <cassert>
#include <cmath>
#include <cstdint>
#include <grid/grid.h>
#include <grid/indexing.h>

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
  Grid2D b;
  Grid2D tmp;
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
    , p(Grid2D(size_x, size_y))
    , u(Grid2D(size_x, size_y))
    , v(Grid2D(size_x, size_y))
    , F(Grid2D(size_x, size_y))
    , G(Grid2D(size_x, size_y))
    , b(Grid2D(size_x, size_y))
    , tmp(Grid2D(size_x, size_y))
    , size_x(size_x)
    , size_y(size_y)
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

double interpolate_at(const PDESystem& sys, const Grid2D& field, double x, double y);

#endif // SYSTEM_H_
