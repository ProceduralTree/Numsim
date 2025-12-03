#ifndef SYSTEM_H_
#define SYSTEM_H_
#include <cassert>
#include <cmath>
#include <cstdint>
#include <grid/grid.h>
#include <utils/index.h>
#include <utils/partitioning.h>
#include <utils/settings.h>

struct Gridsize
{
  const double x;
  const double y;
  const double x_squared;
  const double y_squared;

  Gridsize(const Settings& settings)
    : x(settings.physicalSize[0] / static_cast<double>(settings.nCells[0]))
    , y(settings.physicalSize[1] / static_cast<double>(settings.nCells[1]))
    , x_squared(x * x)
    , y_squared(y * y)
  {
  }
};

struct PDESystem
{
  const Settings& settings;
  double residual = 0;
  const Index begin;
  const Index end;
  double dt;
  Grid2D p;
  Grid2D u;
  Grid2D v;
  Grid2D F;
  Grid2D G;
  Grid2D rhs;
  const Gridsize h;
  MPIInfo partitioning;

  PDESystem(const Settings& settings, const Partitioning::MPIInfo& mpiInfo)
    : settings(settings)
    , begin({ 2, 2 })
    , end(Index(mpiInfo.nCells[0] + 1, mpiInfo.nCells[1] + 1))
    , p(Grid2D(begin, end))
    , u(Grid2D((mpiInfo.left_neighbor >= 0) ? begin - Ix : begin, (mpiInfo.right_neighbor >= 0) ? end : end - Ix))
    , v(Grid2D((mpiInfo.bottom_neighbor >= 0) ? begin - Iy : begin, (mpiInfo.Top_neighbor >= 0) ? end : end - Iy))
    , F(Grid2D((mpiInfo.left_neighbor >= 0) ? begin - Ix : begin, (mpiInfo.right_neighbor >= 0) ? end : end - Ix))
    , G(Grid2D((mpiInfo.bottom_neighbor >= 0) ? begin - Iy : begin, (mpiInfo.Top_neighbor >= 0) ? end : end - Iy))
    , rhs(Grid2D(begin, end))
    , h(Gridsize(settings))
    , partitioning(mpiInfo) {
    };
  PDESystem(const PDESystem&) = delete;
  PDESystem& operator=(const PDESystem&) = delete;
};

void step(PDESystem& system, double time);
void print_pde_system(const PDESystem& sys);

double interpolate_u(const PDESystem& sys, const Grid2D& field, Index I);
double interpolate_v(const PDESystem& sys, const Grid2D& field, Index I);
double interpolate_p(const PDESystem& sys, const Grid2D& field, Index I);

#endif // SYSTEM_H_
