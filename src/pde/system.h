#ifndef SYSTEM_H_
#define SYSTEM_H_
#include "grid/boundary.h"
#include "grid/densetree.h"
#include "grid/sparsegrid.h"
#include <cassert>
#include <cmath>
#include <cstdint>
#include <grid/grid.h>
#include <utils/index.h>
#include <utils/settings.h>

struct CGSolver;

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
  // TreeIndices indices;
  // QuadTree<Index> local_to_global_chache;
  // QuadTree<std::array<Index, 4>> neighbours;

  double dt;
  const BoundaryFlags& boundary;
  SparseGrid2D<double> p;
  SparseGrid2D<double> u;
  SparseGrid2D<double> v;
  SparseGrid2D<double> F;
  SparseGrid2D<double> G;
  SparseGrid2D<double> rhs;
  const Gridsize h;

  PDESystem(const Settings& settings, const BoundaryFlags& flags)
    : settings(settings)
    , boundary(flags)
    , p(flags.tree)
    , u(flags.tree)
    , v(flags.tree)
    , F(flags.tree)
    , G(flags.tree)
    , rhs(flags.tree)
    , h(Gridsize(settings)) { };

  PDESystem(const PDESystem&) = delete;
  PDESystem& operator=(const PDESystem&) = delete;

  PDESystem(PDESystem&&) = default;
  PDESystem& operator=(PDESystem&&) = delete;
};

void step(PDESystem& system, CGSolver& solver, double time);
void print_pde_system(const PDESystem& sys);

double interpolate_u(const PDESystem& sys, const SparseGrid2D<double>& field, Index I);
double interpolate_v(const PDESystem& sys, const SparseGrid2D<double>& field, Index I);
double interpolate_p(const PDESystem& sys, const SparseGrid2D<double>& field, Index I);
void set_uv_boundary(PDESystem& system);

#endif // SYSTEM_H_
