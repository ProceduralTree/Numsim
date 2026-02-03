#ifndef SYSTEM_H_
#define SYSTEM_H_
#include "grid/adjacencymap.h"
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
struct Jacoby;

struct Gridsize
{
  const double _x;
  const double _y;

  Gridsize(const Settings& settings)
    : _x(settings.physicalSize[0] / static_cast<double>(settings.nCells[0]))
    , _y(settings.physicalSize[1] / static_cast<double>(settings.nCells[1])) { };
  constexpr double x(uint16_t depth) const
  {

    double local_cell_size = static_cast<double>(1ULL << depth);
    return _x * static_cast<double>(local_cell_size);
  };
  constexpr double y(uint16_t depth) const
  {
    double local_cell_size = static_cast<double>(1ULL << depth);
    return _y * static_cast<double>(local_cell_size);
  };
  constexpr double x_squared(uint16_t depth) const
  {
    return x(depth) * x(depth);
  };
  constexpr double y_squared(uint16_t depth) const
  {
    return y(depth) * y(depth);
  };
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
  const AdjMap adjacency_map;
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
    , adjacency_map(flags)
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

void update_velocity(PDESystem& system);
void compute_dt(PDESystem& system);
void set_uv_boundary(PDESystem& system);

void print_pde_system(const PDESystem& sys);

double interpolate_u(const PDESystem& sys, const SparseGrid2D<double>& field, Index I);
double interpolate_v(const PDESystem& sys, const SparseGrid2D<double>& field, Index I);
double interpolate_p(const PDESystem& sys, const SparseGrid2D<double>& field, Index I);
void set_uv_boundary(PDESystem& system);

template <typename Solver>
void step(PDESystem& system, Solver& solver, double time);

#endif // SYSTEM_H_
