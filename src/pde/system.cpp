#include "grid/grid.h"
#include "output/vtk.h"
#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <pde/derivatives.h>
#include <pde/pressuresolvers.h>
#include <pde/system.h>
#include <type_traits>
#include <utils/broadcast.h>
#include <utils/index.h>
#include <utils/settings.h>

#define U_RANGE system.begin - Ix, system.end
#define V_RANGE system.begin - Iy, system.end

void calculate_F(Index I, PDESystem& system)
{
  auto& u = system.u;
  auto& v = system.v;
  auto& h = system.h;
  auto& alpha = Settings::get().alpha;
  system.F[I] = u[I] + system.dt * (1 / system.settings.re * (dd(Ix, u, I, h.x_squared) + dd(Iy, u, I, h.y_squared)) - dxx(Ix, u, u, I, h.x, alpha) - duv(Iy, u, v, I, h.y, alpha));
};

void calculate_G(Index I, PDESystem& system)
{
  auto& u = system.u;
  auto& v = system.v;
  auto& h = system.h;
  auto& alpha = Settings::get().alpha;
  system.G[I] = v[I] + system.dt * (1 / system.settings.re * (dd(Ix, v, I, h.x_squared) + dd(Iy, v, I, h.y_squared)) - dxx(Iy, v, v, I, h.y, alpha) - duv(Ix, u, v, I, h.x, alpha));
};

void update_u(Index I, PDESystem& system)
{
  system.u[I] = system.F[I] - system.dt * d(Ix, system.p, I, system.h.x);
};
void update_v(Index I, PDESystem& system)
{
  system.v[I] = system.G[I] - system.dt * d(Iy, system.p, I, system.h.y);
};

void solve_pressure(PDESystem& system)
{
  if (Settings::get().pressureSolver == Settings::SOR)
  {
    auto solver = SORSolver();
    solve(solver, system);
  } else
  {
    auto solver = GaussSeidelSolver();
    solve(solver, system);
  }
};

void calculate_pressure_rhs(Index I, PDESystem& system)
{
  auto& F = system.F;
  auto& G = system.G;
  auto& h = system.h;
  system.rhs[I] = (1 / system.dt) * (d(Ix, F, I - Ix, h.x) + d(Iy, G, I - Iy, h.y));
};

void set(Index I, Offset O, Grid2D& array, double value)
{
  array[I] = value;
};

void set_with_neighbour(Index I, Offset O, Grid2D& array, double value)
{
  array[I] = 2 * value - array[I - O];
}

void set_uv_boundary(PDESystem& system)
{
  broadcast(set, system.u.boundary.left, -Ix, system.u, system.settings.dirichletBcLeft[0]);
  broadcast(set, system.u.boundary.right, Ix, system.u, system.settings.dirichletBcRight[0]);
  broadcast(set_with_neighbour, system.u.boundary.top, Iy, system.u, system.settings.dirichletBcTop[0]);
  broadcast(set_with_neighbour, system.u.boundary.bottom, -Iy, system.u, system.settings.dirichletBcBottom[0]);

  broadcast(set_with_neighbour, system.v.boundary.left, -Ix, system.v, system.settings.dirichletBcLeft[1]);
  broadcast(set_with_neighbour, system.v.boundary.right, Ix, system.v, system.settings.dirichletBcRight[1]);
  broadcast(set, system.v.boundary.top, Iy, system.v, system.settings.dirichletBcTop[1]);
  broadcast(set, system.v.boundary.bottom, -Iy, system.v, system.settings.dirichletBcBottom[1]);
};

void compute_dt(PDESystem& system)
{
  double umax = 0;
  double vmax = 0;
  umax = std::max(system.u.max(), (-system.u.min()));
  vmax = std::max(system.v.max(), (-system.v.min()));
  double dt1 = (system.settings.re / 2) * ((system.h.x_squared * system.h.y_squared) / ((system.h.x_squared) + (system.h.y_squared)));
  double dt2 = system.h.x / umax;
  double dt3 = system.h.y / vmax;
  system.dt = std::min(dt1, std::min(dt2, dt3)) * Settings::get().tau;
};

void step(PDESystem& system, uint16_t i)
{

  set_uv_boundary(system);

  compute_dt(system);

  broadcast_boundary(copy, system.F.boundary, system.u, system.F);
  broadcast_boundary(copy, system.G.boundary, system.v, system.G);

  broadcast(calculate_F, system.u.range, system);
  broadcast(calculate_G, system.v.range, system);

  broadcast(calculate_pressure_rhs, system.p.range, system);

  solve_pressure(system);

  broadcast(update_u, system.u.range, system);
  broadcast(update_v, system.v.range, system);

  set_uv_boundary(system);
};

void print_pde_system(const PDESystem& sys)
{
  printf("╔═══════════════════════════════════════════════╗\n");
  printf("║              PDE System Summary               ║\n");
  printf("╚═══════════════════════════════════════════════╝\n");
  Settings::get().printSettings();
}
double interpolate_u(const PDESystem& sys, const Grid2D& field, Index I)
{
  return (field[I] + field[I + Iy]) / 2;
};

double interpolate_v(const PDESystem& sys, const Grid2D& field, Index I)
{
  return (field[I] + field[I + Ix]) / 2;
};

double interpolate_p(const PDESystem& sys, const Grid2D& field, Index I)
{
  return (field[I] + field[I + Ix] + field[I + Iy] + field[I + Iy + Ix]) / 4;
};
