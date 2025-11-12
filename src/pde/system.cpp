#include "grid/grid.h"
#include "output/vtk.h"
#include <cstdint>
#include <pde/derivatives.h>
#include <pde/pressuresolvers.h>
#include <pde/system.h>
#include <type_traits>
#include <utils/broadcast.h>
#include <utils/index.h>
#include <utils/settings.h>

#define U_RANGE system.begin - Ix, system.end
#define V_RANGE system.begin - Iy, system.end

void calculate_F(PDESystem& system, Index I)
{
  auto& u = system.u;
  auto& v = system.v;
  auto& h = system.h;
  system.F[I] = u[I] + system.dt * (1 / system.Re * (dd(Ix, u, I, h.x_squared) + dd(Iy, u, I, h.y_squared)) - duv(Ix, u, u, I, h.x) - duv(Iy, u, v, I, h.y));
};
void calculate_G(PDESystem& system, Index I)
{
  auto& u = system.u;
  auto& v = system.v;
  auto& h = system.h;
  system.G[I] = v[I] + system.dt * (1 / system.Re * (dd(Ix, v, I, h.x_squared) + dd(Iy, v, I, h.y_squared)) - duv(Iy, v, v, I, h.y) - duv(Ix, u, v, I, h.x));
};

void update_u(PDESystem& system, Index index)
{
  system.u[index] = system.F[index] - system.dt * d(Ix, system.p, index, system.h.x);
};
void update_v(PDESystem& system, Index index)
{
  system.v[index] = system.G[index] - system.dt * d(Iy, system.p, index, system.h.y);
};

void solve_pressure(PDESystem& system)
{
  CGSolver solver = CGSolver(system);
  // SORSolver solver = SORSolver();
  //        GaussSeidelSolver solver = GaussSeidelSolver();
  solve(solver, system);
};

void calculate_pressure_rhs(PDESystem& system, Index I)
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
  broadcast(set, system.u.boundary.left, -Ix, system.u, system.boundaryLeft[0]);
  broadcast(set, system.u.boundary.right, Ix, system.u, system.boundaryRight[0]);
  broadcast(set_with_neighbour, system.u.boundary.top, Iy, system.u, system.boundaryTop[0]);
  broadcast(set_with_neighbour, system.u.boundary.bottom, -Iy, system.u, system.boundaryBottom[0]);

  broadcast(set_with_neighbour, system.v.boundary.left, -Ix, system.v, system.boundaryLeft[1]);
  broadcast(set_with_neighbour, system.v.boundary.right, Ix, system.v, system.boundaryRight[1]);
  broadcast(set, system.v.boundary.top, Iy, system.v, system.boundaryTop[1]);
  broadcast(set, system.v.boundary.bottom, -Iy, system.v, system.boundaryBottom[1]);
};

void step(PDESystem& system, uint16_t i)
{

  set_uv_boundary(system);

  broadcast_boundary(copy, system.F.boundary, system.u, system.F);
  broadcast_boundary(copy, system.G.boundary, system.v, system.G);

  broadcast(calculate_F, system, system.u.range);
  broadcast(calculate_G, system, system.v.range);

  broadcast(calculate_pressure_rhs, system, system.p.range);

  solve_pressure(system);

  broadcast(update_u, system, system.u.range);
  broadcast(update_v, system, system.v.range);

  set_uv_boundary(system);
};

void print_pde_system(const PDESystem& sys)
{
  printf("╔═══════════════════════════════════════════════╗\n");
  printf("║              PDE System Summary               ║\n");
  printf("╚═══════════════════════════════════════════════╝\n");
  sys.settings.printSettings();
}
double interpolate_at(const PDESystem& sys, const Grid2D& field, double x, double y)
{
  assert(x >= 0);
  assert(y >= 0);
  uint16_t ix = std::floor(x);
  uint16_t iy = std::floor(y);
  double dx = x - ix;
  double dy = y - iy;
  auto w11 = (1 - dx) * (1 - dy);
  auto w12 = (1 - dx) * dy;
  auto w21 = dx * (1 - dy);
  auto w22 = dx * dy;
  return w11 * field[ix, iy] + w12 * field[ix + 1, iy] + w21 * field[ix, iy + 1] + w22 * field[ix + 1, iy + 1];
};
