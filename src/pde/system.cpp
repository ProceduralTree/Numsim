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

void calculate_F(PDESystem& system, Index I)
{
  auto& u = system.u;
  auto& v = system.v;
  auto& h = system.h;
  auto& alpha = system.settings.alpha;
  system.F[I] = u[I] + system.dt * (1 / system.Re * (dd(Ix, u, I, h.x_squared) + dd(Iy, u, I, h.y_squared)) - dxx(Ix, u, u, I, h.x, alpha) - duv(Iy, u, v, I, h.y, alpha));
};

void calculate_G(PDESystem& system, Index I)
{
  auto& u = system.u;
  auto& v = system.v;
  auto& h = system.h;
  auto& alpha = system.settings.alpha;
  system.G[I] = v[I] + system.dt * (1 / system.Re * (dd(Ix, v, I, h.x_squared) + dd(Iy, v, I, h.y_squared)) - dxx(Iy, v, v, I, h.y, alpha) - duv(Ix, u, v, I, h.x, alpha));
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
  if (system.settings.pressureSolver == Settings::SOR)
  {
    int i;
    for (int iter = 0; iter < system.settings.maximumNumberOfIterations; iter++)
    {
      system.residual = 0;
      broadcast_boundary(
        [&](PDESystem& s, Index I, Offset o) { s.p[I + o] = s.p[I]; },
        system, system.p);
      broadcast(sor_step, system, system.p.range);
      i = iter;
      if (system.residual < system.settings.epsilon)
        break;
    }
    cout << "SOR Steps: " << i << endl;
  } else
  {
    int i;
    for (int iter = 0; iter < system.settings.maximumNumberOfIterations; iter++)
    {
      system.residual = 0;
      broadcast_boundary(
        [&](PDESystem& s, Index I, Offset o) { s.p[I + o] = s.p[I]; },
        system, system.p);
      broadcast(gauss_seidel_step, system, system.p.range);
      i = iter;
      if (system.residual < system.settings.epsilon)
        break;
    }
    cout << "Gauss-Seidel Steps: " << i << endl;
  }
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


void compute_dt(PDESystem& system)
{
  double umax = 0;
  double vmax = 0;
  umax = std::max(system.u.max(), (-system.u.min()));
  vmax = std::max(system.v.max(), (-system.v.min()));
  double dt1 = (system.Re / 2) * ((system.h.x_squared * system.h.y_squared) / ((system.h.x_squared) + (system.h.y_squared)));
  double dt2 = system.h.x / umax;
  double dt3 = system.h.y / vmax;
  system.dt = std::min(dt1, std::min(dt2, dt3)) * system.settings.tau;
};

void step(PDESystem& system, uint16_t i)
{

  set_uv_boundary(system);

  compute_dt(system);

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
double interpolate_at(const PDESystem& sys, const Grid2D& field, Index I, Offset o)
{
  return (field[I] + field[I - o]) / 2;
};
