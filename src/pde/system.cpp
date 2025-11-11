#include "grid/grid.h"
#include "output/vtk.h"
#include <cstdint>
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
  system.F[I] = u[I] + system.dt * (1 / system.Re * (dd(Ix, u, I, h.x_squared) + dd(Iy, u, I, h.y_squared)) - dxx(Ix, u, u, I, h.x) - duv(Iy, u, v, I, h.y));
};

void calculate_G(PDESystem& system, Index I)
{
  auto& u = system.u;
  auto& v = system.v;
  auto& h = system.h;
  system.G[I] = v[I] + system.dt * (1 / system.Re * (dd(Ix, v, I, h.x_squared) + dd(Iy, v, I, h.y_squared)) - dxx(Iy, v, v, I, h.y) - duv(Ix, u, v, I, h.x));
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

void set_boundary_uv(PDESystem& system)
{
  broadcast_x_boundary(
    [&](PDESystem& s, Index I, Offset o) { s.u[I + o] = 2 * boundary(s, o)[0] - s.u[I]; },
    system, system.u);
  broadcast_x_boundary(
    [&](PDESystem& s, Index I, Offset o) {
      s.v[I + o] = boundary(s, o)[1];
    },
    system, system.v);
  broadcast_y_boundary(
    [&](PDESystem& s, Index I, Offset o) { s.u[I + o] = boundary(s, o)[0]; },
    system, system.u);
  broadcast_y_boundary(
    [&](PDESystem& s, Index I, Offset o) { s.v[I + o] = 2 * boundary(s, o)[1] - s.v[I]; },
    system, system.v);
};

void set_boundary_FG(PDESystem& system)
{
  broadcast_x_boundary(
    [&](PDESystem& s, Index I, Offset o) { s.F[I + o] = s.u[I + o]; },
    system, system.u);
  broadcast_x_boundary(
    [&](PDESystem& s, Index I, Offset o) { s.G[I + o] = s.v[I + o]; },
    system, system.v);
  broadcast_y_boundary(
    [&](PDESystem& s, Index I, Offset o) { s.F[I + o] = s.u[I + o]; },
    system, system.u);
  broadcast_y_boundary(
    [&](PDESystem& s, Index I, Offset o) { s.G[I + o] = s.v[I + o]; },
    system, system.v);
};

void compute_dt(PDESystem& system)
{
  double umax = 0;
  double vmax = 0;
  for (uint16_t j = system.begin.y; j < system.end.y; j++)
  {
    for (uint16_t i = system.begin.x; i < system.end.x; i++)
    {
      umax = std::max(std::abs(system.u[i, j]), umax);
      vmax = std::max(std::abs(system.v[i, j]), vmax);
    }
  }
  double dt1 = (system.Re / 2) * ((system.h.x_squared * system.h.y_squared) / ((system.h.x_squared) + (system.h.y_squared)));
  double dt2 = system.h.x / umax;
  double dt3 = system.h.y / vmax;
  system.dt = std::min(dt1, std::min(dt2, dt3)) * system.settings.tau;
};

auto copy_boundary(const Grid2D& from, Grid2D& to)
{

  return [&](PDESystem& s, Index I, Offset o) { to[I + o] = from[I + o]; };
}

void step(PDESystem& system, uint16_t i)
{

  set_boundary_uv(system);
  compute_dt(system);
  set_boundary_FG(system);

  broadcast(calculate_F, system, system.u.range);
  broadcast(calculate_G, system, system.v.range);

  broadcast(calculate_pressure_rhs, system, system.p.range);

  solve_pressure(system);

  // u update
  broadcast(update_u, system, system.u.range);
  // v update
  broadcast(update_v, system, system.v.range);

  write_vtk(system, static_cast<double>(i));
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
