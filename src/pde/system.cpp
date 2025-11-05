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
  system.F[I] = u[I] + system.dt * (1 / system.Re * (dd(Ix, u, I, h.x_squared) + dd(Iy, u, I, h.y_squared)) - duv(Ix, u, u, I, h.x) - duv(Ix, u, v, I, h.x));
};
void calculate_G(PDESystem& system, Index I)
{
  auto& u = system.u;
  auto& v = system.v;
  auto& h = system.h;
  system.G[I] = v[I] + system.dt * (1 / system.Re * (dd(Ix, v, I, h.x_squared) + dd(Ix, v, I, h.y_squared)) - duv(Iy, v, v, I, h.y) - duv(Iy, u, v, I, h.y));
};

void calculate_FG(PDESystem& system)
{
  auto& u = system.u;
  auto& v = system.v;
  auto& h = system.h;
  for (uint16_t i = 1; i < system.size_x; i++)
  {
    for (uint16_t j = 1; j < system.size_x + 1; j++)
    {
      system.F[i, j] = u[i, j] + system.dt * (1 / system.Re * (ddx(u, i, j, h) + ddy(u, i, j, h)) - dx_interpolated(u, u, i, j, h) - dx_interpolated(u, v, i, j, h));
    }
  }
  for (uint16_t i = 1; i < system.size_x + 1; i++)
  {
    for (uint16_t j = 1; j < system.size_x; j++)
    {
      system.G[i, j] = v[i, j] + system.dt * (1 / system.Re * (ddy(v, i, j, h) + ddx(v, i, j, h)) - dy_interpolated(v, v, i, j, h) - dy_interpolated(u, v, i, j, h));
    }
  }
}

void update_u(PDESystem& system, Index index)
{
  system.u[index] = system.F[index] - system.dt * d(Ix, system.p, index, system.h.x);
};
void update_v(PDESystem& system, Index index)
{
  system.v[index] = system.G[index] - system.dt * d(Iy, system.p, index, system.h.y);
};

void update_uv(PDESystem& system)
{
  auto& p = system.p;
  auto& h = system.h;
  for (uint16_t i = 1; i < system.size_x; i++)
  {
    for (uint16_t j = 1; j < system.size_x + 1; j++)
    {
      system.u[i, j] = system.F[i, j] - system.dt * dx(p, i, j, h);
    }
  }
  for (uint16_t i = 1; i < system.size_x + 1; i++)
  {
    for (uint16_t j = 1; j < system.size_x; j++)
    {
      system.v[i, j] = system.G[i, j] - system.dt * dy(p, i, j, h);
    }
  }
};

void solve_pressure(PDESystem& system)
{
  system.residual = 0;
  std::cout << std::scientific << std::endl;
  std::cout << std::endl;
  std::cout << "Max Pressure: \t" << system.p.max() << "\n";
  std::cout << "Min Pressure: \t" << system.p.min() << "\n";
  std::cout << "Max Velocity x: \t" << system.u.max() << "\n";
  std::cout << "Min Velocity x: \t" << system.u.min() << "\n";
  std::cout << "Max Velocity y: \t" << system.v.max() << "\n";
  std::cout << "Min Velocity y: \t" << system.v.min() << "\n";
  std::cout << "RHS: \t" << system.rhs << "\n";
  std::cout << "G: \t" << system.F << "\n";
  std::cout << "G: \t" << system.G << "\n";
  std::cout << "v: \t" << system.v << "\n";
  std::cout << "u: \t" << system.u << "\n";
  std::cout << "p: \t" << system.p << "\n";
  std::cout << std::endl;
  for (int iter = 0; iter < 10000; iter++)
  {
    system.residual = 0;
    broadcast_boundary(
      [&](PDESystem& s, Index I, Offset o) { s.p[I + o] = s.p[I]; },
      system, system.p);
    broadcast(gauss_seidel_step, system, system.p.range);
    std::cout << "Residual : \t" << system.residual << "\t\r"
              << std::flush;
    // if (system.residual < 1e-4)
    //   break;
  }
  // gauss_seidel(system);
};

void calculate_pressure_rhs(PDESystem& system, Index I)
{
  auto& F = system.F;
  auto& G = system.G;
  auto& h = system.h;
  system.rhs[I] = (1 / system.dt) * (d(Ix, F, I - Ix, h.x) + d(Iy, G, I - Iy, h.y));
};

void calculate_rhs(PDESystem& system)
{
  auto& F = system.F;
  auto& G = system.G;
  auto& h = system.h;
  for (uint16_t i = 1; i < system.size_x + 1; i++)
  {
    for (uint16_t j = 1; j < system.size_y + 1; j++)
    {
      system.rhs[i, j] = 1 / system.dt * (dx(F, i - 1, j, h) + dy(G, i, j - 1, h));
    }
  }
};

void set_boundary_uv(PDESystem& system)
{
  for (uint16_t j = system.u.begin.y; j < system.u.end.y; j++)
  {
    system.u[system.u.begin.y - 1, j] = system.boundaryLeft[0];
    system.u[system.size_x, j] = system.boundaryRight[0];
  }
  for (uint16_t i = system.u.begin.x; i < system.u.end.x; i++)
  {
    system.u[i, 0] = 2 * system.boundaryBottom[0] - system.u[i, 1];
    system.u[i, system.size_y + 1] = 2 * system.boundaryTop[0] - system.u[i, system.size_y];
  }
  for (uint16_t i = 1; i < system.size_x + 1; i++)
  {
    system.v[i, 0] = system.boundaryBottom[1];
    system.v[i, system.size_y] = system.boundaryTop[1];
  }
  for (uint16_t j = 0; j < system.size_y + 1; j++)
  {
    system.v[0, j] = 2 * system.boundaryLeft[1] - system.v[1, j];
    system.v[system.size_x + 1, j] = 2 * system.boundaryRight[1] - system.v[system.size_x, j];
  }
};

void set_boundary_FG(PDESystem& system)
{
  for (uint16_t i = 0; i < system.size_x + 2; i++)
  {
    system.F[i, 0] = system.u[i, 0];
    system.F[i, system.size_y + 1] = system.u[i, system.size_y + 1];
    system.G[i, 0] = system.v[i, 0];
    system.G[i, system.size_y] = system.v[i, system.size_y];
  }
  for (uint16_t j = 0; j < system.size_y + 1; j++)
  {
    system.F[0, j] = system.u[0, j];
    system.F[system.size_x, j] = system.u[system.size_x, j];
    system.G[0, j] = system.v[0, j];
    system.G[system.size_y + 1, j] = system.v[system.size_x + 1, j];
  }
};

auto copy_boundary(const Grid2D& from, Grid2D& to)
{

  return [&](PDESystem& s, Index I, Offset o) { to[I + o] = from[I + o]; };
}

void step(PDESystem& system)
{

  broadcast_x_boundary(
    [&](PDESystem& s, Index I, Offset o) { s.u[I + o] = boundary(s, o)[0]; },
    system, system.u);
  broadcast_x_boundary(
    [&](PDESystem& s, Index I, Offset o) {
      s.v[I + o] = 2 * boundary(s, o)[1] - s.v[I];
    },
    system, system.v);
  broadcast_y_boundary(
    [&](PDESystem& s, Index I, Offset o) { s.u[I + o] = 2 * boundary(s, o)[0] - s.u[I]; },
    system, system.u);
  broadcast_y_boundary(
    [&](PDESystem& s, Index I, Offset o) { s.v[I + o] = boundary(s, o)[1]; },
    system, system.v);

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

  broadcast(calculate_F, system, system.u.range);
  broadcast(calculate_G, system, system.v.range);

  broadcast(calculate_pressure_rhs, system, system.p.range);

  solve_pressure(system);

  // u update
  broadcast(update_u, system, system.u.range);
  // v update
  broadcast(update_v, system, system.v.range);
};

void timestep(PDESystem& system)
{
  set_boundary_uv(system);
  set_boundary_FG(system);
  calculate_FG(system);
  calculate_rhs(system);
  solve_pressure(system);
  update_uv(system);
}
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
