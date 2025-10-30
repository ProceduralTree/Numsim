#include <cstdint>
#include <pde/derivatives.h>
#include <pde/pressuresolvers.h>
#include <pde/system.h>
#include <utils/settings.h>

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
  gauss_seidel(system);
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
      system.rhs[i, j] = (1 / system.dt) * (dx(F, i - 1, j, h) + dy(G, i, j - 1, h));
    }
  }
};

void set_boundary_uv(PDESystem& system)
{
  for (uint16_t j = 1; j < system.size_y + 1; j++)
  {
    system.u[0, j] = system.boundaryLeft[0];
    system.u[system.size_x, j] = system.boundaryRight[0];
  }
  for (uint16_t i = 0; i < system.size_x + 1; i++)
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

void timestep(PDESystem& system)
{
  set_boundary_uv(system);
  set_boundary_FG(system);
  calculate_FG(system);
  calculate_rhs(system);
  solve_pressure(system);
  update_uv(system);
  //
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
  double dx = x - iy;
  double dy = y - iy;
  auto w11 = (1 - dx) * (1 - dy);
  auto w12 = (1 - dx) * dy;
  auto w21 = dx * (1 - dy);
  auto w22 = dx * dy;
  return w11 * field[ix, iy] + w12 * field[ix + 1, iy] + w21 * field[ix, iy + 1] + w22 * field[ix + 1, iy + 1];
};
