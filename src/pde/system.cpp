#include "grid/boundary.h"
#include "grid/sparsegrid.h"
#include "utils/profiler.h"
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

inline void calculate_F(Index I, PDESystem& system)
{
  auto& u = system.u;
  auto& v = system.v;
  auto& h = system.h;
  auto& alpha = Settings::get().alpha;
  system.F[I] = u[I] + system.dt * (1 / system.settings.re * (dd(Ix, u, I, h.x_squared) + dd(Iy, u, I, h.y_squared)) - dxx(Ix, u, u, I, h.x, alpha) - duv(Iy, u, v, I, h.y, alpha)) + system.settings.g[0];
}

inline void calculate_G(Index I, PDESystem& system)
{
  auto& u = system.u;
  auto& v = system.v;
  auto& h = system.h;
  auto& alpha = Settings::get().alpha;
  system.G[I] = v[I] + system.dt * (1 / system.settings.re * (dd(Ix, v, I, h.x_squared) + dd(Iy, v, I, h.y_squared)) - dxx(Iy, v, v, I, h.y, alpha) - duv(Ix, u, v, I, h.x, alpha)) + +system.settings.g[1];
}

inline void update_u(Index I, PDESystem& system)
{
  system.u[I] = system.F[I] - system.dt * d(Ix, system.p, I, system.h.x);
}
void update_v(Index I, PDESystem& system)
{
  system.v[I] = system.G[I] - system.dt * d(Iy, system.p, I, system.h.y);
}

void solve_pressure(PDESystem& system, CGSolver& solver)
{
  // ProfileScope("Pressure Solver");
  // switch (Settings::get().pressureSolver) // Settings::get().pressureSolver)
  //{
  // case Settings::CG:
  //{
  //   static auto solver = CGSolver(system);
  //   solve(solver, system);
  //   break;
  // }
  // default:
  //{
  //   auto solver = SORSolver();
  //   solve(solver, system);
  //   break;
  // }
  // }
  solve(solver, system);
}

inline void calculate_pressure_rhs(Index I, PDESystem& system)
{
  auto& F = system.F;
  auto& G = system.G;
  auto& h = system.h;
  system.rhs[I] = (1 / system.dt) * (d(Ix, F, I - Ix, h.x) + d(Iy, G, I - Iy, h.y));
}

inline void set_with_neighbour(Index I, Offset O, SparseGrid2D<double>& array, double value)
{
  array[I] = 2 * value - array[I - O];
}

void set_uv_boundary(PDESystem& system)
{
  // clang-format off
  broadcast(set_with_neighbour , system.boundary , static_cast<uint16_t>(BoundaryType::U_TOP)    ,  Iy , system.u , system.settings.dirichletBcTop[0]);
  broadcast(set_with_neighbour , system.boundary , static_cast<uint16_t>(BoundaryType::U_BOTTOM) , -Iy , system.u , system.settings.dirichletBcBottom[0]);
  broadcast(set                , system.boundary , static_cast<uint16_t>(BoundaryType::U_LEFT)   , -Ix , system.u , system.settings.dirichletBcLeft[0]);
  broadcast(set                , system.boundary , static_cast<uint16_t>(BoundaryType::U_RIGHT)  ,  Ix , system.u , system.settings.dirichletBcRight[0]);
  // clang-format on
  // clang-format off
  broadcast(set_with_neighbour , system.boundary , static_cast<uint16_t>(BoundaryType::V_TOP)    ,  Iy , system.v , system.settings.dirichletBcTop[1]);
  broadcast(set_with_neighbour , system.boundary , static_cast<uint16_t>(BoundaryType::V_BOTTOM) , -Iy , system.v , system.settings.dirichletBcBottom[1]);
  broadcast(set                , system.boundary , static_cast<uint16_t>(BoundaryType::V_LEFT)   , -Ix , system.v , system.settings.dirichletBcLeft[1]);
  broadcast(set                , system.boundary , static_cast<uint16_t>(BoundaryType::V_RIGHT)  ,  Ix , system.v , system.settings.dirichletBcRight[1]);
  // clang-format on
};

void compute_dt(PDESystem& system)
{
  //   ProfileScope("Compute dt");
  //   double umax = 0;
  //   double vmax = 0;
  //   umax = std::max(system.u.max(), (-system.u.min()));
  //   vmax = std::max(system.v.max(), (-system.v.min()));
  //   double dt1 = (system.settings.re / 2) * ((system.h.x_squared * system.h.y_squared) / ((system.h.x_squared) + (system.h.y_squared)));
  //   double dt2 = system.h.x / umax;
  //   double dt3 = system.h.y / vmax;
  //   system.dt = std::min(dt1, std::min(dt2, dt3)) * system.settings.tau;
  //   system.dt = std::min(system.settings.maximumDt, system.dt);
  //   system.dt = std::max(1e-10, system.dt);
  system.dt = 1e-5;
};

void update_velocity(PDESystem& system)
{
  broadcast(update_u, system.boundary, static_cast<uint16_t>(BoundaryType::U_Inside), system);
  broadcast(update_v, system.boundary, static_cast<uint16_t>(BoundaryType::V_Inside), system);
  ////MPI_COMM_BUFFER* u_comm_buffer = new MPI_COMM_BUFFER(system.u, system.u.boundary.u_ghosts(), MPI_COMM_WORLD, system.partitioning, 16);
  ////MPI_COMM_BUFFER* v_comm_buffer = new MPI_COMM_BUFFER(system.v, system.v.boundary.v_ghosts(), MPI_COMM_WORLD, system.partitioning, 32);
  // delete u_comm_buffer;
  // delete v_comm_buffer;
}

void step(PDESystem& system, CGSolver solver, double time)
{
  ProfileScope("Time Step");

  set_uv_boundary(system);

  compute_dt(system);
  //
  broadcast(copy, system.boundary, static_cast<uint16_t>(BoundaryType::U_BOUNDARY), Offset { 0, 0 }, system.u, system.F);
  broadcast(copy, system.boundary, static_cast<uint16_t>(BoundaryType::U_BOUNDARY), Offset { 0, 0 }, system.v, system.G);

  broadcast(calculate_F, system.boundary, static_cast<uint16_t>(BoundaryType::U_Inside), system);
  broadcast(calculate_G, system.boundary, static_cast<uint16_t>(BoundaryType::V_Inside), system);

  broadcast(calculate_pressure_rhs, system.boundary, static_cast<uint16_t>(BoundaryType::P_Inside), system);

  solve_pressure(system, solver);

  update_velocity(system);

  set_uv_boundary(system);
}

void print_pde_system(const PDESystem& sys)
{
  printf("╔═══════════════════════════════════════════════╗\n");
  printf("║              PDE System Summary               ║\n");
  printf("╚═══════════════════════════════════════════════╝\n");
  Settings::get().printSettings();
}
double interpolate_u(const PDESystem& sys, const SparseGrid2D<double>& field, Index I)
{
  return (field[I] + field[I + Iy]) / 2;
}

double interpolate_v(const PDESystem& sys, const SparseGrid2D<double>& field, Index I)
{
  return (field[I] + field[I + Ix]) / 2;
}

double interpolate_p(const PDESystem& sys, const SparseGrid2D<double>& field, Index I)
{
  return (field[I] + field[I + Ix] + field[I + Iy] + field[I + Iy + Ix]) / 4;
}
