#ifndef PRESSURESOLVERS_H_
#define PRESSURESOLVERS_H_

#include "grid/grid.h"
#include "linalg/matrix.h"
#include "linalg/vector.h"
#include <pde/system.h>
#include <utils/index.h>

struct CGSolver
{
  Grid2D residual;
  Grid2D search_direction;
  CGSolver(PDESystem& system)
    : residual(system.begin, system.end)
    , search_direction(system.begin, system.end) {
    };
};

struct GaussSeidelSolver
{
  // Grid2D residual;
  // GaussSeidelSolver(PDESystem& system)
  //   : residual(system.begin, system.end) { };
};

struct SORSolver
{
  // Grid2D residual;
  // SORSolver(PDESystem& system)
  //   : residual(system.begin, system.end) { };
};
struct BlackRedSolver
{
  Grid2D residual;
  BlackRedSolver(PDESystem& system)
    : residual(system.p.begin, system.p.end) { };
};
struct Jacoby
{
  Grid2D residual;
  Grid2D tmp;
  Jacoby(PDESystem& system)
    : residual(system.p.begin, system.p.end)
    , tmp(system.p.begin, system.p.end) { };
};

void solve(GaussSeidelSolver& S, PDESystem& system);
void solve(SORSolver& S, PDESystem& system);
void solve(CGSolver& S, PDESystem& system);
void solve(BlackRedSolver& S, PDESystem& system);
void solve(Jacoby& S, PDESystem& system);

constexpr void copy_with_offset(Index I, Offset O, Grid2D& array) { array[I] = array[I + O]; };

constexpr std::pair<double, double> jacoby_update(Index I, const PDESystem& system)
{
  auto& p = system.p;
  auto& h = system.h;
  double sum_of_neighbours = ((p[I - Ix] + p[I + Ix]) / h.x_squared) + ((p[I - Iy] + p[I + Iy]) / h.y_squared);
  double a_ij = -2. * (1. / h.y_squared) - 2. * (1. / h.x_squared);
  double residual = std::abs(sum_of_neighbours + a_ij * p[I] - system.rhs[I]);
  double update = (system.rhs[I] - sum_of_neighbours) / a_ij;

  return { update, residual };
}

inline void gauss_seidel_step(Index I, PDESystem& system, GaussSeidelSolver& S)
{
  auto& p = system.p;
  auto& h = system.h;
  double sum_of_neighbours = ((p[I - Ix] + p[I + Ix]) / h.x_squared) + ((p[I - Iy] + p[I + Iy]) / h.y_squared);
  double a_ij = -2 * (1 / h.y_squared) - 2 * (1 / h.x_squared);
  system.residual = std::max(std::abs(sum_of_neighbours + a_ij * p[I] - system.rhs[I]), system.residual);
  p[I] = (system.rhs[I] - sum_of_neighbours) / a_ij;
};

inline void sor_step(Index I, PDESystem& system)
{
  auto& p = system.p;
  auto& h = system.h;
  double sum_of_neighbours = ((p[I - Ix] + p[I + Ix]) / h.x_squared) + ((p[I - Iy] + p[I + Iy]) / h.y_squared);
  double a_ij = -2 * (1 / h.y_squared) - 2 * (1 / h.x_squared);
  double residual = std::abs(sum_of_neighbours + a_ij * p[I] - system.rhs[I]);
  system.residual = std::max(residual, system.residual);
  p[I] = (1 - Settings::get().omega) * p[I] + Settings::get().omega * (system.rhs[I] - sum_of_neighbours) / a_ij;
};
inline void black_red_step(Index I, PDESystem& system, BlackRedSolver& solver)
{
  auto [up, res] = jacoby_update(I, system);
  // DebugF("Update {} , residual {}", up, res);
  solver.residual[I]
    = res;
  system.p[I] = up;
};
constexpr void jacoby_step(Index I, PDESystem& system, Jacoby& solver)
{
  auto [up, res] = jacoby_update(I, system);
  // DebugF("Update {} , residual {}", up, res);
  solver.residual[I]
    = res;
  solver.tmp[I] = up;
};

#endif // PRESSURESOLVERS_H_
