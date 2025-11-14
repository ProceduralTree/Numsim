#include "pde/pressuresolvers.h"
#include "linalg/matrix.h"
#include "linalg/vector.h"
#include "output/vtk.h"
#include "utils/broadcast.h"
#include "utils/settings.h"
#include <algorithm>
#include <cmath>
#include <grid/grid.h>
#include <grid/indexing.h>
#include <ios>
#include <pde/system.h>
#include <set>

void copy_with_offset(Index I, Offset O, Grid2D& array) { array[I] = array[I - O]; }

void gauss_seidel_step(Index I, PDESystem& system, GaussSeidelSolver& S)
{
  auto& p = system.p;
  auto& h = system.h;
  double sum_of_neighbours = ((p[I - Ix] + p[I + Ix]) / h.x_squared) + ((p[I - Iy] + p[I + Iy]) / h.y_squared);
  double a_ij = -2 * (1 / h.y_squared) - 2 * (1 / h.x_squared);
  system.residual = std::max(std::abs(sum_of_neighbours + a_ij * p[I] - system.rhs[I]), system.residual);
  p[I] = (system.rhs[I] - sum_of_neighbours) / a_ij;
}

void sor_step(Index I, PDESystem& system)
{
  auto& p = system.p;
  auto& h = system.h;
  double sum_of_neighbours = ((p[I - Ix] + p[I + Ix]) / h.x_squared) + ((p[I - Iy] + p[I + Iy]) / h.y_squared);
  double a_ij = -2 * (1 / h.y_squared) - 2 * (1 / h.x_squared);
  double residual = std::abs(sum_of_neighbours + a_ij * p[I] - system.rhs[I]);
  system.residual = std::max(residual, system.residual);
  p[I] = (1 - Settings::get().omega) * p[I] + Settings::get().omega * (system.rhs[I] - sum_of_neighbours) / a_ij;
}

void cg_iteration(PDESystem& system, CGSolver& cg)
{

  double residual_norm = INFINITY;
  double old_residual_norm = INFINITY;

  LaplaceMatrixOperator A = LaplaceMatrixOperator(system.h);

  // cg.residual = system.rhs - A*system.p;
  //  broadcast_boundary(copy_with_offset, system.p.boundary, system.p);
  broadcast_boundary(copy_with_offset, system.p.boundary, system.p);
  broadcast_boundary(copy, system.p.boundary, system.p, cg.search_direction);
  broadcast([&](Index I, CGSolver& cg, PDESystem& s, LaplaceMatrixOperator A) { cg.residual[I] = s.rhs[I] - A(s.p, I); }, system.p.range, cg, system, A);

  residual_norm = dot(cg.residual, cg.residual);

  // cg.search_direction = cg.residual;
  broadcast([&](Index I, CGSolver& cg, LaplaceMatrixOperator A) { cg.search_direction[I] = cg.residual[I]; }, system.p.range, cg, A);

  for (int iter = 0; iter < Settings::get().maximumNumberOfIterations; iter++)
  {

    old_residual_norm = residual_norm;

    broadcast_boundary(copy_with_offset, system.p.boundary, system.p);
    broadcast_boundary(copy_with_offset, cg.search_direction.boundary, cg.search_direction);
    // broadcast_boundary(copy, system.p.boundary, system.p, cg.search_direction);
    double alpha = residual_norm / Adot(A, cg.search_direction, cg.search_direction);

    // system.p = system.p + a * cg.search_direction;
    broadcast(
      [&](Index I, CGSolver& cg, PDESystem& s, double a) {
        system.p[I] = system.p[I] + a * cg.search_direction[I];
      },
      system.p.range, cg, system, alpha);

    // cg.residual = cg.residual - a * A * cg.search_direction;
    broadcast(
      [&](Index I, CGSolver& cg, PDESystem& s, double a, LaplaceMatrixOperator A) {
        cg.residual[I] = cg.residual[I] - a * A(cg.search_direction, I);
      },
      system.p.range, cg, system, alpha, A);

    // std::cout << "\rResidual:\t" << residual_norm << " Iterations:\t" << iter << std::flush;
    if (cg.residual.max() < Settings::get().epsilon)
    {
      std::cout << std::scientific << std::setprecision(14) << "Residual: " << cg.residual.max() << std::endl;
      std::cout << "converged after n=" << iter << " Iterations" << std::endl;
      break;
    }
    residual_norm = dot(cg.residual, cg.residual);
    double beta = residual_norm / old_residual_norm;

    broadcast(
      [&](Index I, CGSolver& cg, double beta) {
        cg.search_direction[I] = cg.residual[I] + beta * cg.search_direction[I];
      },
      system.p.range,
      cg, beta);
  }
}

void solve(CGSolver& cg, PDESystem& system)
{
  cg_iteration(system, cg);
}

void solve(GaussSeidelSolver& S, PDESystem& system)
{
  system.residual = 0;
  for (int iter = 0; iter < Settings::get().maximumNumberOfIterations; iter++)
  {
    system.residual = 0;
    broadcast_boundary(copy_with_offset, system.p.boundary, system.p);
    broadcast(gauss_seidel_step, system.p.range, system, S);
    if (system.residual < Settings::get().epsilon)
    {

      std::cout << std::scientific << std::setprecision(14) << "Residual: " << system.residual << std::endl;
      std::cout << "converged after n=" << iter << " Iterations" << std::endl;
      break;
    }
  }
}

void solve(SORSolver& S, PDESystem& system)
{
  for (int iter = 0; iter < Settings::get().maximumNumberOfIterations; iter++)
  {
    system.residual = 0;
    broadcast_boundary(copy_with_offset, system.p.boundary, system.p);
    broadcast(sor_step, system.p.range, system);
    if (iter % 10 && system.residual < Settings::get().epsilon)
    {

      // std::cout << std::scientific << std::setprecision(14) << "Residual: " << system.residual << std::endl;
      // std::cout << "\nSOR converged after n=" << iter << " Iterations" << std::endl;
      break;
    }
  }
}
// void solve(BlackRedSolver S, PDESystem& system)
//{
//   system.residual = 0;
//   for (int iter = 0; iter < Settings::get().maximumNumberOfIterations; iter++)
//   {
//     system.residual = 0;
//     broadcast_boundary(copy_with_offset, system.p.boundary, system.p);
//     broadcast_blackred(gauss_seidel_step, system.p.range, system, S);
//     if (S.residual.max() < Settings::get().epsilon)
//     {
//
//       std::cout << std::scientific << std::setprecision(14) << "Residual: " << system.residual << std::endl;
//       std::cout << "converged after n=" << iter << " Iterations" << std::endl;
//       break;
//     }
//   }
// }
