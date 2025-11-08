#include "pde/pressuresolvers.h"
#include "linalg/matrix.h"
#include "linalg/vector.h"
#include "output/vtk.h"
#include "utils/broadcast.h"
#include "utils/settings.h"
#include <cmath>
#include <grid/grid.h>
#include <grid/indexing.h>
#include <ios>
#include <pde/system.h>

void gauss_seidel_step(PDESystem& system, Index I)
{
  auto& p = system.p;
  auto& h = system.h;
  double sum_of_neighbours = ((p[I - Ix] + p[I + Ix]) / h.x_squared) + ((p[I - Iy] + p[I + Iy]) / h.y_squared);
  double a_ij = -2 * (1 / h.y_squared) - 2 * (1 / h.x_squared);
  double residual = std::abs(sum_of_neighbours + a_ij * p[I] - system.rhs[I]);
  system.residual = std::max(residual, system.residual);
  p[I] = (system.rhs[I] - sum_of_neighbours) / a_ij;
};

void cg_iteration(PDESystem& system, CGSolver& cg)
{

  double residual_norm = INFINITY;
  double old_residual_norm = INFINITY;

  LaplaceMatrixOperator A = LaplaceMatrixOperator(system.h);

  // cg.residual = rhs - A*p;
  broadcast(
    [&](Index I, CGSolver& cg, PDESystem& s, LaplaceMatrixOperator A) {
      cg.residual[I] = s.rhs[I] - A(s.p, I);
    },
    system.p.range, cg, system, A);

  residual_norm = dot(cg.residual, cg.residual);

  // cg.search_direction = cg.residual;
  broadcast(
    [&](Index I, CGSolver& cg, LaplaceMatrixOperator A) {
      cg.search_direction[I] = cg.residual[I];
    },
    system.p.range, cg, A);

  for (int iter = 0; iter < system.settings.maximumNumberOfIterations; iter++)
  {
    old_residual_norm = residual_norm;

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

    if (cg.residual.max() < system.settings.epsilon)
    {

      std::cout << "Residual: " << cg.residual.max() << std::endl;
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

void solve(GaussSeidelSolver gs, PDESystem& system)
{
  system.residual = 0;
  // std::cout << std::scientific << std::endl;
  // std::cout << std::endl;
  // std::cout << "Max Pressure: \t" << system.p.max() << "\n";
  // std::cout << "Min Pressure: \t" << system.p.min() << "\n";
  // std::cout << "Max Velocity x: \t" << system.u.max() << "\n";
  // std::cout << "Min Velocity x: \t" << system.u.min() << "\n";
  // std::cout << "Max Velocity y: \t" << system.v.max() << "\n";
  // std::cout << "Min Velocity y: \t" << system.v.min() << "\n";
  // std::cout << "RHS: \t" << system.rhs << "\n";
  // std::cout << "F: \t" << system.F << "\n";
  // std::cout << "G: \t" << system.G << "\n";
  // std::cout << "v: \t" << system.v << "\n";
  // std::cout << "u: \t" << system.u << "\n";
  // std::cout << "p: \t" << system.p << "\n";
  // std::cout << std::endl;
  for (int iter = 0; iter < 10000; iter++)
  {
    system.residual = 0;
    broadcast_boundary(
      [&](PDESystem& s, Index I, Offset o) { s.p[I + o] = s.p[I]; },
      system, system.p);
    broadcast(gauss_seidel_step, system, system.p.range);
    // std::cout << "Residual : \t" << system.residual << "\t\r"
    //           << std::flush;
    if (system.residual < 1e-4)
      break;
  }
  // gauss_seidel(system);
}
