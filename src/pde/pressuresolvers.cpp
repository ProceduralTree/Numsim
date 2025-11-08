#include "pde/pressuresolvers.h"
#include "linalg/matrix.h"
#include "linalg/vector.h"
#include "output/vtk.h"
#include "utils/broadcast.h"
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
    system.p.range, system, cg, A);

  residual_norm = dot(cg.residual, cg.residual);

  // cg.search_direction = cg.residual;
  broadcast(
    [&](Index I, CGSolver& cg, LaplaceMatrixOperator A) {
      cg.search_direction[I] = cg.residual[I];
    },
    system.p.range, system, cg, A);

  for (int iter = 0; iter < system.settings.maximumNumberOfIterations; iter++)
  {
    old_residual_norm = residual_norm;

    double alpha = residual_norm / Adot(A, cg.search_direction, cg.search_direction);

    // system.p = system.p + a * cg.search_direction;
    broadcast(
      [&](Index I, CGSolver& cg, PDESystem& s, double a) {
        system.p[I] = system.p[I] + a * cg.search_direction[I];
      },
      system.p.range, system, cg, alpha);

    // cg.residual = cg.residual - a * A * cg.search_direction;
    broadcast(
      [&](Index I, CGSolver& cg, PDESystem& s, double a, LaplaceMatrixOperator A) {
        cg.residual[I] = cg.residual[I] - a * A(cg.search_direction, I);
      },
      system.p.range, system, cg, alpha, A);

    if (cg.residual.max() < system.settings.epsilon)
      break;
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
