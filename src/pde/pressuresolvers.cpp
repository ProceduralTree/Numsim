#include "output/vtk.h"
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

void sor_step(PDESystem& system, Index I)
{
  auto& p = system.p;
  auto& h = system.h;
  double sum_of_neighbours = ((p[I - Ix] + p[I + Ix]) / h.x_squared) + ((p[I - Iy] + p[I + Iy]) / h.y_squared);
  double a_ij = -2 * (1 / h.y_squared) - 2 * (1 / h.x_squared);
  double residual = std::abs(sum_of_neighbours + a_ij * p[I] - system.rhs[I]);
  system.residual = std::max(residual, system.residual);
  p[I] = (1 - system.settings.omega) * p[I] + system.settings.omega * (system.rhs[I] - sum_of_neighbours) / a_ij;
};