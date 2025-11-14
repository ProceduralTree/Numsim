#include "pde/pressuresolvers.h"
#include "linalg/matrix.h"
#include "linalg/vector.h"
#include "utils/Logger.h"
#include "utils/broadcast.h"
#include "utils/index.h"
#include "utils/settings.h"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <grid/grid.h>
#include <grid/indexing.h>
#include <ios>
#include <pde/system.h>

void solve(CGSolver& cg, PDESystem& system)
{
  double residual_norm = INFINITY;
  double old_residual_norm = INFINITY;

  LaplaceMatrixOperator A = LaplaceMatrixOperator(system.h);

  // cg.residual = system.rhs - A*system.p;
  //  broadcast_boundary(copy_with_offset, system.p.boundary, system.p);
  broadcast_boundary(copy_with_offset, system.p.boundary, system.p);
  broadcast_boundary(copy, system.p.boundary, system.p, cg.search_direction);
  // cg.residual[I] = s.rhs[I] - A(s.p, I);
  broadcast(aAxpy, system.p.range, cg.residual, -1., A, system.p, system.rhs);

  residual_norm = dot(cg.residual, cg.residual);

  // cg.search_direction = cg.residual;
  broadcast(copy, system.p.range, Offset { 0, 0 }, cg.residual, cg.search_direction);

  for (int iter = 0; iter < Settings::get().maximumNumberOfIterations; iter++)
  {

    old_residual_norm = residual_norm;

    broadcast_boundary(copy_with_offset, cg.search_direction.boundary, cg.search_direction);
    // broadcast_boundary(copy, system.p.boundary, system.p, cg.search_direction);
    double alpha = residual_norm / Adot(A, cg.search_direction, cg.search_direction);

    // system.p = system.p + a * cg.search_direction;
    broadcast(axpy, system.p.range, system.p, alpha, cg.search_direction, system.p);

    // cg.residual = cg.residual - a * A * cg.search_direction;
    broadcast(aAxpy, system.p.range, cg.residual, -alpha, A, cg.search_direction, cg.residual);

    // std::cout << "\rResidual:\t" << residual_norm << " Iterations:\t" << iter << std::flush;
    if (cg.residual.max() < Settings::get().epsilon)
    {
      std::cout << std::scientific << std::setprecision(14) << "Residual: " << cg.residual.max() << std::endl;
      std::cout << "converged after n=" << iter << " Iterations" << std::endl;
      break;
    }
    residual_norm = dot(cg.residual, cg.residual);
    double beta = residual_norm / old_residual_norm;

    // cg.search_direction[I] = cg.residual[I] + beta * cg.search_direction[I];
    broadcast(axpy, system.p.range, cg.search_direction, beta, cg.search_direction, cg.residual);
  }
  broadcast_boundary(copy_with_offset, system.p.boundary, system.p);
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

      std::cout << std::scientific << std::setprecision(14) << "Residual: " << system.residual << std::endl;
      std::cout << "\nSOR converged after n=" << iter << " Iterations" << std::endl;
      break;
    }
  }
}
void solve(BlackRedSolver& S, PDESystem& system)
{
  system.residual = 0;
  for (int iter = 0; iter < Settings::get().maximumNumberOfIterations; iter++)
  {
    parallel_broadcast(set, system.p.range, Offset { 0, 0 }, S.residual, INFINITY);
    broadcast_boundary(copy_with_offset, system.p.boundary, system.p);
    broadcast_blackred(black_red_step, system.p.range, system, S);
    if (iter % 100 && S.residual.max() < Settings::get().epsilon)
    {

      std::cout << std::scientific << std::setprecision(14) << "Residual: " << system.residual << std::endl;
      std::cout << "Black Red converged after n=" << iter << " Iterations" << std::endl;
      break;
    }
  }
}

void solve(Jacoby& S, PDESystem& system)
{
  system.residual = 0;
  parallel_broadcast(set, system.p.range, Offset { 0, 0 }, S.residual, INFINITY);
  for (int iter = 0; iter < Settings::get().maximumNumberOfIterations; iter++)
  {
    broadcast_boundary(copy_with_offset, system.p.boundary, system.p);
    auto r = system.p.range;
    auto& h = system.h;
    auto& p = system.p;
    auto hx2inv = 1 / h.x_squared;
    auto hy2inv = 1 / h.y_squared;
    double a_ij = -2. * (1. / h.y_squared) - 2. * (1. / h.x_squared);
    double a_ijinv = 1 / a_ij;
    // #pragma omp parallel for collapse(2)
    //     for (uint16_t j = 1; j <= 100; j++)
    //     {
    //       for (uint16_t i = 1; i <= 100; i++)
    //       {
    //         uint32_t I = i + p.size_x * j;
    //         uint32_t Ixp1 = (i + 1) + p.size_x * j;
    //         uint32_t Ixm1 = (i - 1) + p.size_x * j;
    //         uint32_t Iyp1 = i + p.size_x * (j + 1);
    //         uint32_t Iym1 = i + p.size_x * (j - 1);
    //         double sum_of_neighbours = ((p[Ixm1] + p[Ixp1]) * hx2inv) + ((p[Iym1] + p[Iyp1]) * hy2inv);
    //         double residual = std::abs(sum_of_neighbours + a_ij * p[I] - system.rhs[I]);
    //         double update = (system.rhs[I] - sum_of_neighbours) * a_ijinv;
    //         S.residual[I] = residual;
    //         S.tmp[I] = update;
    //       }
    //     }
    parallel_broadcast(jacoby_step, system.p.range, system, S);
    std::swap(system.p, S.tmp);
    if (iter % 100 && S.residual.max() < Settings::get().epsilon)
    {

      std::cout << std::scientific << std::setprecision(14) << "Residual: " << S.residual.max() << std::endl;
      std::cout << "Jacoby converged after n=" << iter << " Iterations" << std::endl;
      break;
    }
  }
}
