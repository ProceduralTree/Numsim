#include "pde/pressuresolvers.h"
#include "linalg/matrix.h"
#include "linalg/vector.h"
#include "utils/Logger.h"
#include "utils/broadcast.h"
#include "utils/distributed.h"
#include "utils/index.h"
#include "utils/profiler.h"
#include "utils/settings.h"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <grid/grid.h>
#include <grid/indexing.h>
#include <ios>
#include <mpi.h>
#include <pde/system.h>

void solve(CGSolver& cg, PDESystem& system)
{
  double residual_norm = INFINITY;
  double old_residual_norm = INFINITY;

  LaplaceMatrixOperator A = LaplaceMatrixOperator(system.h);

  // cg.residual = system.rhs - A*system.p;
  //  broadcast_boundary(copy_with_offset, system.p.boundary, system.p);
  broadcast_boundary(copy_with_offset, system.partitioning, system.p.boundary, system.p);
  broadcast_halo(copy, system.p.boundary, system.p, cg.search_direction);
  // cg.residual[I] = s.rhs[I] - A(s.p, I);
  parallel_broadcast(aAxpy, system.p.range, cg.residual, -1., A, system.p, system.rhs);

  residual_norm = dot(cg.residual, cg.residual);

  // cg.search_direction = cg.residual;
  parallel_broadcast(copy, system.p.range, Offset { 0, 0 }, cg.residual, cg.search_direction);

  for (int iter = 0; iter < Settings::get().maximumNumberOfIterations; iter++)
  {
    ProfileScope("CG Iteration");
    old_residual_norm = residual_norm;

    broadcast_boundary(copy_with_offset, system.partitioning, cg.search_direction.boundary, cg.search_direction);
    // broadcast_boundary(copy, system.p.boundary, system.p, cg.search_direction);
    double alpha = residual_norm / Adot(A, cg.search_direction, cg.search_direction);
    DebugF("Alpha: {}", alpha);

    // system.p = system.p + a * cg.search_direction;
    broadcast(axpy, system.p.range, system.p, alpha, cg.search_direction, system.p);

    // cg.residual = cg.residual - a * A * cg.search_direction;
    broadcast(aAxpy, system.p.range, cg.residual, -alpha, A, cg.search_direction, cg.residual);

    // std::cout << "\rResidual:\t" << residual_norm << " Iterations:\t" << iter << std::flush;
    if (cg.residual.max() < Settings::get().epsilon)
    {
      // update Pressure ghosts
      auto comm_buffer = new MPI_COMM_BUFFER(system.p, system.p.boundary.all, MPI_COMM_WORLD, system.partitioning);
      delete comm_buffer;
      DebugF("Residual {:.14e} \nconverged after n={}", cg.residual.max(), iter);
      break;
    }
    residual_norm = dot(cg.residual, cg.residual);
    double beta = residual_norm / old_residual_norm;

    // TODO Update Ghosts
    // cg.search_direction[I] = cg.residual[I] + beta * cg.search_direction[I];
    distributed_broadcast(axpy, system.partitioning, system.p.range, cg.search_direction, cg.search_direction, beta, cg.search_direction, cg.residual);
    // broadcast(axpy, system.p.range, cg.search_direction, beta, cg.search_direction, cg.residual);
  }
  broadcast_boundary(copy_with_offset, system.partitioning, system.p.boundary, system.p);
}

void solve(GaussSeidelSolver& S, PDESystem& system)
{
  system.residual = 0;
  for (int iter = 0; iter < Settings::get().maximumNumberOfIterations; iter++)
  {
    system.residual = 0;
    broadcast_boundary(copy_with_offset, system.partitioning, system.p.boundary, system.p);
    broadcast(gauss_seidel_step, system.p.range, system, S);
    if (system.residual < Settings::get().epsilon)
    {
      DebugF("Residual {:.14e} \nconverged after n={}", system.residual, iter);
      break;
    }
  }
}

void solve(SORSolver& S, PDESystem& system)
{
  for (int iter = 0; iter < Settings::get().maximumNumberOfIterations; iter++)
  {
    ProfileScope("SOR Iteration");
    system.residual = 0;
    broadcast_boundary(copy_with_offset, system.partitioning, system.p.boundary, system.p);
    broadcast(sor_step, system.p.range, system);
    if (iter % 10 && system.residual < Settings::get().epsilon)
    {

      // std::cout << std::scientific << std::setprecision(14) << "Residual: " << system.residual << std::endl;
      // std::cout << "\nSOR converged after n=" << iter << " Iterations" << std::endl;
      break;
    }
  }
}
void solve(BlackRedSolver& S, PDESystem& system)
{
  system.residual = 0;
  parallel_broadcast(set, system.p.range, Offset { 0, 0 }, S.residual, INFINITY);
  for (int iter = 0; iter < Settings::get().maximumNumberOfIterations; iter++)
  {
    broadcast_boundary(copy_with_offset, system.partitioning, system.p.boundary, system.p);
    broadcast_black(black_red_step, system.p.range, system, S);
    broadcast_red(black_red_step, system.p.range, system, S);
    if (iter % 100 && S.residual.max() < Settings::get().epsilon)
    {

      DebugF("Residual {:.14e} \nBlack Red converged after n={}", S.residual.max(), iter);
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
    broadcast_boundary(copy_with_offset, system.partitioning, system.p.boundary, system.p);
    test_broadcast(jacoby_step, system.p.range, system, S);
    std::swap(system.p, S.tmp);
    if (iter % 100 && S.residual.max() < Settings::get().epsilon)
    {
      DebugF("Residual {:.14e} \nJacobi converged after n={}", S.residual.max(), iter);
      break;
    }
  }
}
