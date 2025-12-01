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
#include <cstdlib>
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
  // MPI_COMM_BUFFER* comm_p = new MPI_COMM_BUFFER(system.p, system.p.boundary.all, MPI_COMM_WORLD, system.partitioning);
  // delete comm_p;

  broadcast_boundary(copy_with_offset, system.partitioning, system.p.boundary, system.p);
  // broadcast_boundary(copy, system.partitioning, system.p.boundary, system.p, cg.search_direction);
  //
  //  cg.residual[I] = s.rhs[I] - A(s.p, I);
  broadcast(aAxpy, system.p.range, cg.residual, -1., A, system.p, system.rhs);

  residual_norm = dot(cg.residual, cg.residual);

  // ensure correct ghosts
  // cg.search_direction = cg.residual;
  distributed_broadcast(copy, system.partitioning, system.p.range, cg.search_direction, Offset { 0, 0 }, cg.residual, cg.search_direction);

  for (int iter = 0; iter < Settings::get().maximumNumberOfIterations; iter++)
  {
    ProfileScope("CG Iteration");
    old_residual_norm = residual_norm;

    broadcast_boundary(copy_with_offset, system.partitioning, cg.search_direction.boundary, cg.search_direction);

    double alpha = residual_norm / Adot(A, cg.search_direction, cg.search_direction);
    // DebugF("Alpha: {}", alpha);

    // system.p = system.p + a * cg.search_direction;
    broadcast(axpy, { system.p.begin, system.p.end }, system.p, alpha, cg.search_direction, system.p);
    // broadcast_ghosts(axpy, system.partitioning, system.p.boundary, system.p, alpha, cg.search_direction, system.p);
    //  broadcast(axpy, system.p.boundary.all, system.p, alpha, cg.search_direction, system.p);

    // cg.residual = cg.residual - a * A * cg.search_direction;
    broadcast(aAxpy, system.p.range, cg.residual, -alpha, A, cg.search_direction, cg.residual);

    double residual = cg.residual.max();
    if (residual > 1e5 || residual == -NAN || residual == NAN)
    {
      ErrorF("residual exploded {}", residual);

      for (int i = 0; i < system.partitioning.size; i++)
      {
        MPI_Barrier(MPI_COMM_WORLD);
        if (system.partitioning.rank == i)
        {
          std::cout << "Hello from Rank " << system.partitioning.rank << " of " << system.partitioning.size << std::endl;
          DebugF("Alpha: {}", alpha);
          std::cout << "CG Solver res" << cg.residual << std::endl;
          std::cout << "CG Solver p" << cg.search_direction << std::endl;
          std::cout << "Pressure: " << system.p << std::endl;
          std::cout << "RHS: " << system.rhs << std::endl;
        }
      }
      MPI_Barrier(MPI_COMM_WORLD);
      abort();
    }

    if (residual < Settings::get().epsilon)
    {
      // update Pressure ghosts
      auto comm_buffer = new MPI_COMM_BUFFER(system.p, system.p.boundary.all, MPI_COMM_WORLD, system.partitioning);
      delete comm_buffer;
      DebugF("COnverged after {} Iterations", iter);
      break;
    }
    residual_norm = dot(cg.residual, cg.residual);
    // DebugF("Residual Norm : {}", residual_norm);
    double beta = residual_norm / old_residual_norm;
    // DebugF("Beta: {}", beta);

    // TODO Update Ghosts
    // cg.search_direction[I] = cg.residual[I] + beta * cg.search_direction[I];
    distributed_broadcast(axpy, system.partitioning, system.p.range, cg.search_direction, cg.search_direction, beta, cg.search_direction, cg.residual);
    MPI_Barrier(MPI_COMM_WORLD);
  }
  MPI_Barrier(MPI_COMM_WORLD);
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
    broadcast_black(sor_step, system.p.range, system);
    MPI_COMM_BUFFER* comm_black = new MPI_COMM_BUFFER(system.p, system.p.boundary.all, MPI_COMM_WORLD, system.partitioning);
    delete comm_black;
    broadcast_red(sor_step, system.p.range, system);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_COMM_BUFFER* comm_red = new MPI_COMM_BUFFER(system.p, system.p.boundary.all, MPI_COMM_WORLD, system.partitioning);
    delete comm_red;

    double local_residual = system.residual;
    double global_residual = 0.;
    MPI_Allreduce(&local_residual, &global_residual, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    if (global_residual > 1e16)
    {
      ErrorF("residual exploded {}", global_residual);

      for (int i = 0; i < system.partitioning.size; i++)
      {
        MPI_Barrier(MPI_COMM_WORLD);
        if (system.partitioning.rank == i)
        {
          std::cout << "Hello from Rank " << system.partitioning.rank << " of " << system.partitioning.size << std::endl;
          std::cout << "Pressure: " << system.p << std::endl;
          std::cout << "RHS: " << system.rhs << std::endl;
        }
      }
      abort();
    }

    if (iter % 10 && global_residual < Settings::get().epsilon)
    {

      // std::cout << std::scientific << std::setprecision(14) << "Residual: " << system.residual << std::endl;
      //  std::cout << "\nSOR converged after n=" << iter << " Iterations" << std::endl;
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
