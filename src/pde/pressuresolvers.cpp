#include "pde/pressuresolvers.h"
#include "grid/boundary.h"
#include "linalg/matrix.h"
#include <cstdint>
#include <grid/grid.h>
#include <linalg/sparsevector.h>
#include <mpi.h>
#include <pde/system.h>
#include <utils/broadcast.h>

void solve(CGSolver& cg, PDESystem& system)
{
  double residual_norm = INFINITY;
  double old_residual_norm = INFINITY;

  SparseMatrixOperator A = SparseMatrixOperator(system.h);

  // cg.residual = system.rhs - A*system.p;
  broadcast_boundary(copy_with_offset, system.boundary, static_cast<uint16_t>(BoundaryType::P_BOUNDARY), system.p);
  //  cg.residual[I] = s.rhs[I] - A(s.p, I);
  mipmap(_mean<double>, system.p);
  tree_broadcast(SparseVector::aAxpy, system.boundary, static_cast<uint16_t>(BoundaryType::P_Inside), cg.residual, -1., A, system.p, system.rhs);
  mipmap(_mean<double>, cg.residual);
  // A.a_ij modification for diagonal jacoby preconditioner
  tree_broadcast(SparseVector::axpy, system.boundary, static_cast<uint16_t>(BoundaryType::P_Inside), cg.residual, (1 / A.a_ij - 1.), cg.residual, cg.residual);
  mipmap(_mean<double>, cg.residual);
  residual_norm = SparseVector::dot(cg.residual, cg.residual, system.boundary);

  // ensure correct ghosts
  // cg.search_direction = cg.residual;
  broadcast(copy, system.boundary, static_cast<uint16_t>(BoundaryType::P_Inside), Offset { 0, 0 }, cg.residual, cg.search_direction);

  for (int iter = 0; iter < system.settings.maximumNumberOfIterations; iter++)
  {
    ProfileScope("CG Iteration");
    old_residual_norm = residual_norm;

    broadcast_boundary(copy_with_offset, system.boundary, static_cast<uint16_t>(BoundaryType::P_BOUNDARY), cg.search_direction);

    double alpha = residual_norm / SparseVector::dot(cg.search_direction, A, cg.search_direction, system.boundary);

    // A.a_ij modification for pcg mit diagonal jacoby preconditioner
    //  system.p = system.p + a * cg.search_direction;
    tree_broadcast(SparseVector::axpy, system.boundary, static_cast<uint16_t>(BoundaryType::P_Inside), system.p, A.a_ij * alpha, cg.search_direction, system.p);

    // cg.residual = cg.residual - a * A * cg.search_direction;
    tree_broadcast(SparseVector::aAxpy, system.boundary, static_cast<uint16_t>(BoundaryType::P_Inside), cg.residual, -alpha, A, cg.search_direction, cg.residual);
    mipmap(_mean<double>, cg.residual);

    ProfilePush("Residual Calculation");
    double residual = SparseVector::max(cg.residual, system.boundary);
    double min = SparseVector::min(cg.residual, system.boundary);
    ProfilePop();
    // std::cout << std::format("residual {}", residual) << std::endl;
    // std::cout << std::format("min {}", min) << std::endl;
    // std::cout << std::format("residual Norm {}", residual_norm) << std::endl;
    // std::cout << std::format("alpha {}", alpha) << std::endl;
    // std::cout << std::format("A_ij {}", A.a_ij) << std::endl;
    // std::cout << std::format("1/h_x^2 {}", A.h_x_squared_inv) << std::endl;
    if (residual > 1e6 || residual == -NAN || residual == NAN)
    {
      std::cout << std::format("residual exploded {}", residual) << std::endl;
      // ErrorF("residual exploded {}", residual);
      break;
      // abort();
    }

    if (residual < Settings::get().epsilon)
    {
      break;
    }
    residual_norm = SparseVector::dot(cg.residual, cg.residual, system.boundary);
    // DebugF("Residual Norm : {}", residual_norm);
    double beta = residual_norm / old_residual_norm;
    // DebugF("Beta: {}", beta);

    // TODO Update Ghosts
    // cg.search_direction[I] = cg.residual[I] + beta * cg.search_direction[I];
    tree_broadcast(SparseVector::axpy, system.boundary, static_cast<uint16_t>(BoundaryType::P_Inside), cg.search_direction, beta, cg.search_direction, cg.residual);
    mipmap(_mean<double>, cg.search_direction);
  }
  broadcast_boundary(copy_with_offset, system.boundary, static_cast<uint16_t>(BoundaryType::P_BOUNDARY), system.p);
  mipmap(_mean<double>, system.p);
}

// void solve(GaussSeidelSolver& S, PDESystem& system)
//{
//   system.residual = 0;
//   for (int iter = 0; iter < Settings::get().maximumNumberOfIterations; iter++)
//   {
//     system.residual = 0;
//     broadcast_boundary(copy_with_offset, system.partitioning, system.p.boundary, system.p);
//     broadcast(gauss_seidel_step, system.p.range, system, S);
//     if (system.residual < Settings::get().epsilon)
//     {
//       DebugF("Residual {:.14e} \nconverged after n={}", system.residual, iter);
//       break;
//     }
//   }
// }
//
// void solve(SORSolver& S, PDESystem& system)
//{
//   Index gridpos = system.partitioning.getGridPos();
//   int parity = (gridpos.x + gridpos.y) % 2;
//   for (int iter = 0; iter < Settings::get().maximumNumberOfIterations; iter++)
//   {
//     ProfileScope("SOR Iteration");
//     system.residual = 0;
//     broadcast_boundary(copy_with_offset, system.partitioning, system.p.boundary, system.p);
//     broadcast_blackred(sor_step, parity, system.p.range, system);
//     MPI_COMM_BUFFER* comm_black = new MPI_COMM_BUFFER(system.p, system.p.boundary.all, MPI_COMM_WORLD, system.partitioning);
//     delete comm_black;
//     broadcast_blackred(sor_step, !parity, system.p.range, system);
//     MPI_COMM_BUFFER* comm_red = new MPI_COMM_BUFFER(system.p, system.p.boundary.all, MPI_COMM_WORLD, system.partitioning);
//     delete comm_red;
//
//     double local_residual = system.residual;
//     double global_residual = 0.;
//     MPI_Allreduce(&local_residual, &global_residual, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
//
//     if (global_residual > 1e16)
//     {
//       ErrorF("residual exploded {}", global_residual);
//
//       for (int i = 0; i < system.partitioning.size; i++)
//       {
//         if (system.partitioning.rank == i)
//         {
//           std::cout << "Hello from Rank " << system.partitioning.rank << " of " << system.partitioning.size << std::endl;
//           std::cout << "Pressure: " << system.p << std::endl;
//           std::cout << "RHS: " << system.rhs << std::endl;
//         }
//       }
//       abort();
//     }
//
//     if (iter % 10 && global_residual < Settings::get().epsilon)
//     {
//       // std::cout << "COnverged after N=" << iter << " Iterations" << std::endl;
//
//       // std::cout << std::scientific << std::setprecision(14) << "Residual: " << system.residual << std::endl;
//       //  std::cout << "\nSOR converged after n=" << iter << " Iterations" << std::endl;
//       break;
//     }
//   }
// }
// void solve(BlackRedSolver& S, PDESystem& system)
//{
//   system.residual = 0;
//   parallel_broadcast(set, system.p.range, Offset { 0, 0 }, S.residual, INFINITY);
//   for (int iter = 0; iter < Settings::get().maximumNumberOfIterations; iter++)
//   {
//     broadcast_boundary(copy_with_offset, system.partitioning, system.p.boundary, system.p);
//     broadcast_blackred(black_red_step, 0, system.p.range, system, S);
//     broadcast_blackred(black_red_step, 1, system.p.range, system, S);
//     if (iter % 100 && S.residual.max() < Settings::get().epsilon)
//     {
//
//       DebugF("Residual {:.14e} \nBlack Red converged after n={}", S.residual.max(), iter);
//       break;
//     }
//   }
// }
//
// void solve(Jacoby& S, PDESystem& system)
//{
//   system.residual = 0;
//   parallel_broadcast(set, system.p.range, Offset { 0, 0 }, S.residual, INFINITY);
//   for (int iter = 0; iter < Settings::get().maximumNumberOfIterations; iter++)
//   {
//     broadcast_boundary(copy_with_offset, system.partitioning, system.p.boundary, system.p);
//     test_broadcast(jacoby_step, system.p.range, system, S);
//     std::swap(system.p, S.tmp);
//     if (iter % 100 && S.residual.max() < Settings::get().epsilon)
//     {
//       DebugF("Residual {:.14e} \nJacobi converged after n={}", S.residual.max(), iter);
//       break;
//     }
//   }
// }
