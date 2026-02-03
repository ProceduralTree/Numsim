#include "pde/pressuresolvers.h"
#include "grid/boundary.h"
#include "linalg/matrix.h"
#include "utils/profiler.h"
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

  SparseMatrixOperator A = SparseMatrixOperator(system.h, system.adjacency_map);

  // cg.residual = system.rhs - A*system.p;
  broadcast_boundary(copy_with_offset, system.boundary, static_cast<uint16_t>(BoundaryType::P_BOUNDARY), system.p);
  //  cg.residual[I] = s.rhs[I] - A(s.p, I);
  mipmap(_mean<double>, system.p);
  tree_broadcast(SparseVector::aAxpy, system.boundary, static_cast<uint16_t>(BoundaryType::P_Inside), cg.residual, -1., A, system.p, system.rhs);
  mipmap(_mean<double>, cg.residual);
  // A.a_ij modification for diagonal jacoby preconditionerS
  //  cg.residual[I] =1/A.a_ij[I] * cg.residual[I]
  // tree_broadcast(SparseVector::axpy, system.boundary, static_cast<uint16_t>(BoundaryType::P_Inside), cg.residual, (1 / A.a_ij - 1.), cg.residual, cg.residual);
  tree_broadcast(SparseVector::precondition, system.boundary, static_cast<uint16_t>(BoundaryType::P_Inside), cg.residual, A, cg.residual);
  mipmap(_mean<double>, cg.residual);
  residual_norm = SparseVector::dot(cg.residual, cg.residual, system.boundary);

  // ensure correct ghosts
  // cg.search_direction = cg.residual;
  broadcast(copy, system.boundary, static_cast<uint16_t>(BoundaryType::P_Inside), Offset { 0, 0 }, cg.residual, cg.search_direction);
  mipmap(_mean<double>, cg.search_direction);

  for (int iter = 0; iter < system.settings.maximumNumberOfIterations; iter++)
  {
    ProfileScope("CG Iteration");
    old_residual_norm = residual_norm;

    broadcast_boundary(copy_with_offset, system.boundary, static_cast<uint16_t>(BoundaryType::P_BOUNDARY), cg.search_direction);
    mipmap(_mean<double>, cg.search_direction);

    double alpha = residual_norm / SparseVector::dot(cg.search_direction, A, cg.search_direction, system.boundary);

    // A.a_ij modification for pcg mit diagonal jacoby preconditioner
    //  system.p = system.p + a * cg.search_direction;
    tree_broadcast(SparseVector::axpy_with_jacoby_precondition, system.boundary, static_cast<uint16_t>(BoundaryType::P_Inside), system.p, alpha, A, cg.search_direction, system.p);
    // tree_broadcast(SparseVector::axpy, system.boundary, static_cast<uint16_t>(BoundaryType::P_Inside), system.p, A[0] * alpha, cg.search_direction, system.p);

    // cg.residual = cg.residual - a * A * cg.search_direction;
    tree_broadcast(SparseVector::aAxpy, system.boundary, static_cast<uint16_t>(BoundaryType::P_Inside), cg.residual, -alpha, A, cg.search_direction, cg.residual);
    mipmap(_mean<double>, cg.residual);

    ProfilePush("Residual Calculation");
    double residual = SparseVector::max(cg.residual, system.boundary);
    // double min = SparseVector::min(cg.residual, system.boundary);
    ProfilePop();
    // std::cout << std::format("residual {}", residual) << std::endl;
    //  std::cout << std::format("min {}", min) << std::endl;
    // std::cout << std::format("residual Norm {}", residual_norm) << std::endl;
    // std::cout << std::format("alpha {}", alpha) << std::endl;
    // std::cout << std::format("<x,Ax> {}", SparseVector::dot(cg.search_direction, A, cg.search_direction, system.boundary)) << std::endl;
    // std::cout << std::format("<x,x> {}", SparseVector::dot(cg.search_direction, cg.search_direction, system.boundary)) << std::endl;
    if (residual > 1e6 || residual == -NAN || residual == NAN)
    {
      std::cout << std::format("residual exploded {}", residual) << std::endl;
      // ErrorF("residual exploded {}", residual);
      break;
      // abort();
    }
    if (iter >= system.settings.maximumNumberOfIterations - 2)
    {
      std::cerr << std::format("solver did not converge residual: {}", residual) << std::endl;
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
  mipmap(_mean<double>, system.p);
  broadcast_boundary(copy_with_offset, system.boundary, static_cast<uint16_t>(BoundaryType::P_BOUNDARY), system.p);
}

using std::swap;
void solve(Jacoby& S, PDESystem& system)
{
  ProfileScope("Jacoby Solver");
  system.residual = 0;
  broadcast(set, system.boundary, static_cast<uint16_t>(BoundaryType::P_Inside), Offset { 0, 0 }, S.residual, NAN);
  for (int iter = 0; iter < Settings::get().maximumNumberOfIterations; iter++)
  {
    broadcast_boundary(copy_with_offset, system.boundary, static_cast<uint16_t>(BoundaryType::P_BOUNDARY), system.p);
    mipmap(_mean<double>, system.p);
    tree_broadcast(jacoby_step_adj, system.boundary, static_cast<uint16_t>(BoundaryType::P_Inside), S.tmp, S.residual, system.rhs, system.p, system.h, system.adjacency_map);
    mipmap(_mean<double>, S.tmp);
    swap(system.p, S.tmp);
    if (iter % 100 && SparseVector::max(S.residual, system.boundary) < Settings::get().epsilon)
    {
      DebugF("Residual {:.14e} \nJacobi converged after n={}", SparseVector::max(S.residual, system.boundary), iter);
      break;
    }
  }
}
