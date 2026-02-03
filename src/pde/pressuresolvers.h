#ifndef PRESSURESOLVERS_H_
#define PRESSURESOLVERS_H_

#include "grid/adjacencymap.h"
#include "grid/densetree.h"
#include "grid/sparsegrid.h"
#include "grid/zindex.h"
#include <cstddef>
#include <cstdint>
#include <pde/system.h>
#include <utils/index.h>

struct CGSolver
{
  SparseGrid2D<double> residual;
  SparseGrid2D<double> search_direction;
  CGSolver(const DenseTree::DenseTree& tree)
    : residual(tree)
    , search_direction(tree) {
    };
};

// struct GaussSeidelSolver
//{
//   // Grid2D residual;
//   // GaussSeidelSolver(PDESystem& system)
//   //   : residual(system.begin, system.end) { };
// };
//
// struct SORSolver
//{
//   // Grid2D residual;
//   // SORSolver(PDESystem& system)
//   //   : residual(system.begin, system.end) { };
// };
// struct BlackRedSolver
//{
//   Grid2D residual;
//   BlackRedSolver(PDESystem& system)
//     : residual(system.p.begin, system.p.end) { };
// };
struct Jacoby
{
  SparseGrid2D<double> residual;
  SparseGrid2D<double> tmp;
  Jacoby(const DenseTree::DenseTree& tree)
    : residual(tree)
    , tmp(tree) { };
};

void solve(Jacoby& S, PDESystem& system);
void solve(CGSolver& S, PDESystem& system);

inline void copy_with_offset(Index I, Offset O, SparseGrid2D<double>& array) { array[I] = array[I + O]; };

constexpr void jacoby_step(size_t index, uint16_t depth, SparseGrid2D<double>& result, SparseGrid2D<double>& residual, const SparseGrid2D<double>& rhs, const SparseGrid2D<double>& p, Gridsize h)
{
  double sum_of_neighbours = 0;
  double local_cell_size = static_cast<double>(1ULL << (depth));
  Index I = ZorderToIndex(p.tree._index_cache[index]);
  sum_of_neighbours += ((p[I - Ix] + p[I + local_cell_size * Ix]) / h.x_squared(depth)); //+ ((p[I - Iy] + p[I + Iy]) / h.y_squared(depth));
  sum_of_neighbours += ((p[I - Iy] + p[I + local_cell_size * Iy]) / h.y_squared(depth));
  double a_ij = -2. * (1. / h.y_squared(depth)) - 2. * (1. / h.x_squared(depth));
  residual[index] = std::abs(sum_of_neighbours + a_ij * p[index] - rhs[index]);
  result[index] = (rhs[index] - sum_of_neighbours) / a_ij;
}
constexpr void jacoby_step_adj(size_t index, uint16_t depth, SparseGrid2D<double>& result, SparseGrid2D<double>& residual, const SparseGrid2D<double>& rhs, const SparseGrid2D<double>& p, Gridsize h, const AdjMap& adj)
{
  double sum_of_neighbours = 0;
  size_t top = adj._top[index];
  size_t bottom = adj._bottom[index];
  size_t left = adj._left[index];
  size_t right = adj._right[index];
  sum_of_neighbours += ((p[left] + p[right]) / h.x_squared(depth));
  sum_of_neighbours += ((p[bottom] + p[top]) / h.y_squared(depth));
  double a_ij = -2. * (1. / h.y_squared(depth)) - 2. * (1. / h.x_squared(depth));
  residual[index] = std::abs(sum_of_neighbours + a_ij * p[index] - rhs[index]);
  result[index] = (rhs[index] - sum_of_neighbours) / a_ij;
}
//
// inline void gauss_seidel_step(Index I, PDESystem& system, GaussSeidelSolver& S)
//{
//   auto& p = system.p;
//   auto& h = system.h;
//   double sum_of_neighbours = ((p[I - Ix] + p[I + Ix]) / h.x_squared) + ((p[I - Iy] + p[I + Iy]) / h.y_squared);
//   double a_ij = -2 * (1 / h.y_squared) - 2 * (1 / h.x_squared);
//   system.residual = std::max(std::abs(sum_of_neighbours + a_ij * p[I] - system.rhs[I]), system.residual);
//   p[I] = (system.rhs[I] - sum_of_neighbours) / a_ij;
// };
//
// inline void sor_step(Index I, PDESystem& system)
//{
//   auto& p = system.p;
//   auto& h = system.h;
//   double sum_of_neighbours = ((p[I - Ix] + p[I + Ix]) / h.x_squared) + ((p[I - Iy] + p[I + Iy]) / h.y_squared);
//   double a_ij = -2 * (1 / h.y_squared) - 2 * (1 / h.x_squared);
//   double residual = std::abs(sum_of_neighbours + a_ij * p[I] - system.rhs[I]);
//   system.residual = std::max(residual, system.residual);
//   p[I] = (1 - Settings::get().omega) * p[I] + Settings::get().omega * (system.rhs[I] - sum_of_neighbours) / a_ij;
// };
// inline void black_red_step(Index I, PDESystem& system, BlackRedSolver& solver)
//{
//   auto [up, res] = jacoby_update(I, system);
//   // DebugF("Update {} , residual {}", up, res);
//   solver.residual[I]
//     = res;
//   system.p[I] = up;
// };
// inline void jacoby_step(Index I, const PDESystem& system, Jacoby& solver)
//{
//   auto [up, res] = jacoby_update(I, system);
//   // DebugF("Update {} , residual {}", up, res);
//   solver.residual[I] = res;
//   solver.tmp[I] = up;
// };

#endif // PRESSURESOLVERS_H_
