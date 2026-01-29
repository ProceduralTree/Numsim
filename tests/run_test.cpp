
#include "grid/boundary.h"
#include "grid/sparsegrid.h"
#include "linalg/matrix.h"
#include "output/vtk_tree.h"
#include "pde/system.h"
#include "utils/Logger.h"
#include "utils/index.h"
#include "utils/profiler.h"
#include "utils/settings.h"
#include <cassert>
#include <csignal>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <format>
#include <grid/densetree.h>
#include <grid/util.h>
#include <iostream>
#include <linalg/sparsevector.h>
#include <mpi.h>
#include <pde/pressuresolvers.h>

#define ASSERT(condition, message)                               \
  do                                                             \
  {                                                              \
    if (!(condition))                                            \
    {                                                            \
      std::cerr << "Assertion failed: " << message << std::endl; \
      assert(condition);                                         \
    }                                                            \
  } while (0)

void signalInt(int sig)
{
  DebugF("Interrupt from: {}", sig);
  Profiler::Close();
  exit(sig);
};

struct Range get_test_range()
{
  auto begin = Index { 1, 1, 0 };
  auto end = Index { 20, 20, 0 };
  return { begin, end };
}

DenseTree::DenseTree get_test_tree()
{

  return DenseTree::from_range(get_test_range());
};

void test_build_tree()
{

  auto t = get_test_tree();
  auto data_set = init(t);
  // ASSERT(tree.sizes, message)
  // t.print();
  write_depth("Tree Depth", t, data_set);
  save_dataset(data_set);
};

void test_tree_refinement()
{
  auto t = get_test_tree();

  auto updated_tree = DenseTree::build_tree(is_desired_depth, t.maxDepth, t.maxDepth, 4, t);
  updated_tree.print();
  auto data_set = init(updated_tree);
  write_depth("Updated Depth", updated_tree, data_set);
  save_dataset(data_set);
};

void test_set_cartesian_index()
{
  auto t = get_test_tree();
  auto updated_tree = DenseTree::build_tree(is_desired_depth, t.maxDepth, t.maxDepth, 4, t);
  auto sparse_grid = SparseGrid2D<double>(updated_tree);
  auto data_set = init(updated_tree);
  for (uint16_t i = 0; i < 1ULL << t.maxDepth; i++)
  {
    sparse_grid[{ i, i, updated_tree.maxDepth }] = 1. * i + 1.;
  }
  for (uint16_t i = 0; i < 1ULL << t.maxDepth; i++)
  {
    auto val = sparse_grid[{ i, i, updated_tree.maxDepth }];
    ASSERT((val == 1. * i + 1), (std::format("Did not set value at expected point val={} , i={}", val, i)));
  }
  write_field("Grid", sparse_grid.tree, sparse_grid._data, data_set);
  save_dataset(data_set);
};

void test_set_boundary()
{
  auto t = get_test_tree();
  auto updated_tree = DenseTree::build_tree(is_desired_depth, t.maxDepth, t.maxDepth, 4, t);
  Range r = get_test_range();
  BoundaryFlags b = BoundaryFlags(updated_tree, r);
  auto data_set = init(updated_tree);
  write_field("BoundaryFlags", b.tree, b.flags._data, data_set);
  save_dataset(data_set);
};

constexpr void _set(size_t index, uint16_t depth, SparseGrid2D<double>& grid, double value)
{
  grid[index] = value;
};

void test_set_values()
{
  auto t = get_test_tree();
  auto updated_tree = DenseTree::build_tree(is_desired_depth, t.maxDepth, t.maxDepth, t.maxDepth, t);
  Range r = get_test_range();
  BoundaryFlags b = BoundaryFlags(updated_tree, r);
  SparseGrid2D<double> ugrid = SparseGrid2D<double>(updated_tree);
  SparseGrid2D<double> vgrid = SparseGrid2D<double>(updated_tree);
  SparseGrid2D<double> pgrid = SparseGrid2D<double>(updated_tree);
  uint16_t u_boundary = static_cast<uint16_t>(BoundaryType::U_BOUNDARY);
  uint16_t v_boundary = static_cast<uint16_t>(BoundaryType::V_BOUNDARY);
  uint16_t p_boundary = static_cast<uint16_t>(BoundaryType::P_BOUNDARY);
  auto data_set = init(updated_tree);
  tree_broadcast(_set, b, u_boundary, ugrid, 1.);
  tree_broadcast(_set, b, v_boundary, vgrid, 1.);
  tree_broadcast(_set, b, p_boundary, pgrid, 1.);
  tree_broadcast(_set, b, static_cast<uint16_t>(BoundaryType::U_Inside), ugrid, -1.);
  tree_broadcast(_set, b, static_cast<uint16_t>(BoundaryType::V_Inside), vgrid, -1.);
  tree_broadcast(_set, b, static_cast<uint16_t>(BoundaryType::P_Inside), pgrid, -1.);
  write_field("U boundary", b.tree, ugrid._data, data_set);
  write_field("V boundary", b.tree, vgrid._data, data_set);
  write_field("P boundary", b.tree, pgrid._data, data_set);
  save_dataset(data_set);
};

void test_vector_operations()
{
  auto t = get_test_tree();
  auto updated_tree = DenseTree::build_tree(is_desired_depth, t.maxDepth, t.maxDepth, 4, t);
  Range r = get_test_range();
  BoundaryFlags b = BoundaryFlags(updated_tree, r);
  SparseGrid2D<double> u = SparseGrid2D<double>(updated_tree);
  SparseGrid2D<double> v = SparseGrid2D<double>(updated_tree);
  SparseGrid2D<double> p = SparseGrid2D<double>(updated_tree);
  u[{ 4, 2, u.tree.maxDepth }] = 2.;
  v[{ 4, 2, v.tree.maxDepth }] = 2.;
  // ASSERT(SparseVector::dot(u, v, b) == 4, "<a,b> != 4");
  tree_broadcast(SparseVector::axpy, b, static_cast<uint16_t>(BoundaryType::P_Inside), p, 3., u, v);

  // auto A = SparseMatrixOperator();
  ASSERT((p[{ 4, 2, p.tree.maxDepth }] == 8.), "axpy did not succed");
  // tree_broadcast(SparseVector::aAxpy, b, static_cast<uint16_t>(BoundaryType::P_Inside), p, 3., A, u, v);
};

void test_init_system()
{
  auto t = get_test_tree();
  auto updated_tree = DenseTree::build_tree(is_desired_depth, t.maxDepth, t.maxDepth, t.maxDepth, t);
  PDESystem system = PDESystem(Settings::get(), updated_tree);
};

void test_step_system()
{
  auto t = get_test_tree();
  auto updated_tree = DenseTree::build_tree(is_desired_depth, t.maxDepth, t.maxDepth, t.maxDepth, t);
  PDESystem system = PDESystem(Settings::get(), updated_tree);
  CGSolver solver = CGSolver(updated_tree);
  auto data_set = init(updated_tree);
  step(system, solver, 0.);
  write("PDE System", system, data_set);
  save_dataset(data_set);
};
int main()
{
  signal(SIGINT, signalInt);
  signal(SIGTERM, signalInt);
  LOG::Init(LOG::LoggerType::STDOUT);
  Profiler::Init(Profiler::Type::ACCUMULATE);
  auto settings = std::filesystem::path(__FILE__).parent_path() / "settings.txt";
  if (!Settings::loadFromFile(settings))
  {
    LOG::Warning("couldn't parse settings file");
    LOG::Close();
    Profiler::Close();
    return -1;
  }

  test_build_tree();
  test_tree_refinement();
  test_set_cartesian_index();
  test_set_boundary();
  test_set_values();
  test_vector_operations();
  test_init_system();
  test_step_system();

  LOG::Close();
  Profiler::Close();
};
