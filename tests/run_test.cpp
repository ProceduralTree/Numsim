
#include "grid/boundary.h"
#include "grid/sparsegrid.h"
#include "output/vtk_tree.h"
#include "utils/Logger.h"
#include "utils/index.h"
#include "utils/profiler.h"
#include <cassert>
#include <csignal>
#include <cstddef>
#include <cstdint>
#include <grid/densetree.h>
#include <iostream>
#include <mpi.h>

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
  MPI_Barrier(MPI_COMM_WORLD);
  exit(sig);
};

uint16_t intersects_range(size_t index, uint16_t depth, Range r, uint16_t maxDepth)
{
  // size_t x_power = std::bit_width<size_t>(nx - 1);
  // size_t y_power = std::bit_width<size_t>(ny - 1);
  // const size_t size = std::max(x_power, y_power);
  if (depth >= maxDepth)
    return 0;

  size_t local_cell_size = 1 << (maxDepth - depth);
  auto [x, y, _] = ZorderToIndex(index);
  // std::cout << "X:" << x << "Y:" << y << std::endl;
  // std::cout << "local size:" << local_cell_size << std::endl;
  bool on_top_boundary = (x <= r.end.x && r.end.x <= x + local_cell_size);
  bool on_bottom_boundary = (x <= r.begin.x && r.begin.x <= x + local_cell_size);
  bool on_left_boundary = (y <= r.begin.y && r.begin.y <= y + local_cell_size);
  bool on_right_boundary = (y <= r.end.y && r.end.y <= y + local_cell_size);
  bool in_range = x <= r.end.x && y <= r.end.y;
  if ((on_top_boundary || on_bottom_boundary || on_left_boundary || on_right_boundary) && in_range)
    return 1;
  return 0;
};
uint16_t is_desired_depth(size_t index, uint16_t depth, uint16_t maxDepth, uint16_t desired, const DenseTree::DenseTree& old_tree)
{
  if (depth <= desired && depth <= maxDepth)
  {
    return 1;
  }
  size_t old_index = DenseTree::get_dense_index(old_tree, { index, maxDepth });
  auto [_, old_depth] = old_tree._index_cache[old_index];
  if (depth <= old_depth && depth <= maxDepth)
  {
    return 1;
  }
  return 0;
};
Range get_test_range()
{
  auto begin = Index { 1, 1, 0 };
  auto end = Index { 5, 5, 0 };
  return { begin, end };
}

DenseTree::DenseTree get_test_tree()
{
  Range r = get_test_range();

  size_t x_power = std::bit_width<size_t>(r.end.x - 1);
  size_t y_power = std::bit_width<size_t>(r.end.y - 1);
  const uint16_t maxDepth = std::max(x_power, y_power);
  auto t = DenseTree::build_tree(intersects_range, maxDepth, r, maxDepth);
  return t;
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
  auto sparse_grid = SparseGrid2D<double>(t);
  auto data_set = init(updated_tree);
  for (uint16_t i = 0; i < 1ULL << t.maxDepth; i++)
  {
    sparse_grid[{ i, i, t.maxDepth }] = 1. * i + 1.;
  }
  for (uint16_t i = 0; i < 10; i++)
  {
    assert((sparse_grid[{ i, i, t.maxDepth }] == 1. * i && "Did not set value at expected point"));
  }
  write_field("Grid", sparse_grid.tree, sparse_grid._data, data_set);
  // ASSERT(sparse_grid[{ 5, 5, t.maxDepth }] == 1., "Did not set value at expected point");
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
  auto updated_tree = DenseTree::build_tree(is_desired_depth, t.maxDepth, t.maxDepth, 4, t);
  Range r = get_test_range();
  BoundaryFlags b = BoundaryFlags(updated_tree, r);
  SparseGrid2D<double> ugrid = SparseGrid2D<double>(updated_tree);
  SparseGrid2D<double> vgrid = SparseGrid2D<double>(updated_tree);
  SparseGrid2D<double> pgrid = SparseGrid2D<double>(updated_tree);
  uint16_t u_boundary = static_cast<uint16_t>(BoundaryType::U_BOTTOM)
    | static_cast<uint16_t>(BoundaryType::U_TOP)
    | static_cast<uint16_t>(BoundaryType::U_LEFT)
    | static_cast<uint16_t>(BoundaryType::U_RIGHT);
  uint16_t v_boundary = static_cast<uint16_t>(BoundaryType::V_BOTTOM)
    | static_cast<uint16_t>(BoundaryType::V_TOP)
    | static_cast<uint16_t>(BoundaryType::V_LEFT)
    | static_cast<uint16_t>(BoundaryType::V_RIGHT);
  uint16_t p_boundary = static_cast<uint16_t>(BoundaryType::P_BOTTOM)
    | static_cast<uint16_t>(BoundaryType::P_TOP)
    | static_cast<uint16_t>(BoundaryType::P_LEFT)
    | static_cast<uint16_t>(BoundaryType::P_RIGHT);
  auto data_set
    = init(updated_tree);
  tree_broadcast(_set, b, u_boundary, ugrid, 1.);
  tree_broadcast(_set, b, v_boundary, vgrid, 1.);
  tree_broadcast(_set, b, p_boundary, pgrid, 1.);
  write_field("U boundary", b.tree, ugrid._data, data_set);
  write_field("V boundary", b.tree, vgrid._data, data_set);
  write_field("P boundary", b.tree, pgrid._data, data_set);
  save_dataset(data_set);
};

int main()
{
  signal(SIGINT, signalInt);
  signal(SIGTERM, signalInt);
  LOG::Init(LOG::LoggerType::STDOUT);
  Profiler::Init(Profiler::Type::ACCUMULATE);

  test_build_tree();
  test_tree_refinement();
  test_set_cartesian_index();
  test_set_boundary();
  test_set_values();

  LOG::Close();
  Profiler::Close();
};
