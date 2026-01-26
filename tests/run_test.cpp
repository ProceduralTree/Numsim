
#include "grid/boundary.h"
#include "grid/sparsegrid.h"
#include "output/vtk_tree.h"
#include "utils/Logger.h"
#include "utils/index.h"
#include "utils/profiler.h"
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
  auto end = Index { 50, 50, 0 };
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
  write_field("Indices", t, t._depths, data_set);
  save_dataset(data_set);
};

void test_tree_refinement()
{
  auto t = get_test_tree();

  auto updated_tree = DenseTree::build_tree(is_desired_depth, t.maxDepth, t.maxDepth, 4, t);
  updated_tree.print();
  auto data_set = init(updated_tree);
  write_field("Updated Indices", updated_tree, updated_tree._indices, data_set);
  save_dataset(data_set);
};

void test_set_boundary()
{
  auto t = get_test_tree();
  Range r = get_test_range();
  BoundaryFlags b = BoundaryFlags(t, r);
  auto data_set = init(t);
  write_field("BoundaryFlags", b.tree, b.boundary_flags, data_set);
  save_dataset(data_set);
};

// void set_one(Index I, Offset O, Grid2D& array)
//{
//
// array[I] += O.x;
// array[I] += O.y;
//};
//
// void test_boundary(PDESystem& system)
//{
//// TODO add proper boundary tests
// broadcast_boundary(set_one, system.p.boundary, system.p);
// ASSERT(system.p[system.p.end + Ix] == 1, "system boundary was " << system.p[system.p.end + Ix])
//
// std::cout << system.p;
//}
//
// void test_index()
//{
// Index I = { 1, 1 };
// Index Ipx = I + Ix;
// ASSERT(Ipx.x == 2 && Ipx.y == 1, "Plus Failed I.x=" << Ipx.x << " I.y=" << Ipx.y);
// Index Imx = I - Ix;
// ASSERT(Imx.x == 0 && Imx.y == 1, "Minus Failed I.x=" << Imx.x << " I.y=" << Imx.y);
// ASSERT((-5 * Ix).x == -5 && (-5 * Ix).y == 0, "Invert Failed I.x=" << (-5 * Ix).x << " I.y=" << (-5 * Ix).y);
// ASSERT((-Ix).x == -1 && (-Ix).y == 0, "Invert Failed I.x=" << (-Ix).x << " I.y=" << (-Ix).y);
//}

int main()
{
  signal(SIGINT, signalInt);
  signal(SIGTERM, signalInt);
  LOG::Init(LOG::LoggerType::STDOUT);
  Profiler::Init(Profiler::Type::ACCUMULATE);

  test_build_tree();
  test_tree_refinement();
  test_set_boundary();

  LOG::Close();
  Profiler::Close();
};
