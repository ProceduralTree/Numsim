
#include "grid/adjacencymap.h"
#include "grid/boundary.h"
#include "grid/quadtree.h"
#include "grid/rectangle.h"
#include "grid/sparsegrid.h"
#include "linalg/matrix.h"
#include "output/vtk_tree.h"
#include "pde/system.h"
#include "utils/Logger.h"
#include "utils/broadcast.h"
#include "utils/index.h"
#include "utils/profiler.h"
#include "utils/settings.h"
#include <cassert>
#include <csignal>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <format>
#include <fstream>
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
  auto begin = Index { 2, 2, 0 };
  auto end = Index { 51, 51, 0 };
  return { begin, end };
}

DenseTree::DenseTree get_test_tree()
{

  return DenseTree::from_range(get_test_range());
};

void test_build_tree()
{

  auto t = get_test_tree();
  auto data_set = init(t, false);
  // ASSERT(tree.sizes, message)
  // t.print();
  write_depth("Tree Depth", t, data_set);
  save_dataset(data_set);
};

void test_build_tree_from_image()
{
  QuadTree imageTree("input/test32xp.png");
  // std::ofstream file("outputTest.txt");
  // file << imageTree;
  // file.close();
  // std::cout << imageTree;

  auto t = DenseTree::build_tree([&](auto z, auto d) { return imageTree.hasChildrenP(z, d); }, imageTree.getDepth());

  BoundaryFlags b(t, imageTree);

  auto data_set = init(t, false);
  write_depth("depth", t, data_set);
  write_field("BoundaryFlags", b.tree, b.flags._data, data_set);
  save_dataset(data_set);
}

void test_tree_refinement()
{
  auto t = get_test_tree();

  auto updated_tree = DenseTree::build_tree(is_desired_depth, t.maxDepth, t.maxDepth, t.maxDepth - 1, t);
  // updated_tree.print();
  auto data_set = init(updated_tree, false);
  write_depth("Updated Depth", updated_tree, data_set);
  save_dataset(data_set);
};

void test_set_cartesian_index()
{
  auto t = get_test_tree();

  std::cerr << "MaxDepth:" << t.maxDepth << std::endl;
  std::cerr << "Sizes" << std::endl;
  std::ranges::copy(t._sizes, std::ostream_iterator<size_t>(std::cerr, " "));
  std::cerr << std::endl;
  auto updated_tree = DenseTree::build_tree(is_desired_depth, t.maxDepth, t.maxDepth, t.maxDepth - 1, t);
  std::cerr << "MaxDepth:" << updated_tree.maxDepth << std::endl;
  std::cerr << "Sizes" << std::endl;
  std::ranges::copy(updated_tree._sizes, std::ostream_iterator<size_t>(std::cerr, " "));
  std::cerr << std::endl;
  auto sparse_grid = SparseGrid2D<double>(updated_tree);
  auto data_set = init(updated_tree, false);
  for (uint16_t i = 1; i < get_test_range().end.x; i++)
  {
    sparse_grid[{ i, i, updated_tree.maxDepth }] = 1. * i + 1.;
  }
  for (uint16_t i = 1; i < get_test_range().end.x; i++)
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
  auto data_set = init(updated_tree, false);
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
  auto updated_tree = DenseTree::build_tree(is_desired_depth, t.maxDepth, t.maxDepth, t.maxDepth - 1, t);
  Range r = get_test_range();
  BoundaryFlags b = BoundaryFlags(updated_tree, r);
  SparseGrid2D<double> ugrid = SparseGrid2D<double>(updated_tree);
  SparseGrid2D<double> vgrid = SparseGrid2D<double>(updated_tree);
  SparseGrid2D<double> pgrid = SparseGrid2D<double>(updated_tree);
  uint16_t u_boundary = static_cast<uint16_t>(BoundaryType::U_BOUNDARY);
  uint16_t v_boundary = static_cast<uint16_t>(BoundaryType::V_BOUNDARY);
  uint16_t p_boundary = static_cast<uint16_t>(BoundaryType::P_BOUNDARY);
  auto data_set = init(updated_tree, false);
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
  auto updated_tree = DenseTree::build_tree(is_desired_depth, t.maxDepth, t.maxDepth, t.maxDepth - 1, t);
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
  ASSERT((p[{ 4, 2, p.tree.maxDepth }] == 8.), (std::format("axpy did not succed p={}", p[{ 4, 2, p.tree.maxDepth }])));
  // tree_broadcast(SparseVector::aAxpy, b, static_cast<uint16_t>(BoundaryType::P_Inside), p, 3., A, u, v);
};

void test_init_system()
{
  auto r = Range { Index { 1, 1, 0 }, Index { static_cast<uint16_t>(Settings::get().nCells[0] + 1), static_cast<uint16_t>(Settings::get().nCells[1] + 1), 0 } };
  auto t = DenseTree::from_range(r);

  auto tree = DenseTree::build_tree(is_desired_depth, t.maxDepth, t.maxDepth, t.maxDepth, t);
  auto flags = BoundaryFlags(tree, r);
  PDESystem system = PDESystem(Settings::get(), flags);
  CGSolver solver = CGSolver(tree);
  auto data_set = init(tree, false);
  write_field("boundary", system.boundary.tree, system.boundary.flags._data, data_set);
  write_depth("Depth", system.boundary.tree, data_set);
  save_dataset(data_set);
};

void test_step_system()
{
  auto r = Range { Index { 1, 1, 0 }, Index { static_cast<uint16_t>(Settings::get().nCells[0] + 1), static_cast<uint16_t>(Settings::get().nCells[1] + 1), 0 } };
  auto t = DenseTree::from_range(r);

  auto tree = DenseTree::build_tree(is_desired_depth, t.maxDepth, t.maxDepth, t.maxDepth, t);
  auto flags = BoundaryFlags(tree, r);
  PDESystem system = PDESystem(Settings::get(), flags);
  // CGSolver solver = CGSolver(tree);
  Jacoby solver = Jacoby(tree);
  auto data_set = init(tree, false);
  write_field("Boundary Data", system.boundary.tree, system.boundary.flags._data, data_set);
  step(system, solver, 0.);
  write_field("U raw data", system.boundary.tree, system.u._data, data_set);
  write_field("V raw data", system.boundary.tree, system.v._data, data_set);
  write_field("P raw data", system.boundary.tree, system.p._data, data_set);
  write_field("F raw data", system.boundary.tree, system.F._data, data_set);
  write_field("G raw data", system.boundary.tree, system.G._data, data_set);
  write_field("RHS raw data", system.boundary.tree, system.rhs._data, data_set);
  write_field("Residual", system.boundary.tree, solver.residual._data, data_set);
  // write_field("TMP", system.boundary.tree, solver.tmp._data, data_set);
  //  write_field("Search Direction", system.boundary.tree, solver.search_direction._data, data_set);
  save_dataset(data_set);
};

void test_set_pressure_boundary()
{
  auto t = get_test_tree();
  auto updated_tree = DenseTree::build_tree(is_desired_depth, t.maxDepth, t.maxDepth, t.maxDepth, t);
  auto flags = BoundaryFlags(updated_tree, get_test_range());
  PDESystem system = PDESystem(Settings::get(), flags);
  CGSolver solver = CGSolver(updated_tree);
  auto data_set = init(updated_tree, false);
  tree_broadcast(_set, system.boundary, static_cast<double>(BoundaryType::P_Inside), system.p, 1.);
  broadcast_boundary(copy_with_offset, system.boundary, static_cast<uint16_t>(BoundaryType::P_BOUNDARY), system.p);
  write_field("Pressure", flags.tree, system.p._data, data_set);
  // write("PDE System", system, data_set);
  write_field("BoundaryFlags", flags.tree, flags.flags._data, data_set);
  save_dataset(data_set);
};

void test_mat_mult()
{
  auto t = get_test_tree();
  auto updated_tree = DenseTree::build_tree(is_desired_depth, t.maxDepth, t.maxDepth, t.maxDepth, t);
  auto flags = BoundaryFlags(updated_tree, get_test_range());
  PDESystem system = PDESystem(Settings::get(), flags);
  CGSolver solver = CGSolver(updated_tree);
  auto data_set = init(updated_tree, false);
  system.p[{ 25, 25, updated_tree.maxDepth }] = 1.;
  auto A = SparseMatrixOperator(system.h, system.adjacency_map);
  tree_broadcast(SparseVector::aAxpy, system.boundary, static_cast<uint16_t>(BoundaryType::P_Inside), system.u, 1., A, system.p, system.p);
  write_field("Mat Mult", flags.tree, system.u._data, data_set);
  save_dataset(data_set);
};

void test_print()
{
  auto r = Range { Index { 0, 0, 0 }, Index { 40, 40, 0 } };
  auto t = DenseTree::from_range(r);
  auto data_set = init(t, false);
  write_data(set_index, "Indices", 1, t, data_set, false, t);
  write_depth("Depth", t, data_set);
  save_dataset(data_set);
}

void test_build_better_tree()
{
  auto t = get_test_tree();
  auto tree = dilate(t);
  auto flags = BoundaryFlags(tree, get_test_range());

  auto data_set = init(tree, false);
  write_field("Boundary Flags N 0", flags.tree, flags.flags._data, flags.tree.maxDepth, data_set);
  write_field("Boundary Flags N-1", flags.tree, flags.flags._data, flags.tree.maxDepth - 1, data_set);
  write_field("Boundary Flags N-2", flags.tree, flags.flags._data, flags.tree.maxDepth - 2, data_set);
  write_field("Boundary Flags N-3", flags.tree, flags.flags._data, flags.tree.maxDepth - 3, data_set);
  write_depth("Original Boundary", t, data_set);
  write_depth("Boundary", tree, data_set);
  save_dataset(data_set);
}

void test_mipmap()
{
  auto t = get_test_tree();
  auto tree = dilate(t);
  auto flags = BoundaryFlags(tree, get_test_range());
  PDESystem system = PDESystem(Settings::get(), flags);
  system.p[{ 4, 3, tree.maxDepth }] = 5.;
  mipmap(_mean<double>, system.p);

  auto data_set = init(tree, false);
  write_field("P(N-0)  raw data", system.boundary.tree, system.p._data, system.boundary.tree.maxDepth - 0, data_set);
  write_field("P(N-1)  raw data", system.boundary.tree, system.p._data, system.boundary.tree.maxDepth - 1, data_set);
  write_field("P(N-2)  raw data", system.boundary.tree, system.p._data, system.boundary.tree.maxDepth - 2, data_set);
  save_dataset(data_set);
};

void test_heterogenious_aAxpy()
{
  auto t = get_test_tree();
  auto tree = dilate(t);
  auto flags = BoundaryFlags(tree, get_test_range());
  PDESystem system = PDESystem(Settings::get(), flags);
  system.p[{ 4, 3, tree.maxDepth }] = 5.;
  auto A = SparseMatrixOperator(system.h, system.adjacency_map);
  mipmap(_mean<double>, system.p);
  tree_broadcast(SparseVector::aAxpy, system.boundary, static_cast<uint16_t>(BoundaryType::P_Inside), system.u, 1., A, system.p, system.p);
  mipmap(_mean<double>, system.u);

  auto data_set = init(tree, false);
  write_field("P(N-0)  raw data", system.boundary.tree, system.u._data, system.boundary.tree.maxDepth - 0, data_set);
  write_field("P(N-1)  raw data", system.boundary.tree, system.u._data, system.boundary.tree.maxDepth - 1, data_set);
  write_field("P(N-2)  raw data", system.boundary.tree, system.u._data, system.boundary.tree.maxDepth - 2, data_set);
  save_dataset(data_set);
};

void test_time_step_mipmap()
{
  auto r = Range { Index { 1, 1, 0 }, Index { static_cast<uint16_t>(Settings::get().nCells[0] + 1), static_cast<uint16_t>(Settings::get().nCells[1] + 1), 0 } };
  auto t = DenseTree::from_range(r);

  // auto t = dilate(t0);
  auto flags0 = BoundaryFlags(t, r);
  auto tree = dilate(flags0);
  auto flags = BoundaryFlags(tree, r);
  PDESystem system = PDESystem(Settings::get(), flags);
  // CGSolver solver = CGSolver(tree);
  auto solver = CGSolver(tree);
  for (int i = 0; i < 1; i++)
  {
    step(system, solver, 0.);
  }

  auto data_set = init(tree, false);
  write_depth("Boundary", tree, data_set);
  write_field("Boundary Data", system.boundary.tree, system.boundary.flags._data, data_set);
  write_field("U raw data", system.boundary.tree, system.u._data, data_set);
  write_field("V raw data", system.boundary.tree, system.v._data, data_set);
  write_field("P raw data", system.boundary.tree, system.p._data, data_set);
  write_field("F raw data", system.boundary.tree, system.F._data, data_set);
  write_field("G raw data", system.boundary.tree, system.G._data, data_set);
  write_field("RHS raw data", system.boundary.tree, system.rhs._data, data_set);
  write_field("Residual", system.boundary.tree, solver.residual._data, data_set);
  save_dataset(data_set);
  // auto data_set = init(tree, false);
  // write_field("G raw data", system.boundary.tree, system.G._data, data_set);
  // write_field("F raw data", system.boundary.tree, system.F._data, data_set);
  // write_field("P(N-0)  raw data", system.boundary.tree, system.p._data, system.boundary.tree.maxDepth - 0, data_set);
  // write_field("P(N-1)  raw data", system.boundary.tree, system.p._data, system.boundary.tree.maxDepth - 1, data_set);
  // write_field("P(N-2)  raw data", system.boundary.tree, system.p._data, system.boundary.tree.maxDepth - 2, data_set);
  // write_field("P(N-3)  raw data", system.boundary.tree, system.p._data, system.boundary.tree.maxDepth - 3, data_set);
  // write_field("RHS raw data", system.boundary.tree, system.rhs._data, data_set);
  // write_field("Residual", system.boundary.tree, solver.residual._data, data_set);
  //// write_field("Search Direction", system.boundary.tree, solver.search_direction._data, data_set);
  // write_depth("Boundary", tree, data_set);
  // save_dataset(data_set);
};

void test_adjacency_map()
{
  auto r = Range { Index { 1, 1, 0 }, Index { static_cast<uint16_t>(Settings::get().nCells[0] + 1), static_cast<uint16_t>(Settings::get().nCells[1] + 1), 0 } };
  auto t0 = DenseTree::from_range(r);

  auto t = dilate(t0);
  auto flags0 = BoundaryFlags(t, r);
  auto tree = dilate(flags0);
  auto flags = BoundaryFlags(tree, r);
  auto adjacency_map = AdjMap(flags);
  auto grid = SparseGrid2D<double>(flags.tree);
  Index I = { 10, 10, static_cast<uint16_t>(tree.maxDepth) };
  size_t local_index = DenseTree::get_dense_index(tree, I);
  uint16_t local_depth = tree._index_cache[local_index].depth;
  I.depth = local_depth;
  size_t local_cell_size = 1ULL << (tree.maxDepth - local_depth);

  grid[local_index] = 1.;
  size_t left_index = adjacency_map._left[local_index];
  size_t left = DenseTree::get_dense_index(tree, I - local_cell_size * Ix);
  grid[left_index] = 3.;

  size_t right_index = adjacency_map._right[local_index];
  size_t right = DenseTree::get_dense_index(tree, I - local_cell_size * Ix);
  grid[right_index] = 2.;

  size_t top_index = adjacency_map._top[local_index];
  size_t top = DenseTree::get_dense_index(tree, I + local_cell_size * Iy);
  grid[top_index] = 4.;

  size_t bottom_index = adjacency_map._bottom[local_index];
  size_t bottom = DenseTree::get_dense_index(tree, I - local_cell_size * Iy);
  grid[bottom_index] = 5.;

  auto data_set = init(tree, false);
  write_depth("Depth", tree, data_set);
  write_field("Adjacency setter (N-0)", tree, grid._data, tree.maxDepth - 0, data_set);
  write_field("Adjacency setter (N-1)", tree, grid._data, tree.maxDepth - 1, data_set);
  write_field("Adjacency setter (N-2)", tree, grid._data, tree.maxDepth - 2, data_set);
  write_field("Adjacency setter (N-3)", tree, grid._data, tree.maxDepth - 3, data_set);
  write_field("Adjacency setter (N-4)", tree, grid._data, tree.maxDepth - 4, data_set);
  save_dataset(data_set);

  ASSERT(
    (grid[I - local_cell_size * Iy] == 5.),
    (std::format("Failed setting bottom neighbour using Adjacency map \n set_value={},\n expected value={},\n bottom_index={},\n real_index={},\n local_depth={}", grid[I - local_cell_size * Iy], 5., bottom_index, bottom, local_depth)));
  ASSERT(
    (grid[I + local_cell_size * Iy] == 4.),
    (std::format("Failed setting top neighbour using Adjacency map \n set_value={},\n expected value={},\n bottom_index={},\n real_index={},\n local_depth={}", grid[I + local_cell_size * Iy], 4., top_index, top, local_depth)));
  ASSERT(
    (grid[I + local_cell_size * Ix] == 2.),
    (std::format("Failed setting right neighbour using Adjacency map \n set_value={},\n expected value={},\n bottom_index={},\n real_index={},\n local_depth={}", grid[I + local_cell_size * Ix], 2., left_index, right, local_depth)));
  ASSERT(
    (grid[I - local_cell_size * Ix] == 3.),
    (std::format("Failed setting left neighbour using Adjacency map \n set_value={},\n expected value={},\n bottom_index={},\n real_index={},\n local_depth={}", grid[I - local_cell_size * Ix], 3., left_index, left, local_depth)));
}

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

  test_build_tree_from_image();

  // test_build_tree();
  // test_tree_refinement();
  // test_set_cartesian_index();
  // test_set_boundary();
  // test_set_values();
  // // test_vector_operations();
  // test_init_system();
  // test_step_system();
  // test_set_pressure_boundary();
  // test_mat_mult();
  // test_build_better_tree();
  // test_print();

  LOG::Close();
  Profiler::Close();
};
