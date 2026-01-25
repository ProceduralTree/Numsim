
#include "output/vtk_tree.h"
#include <cstddef>
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

uint16_t intersects_boundary(size_t index, uint16_t maxDepth, uint16_t depth, size_t nx, size_t ny)
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
  bool on_x_boundary = (x <= nx && nx <= x + local_cell_size);
  bool on_y_boundary = (y <= ny && ny <= y + local_cell_size);
  bool in_xy = x <= nx && y <= ny;
  if ((on_x_boundary || on_y_boundary) && in_xy)
    return 1;
  return 0;
};
void test_build_tree()
{
  size_t nx = 50;
  size_t ny = 50;
  size_t x_power = std::bit_width<size_t>(nx - 1);
  size_t y_power = std::bit_width<size_t>(ny - 1);
  const uint16_t maxDepth = std::max(x_power, y_power);
  auto t = DenseTree::build_tree(intersects_boundary, maxDepth, nx, ny);
  t.print();
  write_tree(t);
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
  test_build_tree();
};
