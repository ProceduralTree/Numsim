#include "utils/broadcast.h"
#include <chrono>
#include <cstdint>
#include <functional>
#include <grid/grid.h>
#include <iostream>
#include <memory>
#include <mpi.h>
#include <output/vtk.h>
#include <pde/system.h>
#include <vector>
#define ASSERT(condition, message)                               \
  do                                                             \
  {                                                              \
    if (!(condition))                                            \
    {                                                            \
      std::cerr << "Assertion failed: " << message << std::endl; \
      assert(condition);                                         \
    }                                                            \
  } while (0)

void set_one(Index I, Offset O, Grid2D& array)
{

  array[I] += O.x;
  array[I] += O.y;
};

void test_boundary(PDESystem& system)
{
  broadcast_boundary(set_one, system.p.boundary, system.p);

  // broadcast_x_boundary(
  //   [&](PDESystem& s, Index I, Offset o) {
  //     s.v[I + o] = o.y;
  //   },
  //   system, system.v);
  // std::cout << system.v;
  // std::cout << "\n";
  // broadcast_boundary(
  //   [&](PDESystem& s, Index I, Offset o) {
  //     s.p[I + o] = o.x + o.y;
  //   },
  //   system, system.p);
  std::cout << system.p;
  std::cout << "\n";
}

void test_index()
{
  Index I = { 1, 1 };
  Index Ipx = I + Ix;
  ASSERT(Ipx.x == 2 && Ipx.y == 1, "Plus Failed I.x=" << Ipx.x << " I.y=" << Ipx.y);
  Index Imx = I - Ix;
  ASSERT(Imx.x == 0 && Imx.y == 1, "Minus Failed I.x=" << Imx.x << " I.y=" << Imx.y);
  ASSERT((-5 * Ix).x == -5 && (-5 * Ix).y == 0, "Invert Failed I.x=" << (-5 * Ix).x << " I.y=" << (-5 * Ix).y);
  ASSERT((-Ix).x == -1 && (-Ix).y == 0, "Invert Failed I.x=" << (-Ix).x << " I.y=" << (-Ix).y);
}

void fill_inorder(Grid2D& grid)
{
  for (uint64_t i = 0; i < grid.size_x * grid.size_y; i++)
  {
    grid[i] = static_cast<double>(i);
  }
};
void benchmark_laplace(
  int N, std::function<void(const Grid2D&, Grid2D&)> algorithm)
{
  auto grid = std::make_unique<Grid2D>(N, N);
  auto out = std::make_unique<Grid2D>(N, N);
  fill_inorder(*grid);
  std::cout << "\tBegin Timing with a " << N << "x" << N << " Grid"
            << std::endl;
  auto start = std::chrono::high_resolution_clock::now();
  for (int i = 0; i < 1000; i++)
  {
    algorithm(*grid, *out);
    std::swap(grid, out);
  }
  // std::cout << *out;
  auto end = std::chrono::high_resolution_clock::now();

  std::cout << "\tDuration: "
            << std::chrono::duration<double, std::milli>(end - start).count()
            << " ms\n";
}

auto main(int argc, char* argv[]) -> int
{
  test_index();
  MPI_Init(&argc, &argv);
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  // std::cout << "With Z-Order Interleaving" << std::endl;
  // benchmark_laplace(N, laplace_cartesian);
  //
  PDESystem test_system = PDESystem(100., 1e-4, 50, 50, 0.02, 0.02, { 0, 0 }, { 1., 0. }, { 0, 0 }, { 0, 0 });
  test_system.settings.loadFromFile("");

  print_pde_system(test_system);

  // test_boundary(test_system);

  for (int i = 0; i < 300; i++)
  {
    std::cout << "\r[Iteration]: " << i << "\t"
              << std::flush;
    step(test_system, i);
    if (i % 10 == 0)
    {
      std::cout << "Pressure \n"
                << test_system.p << std::endl;
      std::cout << "Velocity \n"
                << test_system.u << std::endl;
    }
    write_vtk(test_system, static_cast<double>(i));
  }

  std::cout << "Hello from Rank " << rank << " of " << size << std::endl;

  MPI_Finalize();
  return 0;
}
