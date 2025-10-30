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
  MPI_Init(&argc, &argv);
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  // std::cout << "With Z-Order Interleaving" << std::endl;
  // benchmark_laplace(N, laplace_cartesian);
  //
  PDESystem test_system = PDESystem(1., 1e-4, 100, 100, 0.01, 0.01, { 0, 0 }, { 0, 0.1 }, { 0, 0 }, { 0, 0 });
  test_system.settings.loadFromFile("");
  print_pde_system(test_system);

  for (int i = 0; i < 10; i++)
  {
    std::cout << "[Iteration]: " << i << "\t\r" << std::flush;
    timestep(test_system);
    write_vtk(test_system, static_cast<double>(i));
  }

  std::cout << "Hello from Rank " << rank << " of " << size << std::endl;

  MPI_Finalize();
  return 0;
}
