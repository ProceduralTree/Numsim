#include "grid/grid.h"
#include "indexing.h"
#include <chrono>
#include <cstdint>
#include <functional>
#include <iostream>
#include <memory>
#include <mpi.h>
#include <vector>
template <Indexing I> void fill_inorder(Grid2D<I> &grid) {
    for (uint64_t i = 0; i < grid.size_x * grid.size_y; i++) {
        grid[i] = static_cast<double>(i);
    }
};
template <Indexing I>
void benchmark_laplace(
    int N, std::function<void(const Grid2D<I> &, Grid2D<I> &)> algorithm) {
    auto grid = std::make_unique<Grid2D<I>>(N, N);
    auto out = std::make_unique<Grid2D<I>>(N, N);
    fill_inorder(*grid);
    std::cout << "\tBegin Timing with a " << N << "x" << N << " Grid"
              << std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < 1000; i++) {
        algorithm(*grid, *out);
        std::swap(grid, out);
    }
    // std::cout << *out;
    auto end = std::chrono::high_resolution_clock::now();

    std::cout << "\tDuration: "
              << std::chrono::duration<double, std::milli>(end - start).count()
              << " ms\n";
}

auto main(int argc, char *argv[]) -> int {
    MPI_Init(&argc, &argv);
    std::vector<int> init;
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    const int N = 1 << 11;
    std::cout << "With Z-Order Index" << std::endl;
    benchmark_laplace<Indexing::ZOrder>(N, laplace<Indexing::ZOrder>);
    std::cout << "With Cartesian Index" << std::endl;
    benchmark_laplace<Indexing::Cartesian>(N, laplace<Indexing::Cartesian>);
    // std::cout << "With Z-Order Interleaving" << std::endl;
    // benchmark_laplace(N, laplace_cartesian);

    std::cout << "Hello from Rank " << rank << " of " << size << std::endl;

    MPI_Finalize();
    return 0;
}
