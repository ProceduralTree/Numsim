#include "grid/grid.h"
#include <iostream>
#include <mpi.h>
#include <vector>

auto main(int argc, char *argv[]) -> int {
    MPI_Init(&argc, &argv);
    std::vector<int> init;
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    Grid2D grid = Grid2D(16, 16);
    Grid2D out = Grid2D(16, 16);
    grid[1, 1] = 1;
    grid[2, 2] = 69.420;
    grid[7, 3] = 2.132;
    grid[13, 8] = 5.3;

    if (rank == 0) {
        std::cout << "Grid: \t\n" << grid << std::endl;
        laplace(grid, out);
        std::cout << "Grid after Laplace: \t\n" << out << std::endl;
    }

    std::cout << "Hello from Rank " << rank << " of " << size << std::endl;
    hello();

    MPI_Finalize();
    return 0;
}
