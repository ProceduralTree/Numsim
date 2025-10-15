#include "grid/grid.h"
#include <iostream>
#include <mpi.h>

auto main(int argc, char *argv[]) -> int {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    std::cout << "Hello from Rank " << rank << " of " << size << std::endl;
    hello();

    MPI_Finalize();
    return 0;
}
