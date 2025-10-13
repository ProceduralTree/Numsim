#include <iostream>
#include <mpi.h>

auto main(int argc, char *argv[]) -> int {
  MPI_INIT(NULL , NULL);
  std::cout << "Hello World" << std::endl;
  return 0;
}
