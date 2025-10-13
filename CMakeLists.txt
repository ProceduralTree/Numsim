cmake_minimum_required(VERSION 3.10)
project (main)

set(CMAKE_CXX_STANDARD 14)
# -fopenmp flag (enables OpenMP)
set(GCC_COVERAGE_COMPILE_FLAGS "-Wall -pedantic -lm -O3 -funroll-loops -fopenmp -g")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}")

add_executable(mpi src/main.cc)

# MPI
find_package(MPI REQUIRED)
target_link_libraries(mpi PUBLIC MPI::MPI_CXX)
