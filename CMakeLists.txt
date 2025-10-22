cmake_minimum_required(VERSION 3.10)
project (Numsim)


set(CMAKE_CXX_STANDARD 23)
# -fopenmp flag (enables OpenMP)
set(GCC_COVERAGE_COMPILE_FLAGS "-march=native -Wall -pedantic -lm -O3  -funroll-loops -ftree-parallelize-loops=8  -fopenmp -g   -ftree-vectorize -fopt-info-vec")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}")

add_subdirectory(src)

