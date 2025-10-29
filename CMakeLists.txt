cmake_minimum_required(VERSION 3.10)
project (Numsim)


set(CMAKE_CXX_STANDARD 23)
# -fopenmp flag (enables OpenMP)
set(GCC_COVERAGE_COMPILE_FLAGS "-march=native -Wall -pedantic -lm -O0  -funroll-loops -ftree-parallelize-loops=8  -fopenmp -g   -ftree-vectorize ")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}")


#set(CMAKE_C_COMPILER "/usr/lib64/openmpi/bin/mpicc")
#set(CMAKE_CC_COMPILER "/usr/lib64/openmpi/bin/mpicxx")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

add_subdirectory(src)
#add_subdirectory(external)

