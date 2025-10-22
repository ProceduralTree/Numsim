
add_subdirectory(grid)
add_executable(mpi main.cc)

# MPI
find_package(MPI REQUIRED)
target_link_libraries(mpi PUBLIC MPI::MPI_CXX)
target_link_libraries(mpi PRIVATE grid)

target_include_directories(mpi PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# Profiling executable
add_executable(mpi_profiling main.cc)
target_compile_options(mpi_profiling PRIVATE -pg)
target_link_libraries(mpi_profiling PUBLIC MPI::MPI_CXX PRIVATE grid)
target_link_options(mpi_profiling PRIVATE -pg)

