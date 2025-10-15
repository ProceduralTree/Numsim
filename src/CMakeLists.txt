
add_subdirectory(grid)
add_executable(mpi main.cc)

# MPI
find_package(MPI REQUIRED)
target_link_libraries(mpi PUBLIC MPI::MPI_CXX)
target_link_libraries(mpi PRIVATE grid)

target_include_directories(mpi PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
