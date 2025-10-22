add_library(grid STATIC grid.cpp indexing.cpp)
target_include_directories(grid PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
