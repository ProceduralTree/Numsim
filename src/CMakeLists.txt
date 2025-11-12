
find_package(VTK REQUIRED COMPONENTS
   CommonCore
    CommonDataModel
    IOXML
)
find_package(MPI REQUIRED)

file(GLOB_RECURSE SRC_FILES *.cpp)
file(GLOB_RECURSE HEADER_FILES *.h)

add_executable(numsim main.cc
  ${SRC_FILES} ${HEADER_FILES}
)

add_executable(numsim_profiling main.cc
  ${SRC_FILES} ${HEADER_FILES}
)
target_link_libraries(numsim PUBLIC MPI::MPI_CXX ${VTK_LIBRARIES})
target_link_libraries(numsim PRIVATE  ${VTK_LIBRARIES})
target_link_libraries(numsim_profiling PUBLIC MPI::MPI_CXX ${VTK_LIBRARIES})
target_link_libraries(numsim_profiling PRIVATE  ${VTK_LIBRARIES})

# MPI

# Profiling executable

target_compile_options(numsim_profiling PRIVATE -pg)
target_link_options(numsim_profiling PRIVATE -pg)




vtk_module_autoinit(
    TARGETS numsim
    MODULES ${VTK_LIBRARIES}
)

target_include_directories(numsim PRIVATE ${VTK_INCLUDE_DIRS})
target_include_directories(numsim PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(numsim_profiling PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

install(TARGETS numsim RUNTIME DESTINATION ${PROJECT_SOURCE_DIR}/build)
