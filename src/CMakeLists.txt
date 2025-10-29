
find_package(VTK REQUIRED COMPONENTS
   CommonCore
    CommonDataModel
    IOXML
)
find_package(MPI REQUIRED)

file(GLOB_RECURSE SRC_FILES *.cpp)
file(GLOB_RECURSE HEADER_FILES *.h)

add_executable(simnum main.cc
  ${SRC_FILES} ${HEADER_FILES}
)

add_executable(simnum_profiling main.cc
  ${SRC_FILES} ${HEADER_FILES}
)
target_link_libraries(simnum PUBLIC MPI::MPI_CXX ${VTK_LIBRARIES})
target_link_libraries(simnum PRIVATE  ${VTK_LIBRARIES})
target_link_libraries(simnum_profiling PUBLIC MPI::MPI_CXX ${VTK_LIBRARIES})
target_link_libraries(simnum_profiling PRIVATE  ${VTK_LIBRARIES})

# MPI

# Profiling executable

target_compile_options(simnum_profiling PRIVATE -pg)
target_link_options(simnum_profiling PRIVATE -pg)




vtk_module_autoinit(
    TARGETS simnum
    MODULES ${VTK_LIBRARIES}
)

target_include_directories(simnum PRIVATE ${VTK_INCLUDE_DIRS})
target_include_directories(simnum PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(simnum_profiling PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
