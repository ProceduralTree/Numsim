find_package(VTK REQUIRED COMPONENTS
    CommonCore
    CommonDataModel
    IOXML
)
find_package(MPI REQUIRED)

file(GLOB_RECURSE SRC_FILES
     CONFIGURE_DEPENDS
     ${CMAKE_CURRENT_LIST_DIR}/../src/*.cpp
)

add_executable(tests
    run_test.cpp
    ${SRC_FILES}
)

target_link_libraries(tests
    PRIVATE
        MPI::MPI_CXX
        ${VTK_LIBRARIES}
)

target_compile_definitions(tests
    PRIVATE $<$<CONFIG:Debug>:DEBUG>
)

vtk_module_autoinit(
    TARGETS tests
    MODULES ${VTK_LIBRARIES}
)

target_include_directories(tests
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/../src/
)
