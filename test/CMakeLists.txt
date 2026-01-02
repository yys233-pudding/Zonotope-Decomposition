cmake_minimum_required(VERSION 3.15)
project(example_operations)

# Using add_subdirectory
add_subdirectory(${CMAKE_SOURCE_DIR}/.. ZonoOpt)

add_executable(example_operations
    src/example_operations.cpp
)
target_link_libraries(example_operations PRIVATE ZonoOpt)

add_executable(decomposition_comparison
    src/decomposition_comparison.cpp
)
target_link_libraries(decomposition_comparison PRIVATE ZonoOpt)