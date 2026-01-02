cmake_minimum_required(VERSION 3.0)

# set mode (Debug or Release)
set(CMAKE_BUILD_TYPE Release)

# set C++ standard
set(CMAKE_CXX_STANDARD 17)

find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

project(HybZonoMPC)
include_directories(include 
    ../mi-solver/include
    ../extern/eigen 
    ../qp-custom/cpp/include
    ../extern/zonocpp/include)

# Build library
add_library(HybZonoMPC SHARED
    src/MPC_QP.cpp
    ../qp-custom/cpp/src/QP_primal_dual_mpc.cpp
)

# Install library
install(TARGETS HybZonoMPC
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    RUNTIME DESTINATION bin
)

# Install headers
install(DIRECTORY include/
    DESTINATION include
    FILES_MATCHING PATTERN "*.hpp"
)

project(MPC_MIQP)
include_directories(include 
    ../mi-solver/include
    ../extern/eigen 
    ../qp-custom/cpp/include
    ../extern/zonocpp/include)
add_executable(DblIntegratorObsAvoid 
    examples/double_integrator_obs_avoid_map_change.cpp 
    src/MPC_MIQP.cpp
    src/MPC_QP.cpp
    ../qp-custom/cpp/src/QP_primal_dual_mpc.cpp
    ../qp-custom/cpp/src/QP_primal_dual.cpp
    ../mi-solver/src/Data.cpp 
    ../mi-solver/src/LeavesQueue.cpp
    ../mi-solver/src/Node.cpp 
    ../mi-solver/src/MI_MPC.cpp 
    ../mi-solver/src/DataVolatile.cpp 
    ../mi-solver/src/VerboseOutputs.cpp
)

# find_package(Threads REQUIRED)
# target_link_libraries(DblIntegratorObsAvoid PRIVATE Threads::Threads)

# project(MPC_MIQP)
# include_directories(include 
#     ../mi-solver/include
#     ../extern/eigen 
#     ../qp-custom/cpp/include)
# add_executable(BicycleObsAvoid 
#     examples/bicycle_example_obs_avoid_cost_map.cpp 
#     src/MPC_MIQP.cpp
#     src/MPC_QP.cpp
#     src/HybZono.cpp 
#     ../qp-custom/cpp/src/QP_primal_dual_mpc.cpp
#     ../qp-custom/cpp/src/QP_primal_dual.cpp
#     ../mi-solver/src/Data.cpp 
#     ../mi-solver/src/LeavesQueue.cpp
#     ../mi-solver/src/Node.cpp 
#     ../mi-solver/src/MI_MPC.cpp 
#     ../mi-solver/src/DataVolatile.cpp 
#     ../mi-solver/src/VerboseOutputs.cpp
#     src/ConZono.cpp)

# find_package(Threads REQUIRED)
# target_link_libraries(BicycleObsAvoid PRIVATE Threads::Threads)
