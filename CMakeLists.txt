cmake_minimum_required(VERSION 3.15)
project(ELBM_Solver VERSION 1.0 LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Compiler flags for optimization
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -Wall -Wextra")

# Default to Release build
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Include directories
include_directories(${PROJECT_SOURCE_DIR}/include)

# Create output directories at build time
file(MAKE_DIRECTORY ${CMAKE_SOURCE_DIR}/output/analytical_validation)
file(MAKE_DIRECTORY ${CMAKE_SOURCE_DIR}/output/channel_cylinder)
file(MAKE_DIRECTORY ${CMAKE_SOURCE_DIR}/output/twophase)
file(MAKE_DIRECTORY ${CMAKE_SOURCE_DIR}/output/rectangular_pipe)
file(MAKE_DIRECTORY ${CMAKE_SOURCE_DIR}/figures/analytical_validation)
file(MAKE_DIRECTORY ${CMAKE_SOURCE_DIR}/figures/channel_cylinder)
file(MAKE_DIRECTORY ${CMAKE_SOURCE_DIR}/figures/twophase)
file(MAKE_DIRECTORY ${CMAKE_SOURCE_DIR}/figures/geometry_diagrams)

# Find OpenMP (optional)
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    message(STATUS "OpenMP found, enabling parallel support")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

# Main executable (header-only library, no separate core library needed)
add_executable(elbm_solver src/main.cpp)
if(OpenMP_CXX_FOUND)
    target_link_libraries(elbm_solver OpenMP::OpenMP_CXX)
endif()

# Basic flow validation (Couette, Poiseuille, Taylor-Green)
add_executable(test_basic_flows src/test_basic_flows.cpp)
if(OpenMP_CXX_FOUND)
    target_link_libraries(test_basic_flows OpenMP::OpenMP_CXX)
endif()

# Benchmark tests (cavity, cylinder)
add_executable(test_benchmark src/test_benchmark.cpp)
if(OpenMP_CXX_FOUND)
    target_link_libraries(test_benchmark OpenMP::OpenMP_CXX)
endif()

# Extended analytical tests (Womersley, Hagen-Poiseuille, Stokes, Kolmogorov)
add_executable(test_extended_analytical src/test_extended_analytical.cpp)
if(OpenMP_CXX_FOUND)
    target_link_libraries(test_extended_analytical OpenMP::OpenMP_CXX)
endif()

# Two-phase: H-theorem validation (BGK vs ELBM)
add_executable(test_twophase src/test_twophase.cpp)
if(OpenMP_CXX_FOUND)
    target_link_libraries(test_twophase OpenMP::OpenMP_CXX)
endif()

# Color-Gradient two-phase flow with ELBM
add_executable(test_colorgradient src/test_colorgradient.cpp)
if(OpenMP_CXX_FOUND)
    target_link_libraries(test_colorgradient OpenMP::OpenMP_CXX)
endif()

# Droplet collision test with Color-Gradient ELBM
add_executable(test_droplet_collision src/test_droplet_collision.cpp)
if(OpenMP_CXX_FOUND)
    target_link_libraries(test_droplet_collision OpenMP::OpenMP_CXX)
endif()

# Benchmark suite comparing LBGK vs ELBM
add_executable(test_benchmark_suite src/test_benchmark_suite.cpp)
if(OpenMP_CXX_FOUND)
    target_link_libraries(test_benchmark_suite OpenMP::OpenMP_CXX)
endif()


# Testing (optional - can be added later with Catch2)
option(BUILD_TESTS "Build test suite" OFF)
if(BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

# Installation
install(TARGETS elbm_solver DESTINATION bin)
install(DIRECTORY include/ DESTINATION include)


# Active Nematic spontaneous flow test
add_executable(test_active_nematic src/test_active_nematic.cpp)
if(OpenMP_CXX_FOUND)
    target_link_libraries(test_active_nematic OpenMP::OpenMP_CXX)
endif()

# Active Particle Swarm test
add_executable(test_active_swarm src/test_active_swarm.cpp src/active_particles.cpp)
if(OpenMP_CXX_FOUND)
    target_link_libraries(test_active_swarm OpenMP::OpenMP_CXX)
endif()

