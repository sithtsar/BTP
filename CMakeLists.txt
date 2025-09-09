cmake_minimum_required(VERSION 3.15)
project(LBM_Poiseuille VERSION 1.0.0 LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set build type if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Compiler flags
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -Wall -Wextra")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g")

# Include directories
include_directories(${CMAKE_SOURCE_DIR}/include)

# Find required packages
find_package(GTest QUIET)

# Create directories for build outputs
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/tests)

# Set output directories
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# Core LBM library
set(LBM_CORE_SOURCES
    src/lbm/solver_base.cpp
    src/lbm/single_phase.cpp
    src/lbm/multiphase.cpp
    src/lbm/stability_monitor.cpp
    src/utils/math_utils.cpp
    src/utils/io_utils.cpp
    src/utils/validation.cpp
)

# Analysis library
set(LBM_ANALYSIS_SOURCES
    src/analysis/h_theorem.cpp
)

# Create core LBM library
add_library(lbm_core STATIC ${LBM_CORE_SOURCES})
target_include_directories(lbm_core PUBLIC ${CMAKE_SOURCE_DIR}/include)

# Create analysis library
add_library(lbm_analysis STATIC ${LBM_ANALYSIS_SOURCES})
target_include_directories(lbm_analysis PUBLIC ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(lbm_analysis lbm_core)

# Legacy single-phase application
add_executable(lbm_sim src/main.cpp)

# Modern single-phase application (when implemented)
if(EXISTS ${CMAKE_SOURCE_DIR}/src/apps/single_phase_app.cpp)
    add_executable(single_phase_app src/apps/single_phase_app.cpp)
    target_link_libraries(single_phase_app lbm_core lbm_analysis)
endif()

# Modern multiphase application (when implemented)
if(EXISTS ${CMAKE_SOURCE_DIR}/src/apps/multiphase_app.cpp)
    add_executable(multiphase_app src/apps/multiphase_app.cpp)
    target_link_libraries(multiphase_app lbm_core lbm_analysis)
endif()

# Testing
if(GTest_FOUND)
    enable_testing()
    
    # Unit tests
    set(UNIT_TEST_SOURCES
        tests/unit/test_solver_base.cpp
        tests/unit/test_single_phase.cpp
        tests/unit/test_multiphase.cpp
        tests/unit/test_stability_monitor.cpp
        tests/unit/test_math_utils.cpp
        tests/unit/test_io_utils.cpp
        tests/unit/test_data_validator.cpp
        tests/unit/test_h_theorem.cpp
        tests/unit/test_collision_operators.cpp
    )
    
    foreach(test_source ${UNIT_TEST_SOURCES})
        get_filename_component(test_name ${test_source} NAME_WE)
        add_executable(${test_name} ${test_source})
        target_link_libraries(${test_name} lbm_core lbm_analysis GTest::gtest GTest::gtest_main)
        target_include_directories(${test_name} PRIVATE ${CMAKE_SOURCE_DIR}/include)
        add_test(NAME ${test_name} COMMAND ${test_name})
        set_target_properties(${test_name} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/tests)
    endforeach()
    
    # Integration tests
    if(EXISTS ${CMAKE_SOURCE_DIR}/tests/integration)
        file(GLOB INTEGRATION_TEST_SOURCES tests/integration/*.cpp)
        foreach(test_source ${INTEGRATION_TEST_SOURCES})
            get_filename_component(test_name ${test_source} NAME_WE)
            add_executable(${test_name} ${test_source})
            target_link_libraries(${test_name} lbm_core lbm_analysis GTest::gtest GTest::gtest_main)
            target_include_directories(${test_name} PRIVATE ${CMAKE_SOURCE_DIR}/include)
            add_test(NAME ${test_name} COMMAND ${test_name})
            set_target_properties(${test_name} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/tests)
        endforeach()
    endif()
    
    # Benchmark tests
    if(EXISTS ${CMAKE_SOURCE_DIR}/tests/benchmarks)
        file(GLOB BENCHMARK_TEST_SOURCES tests/benchmarks/*.cpp)
        foreach(test_source ${BENCHMARK_TEST_SOURCES})
            get_filename_component(test_name ${test_source} NAME_WE)
            add_executable(${test_name} ${test_source})
            target_link_libraries(${test_name} lbm_core lbm_analysis GTest::gtest GTest::gtest_main)
            target_include_directories(${test_name} PRIVATE ${CMAKE_SOURCE_DIR}/include)
            # Benchmarks are not run automatically with ctest
            set_target_properties(${test_name} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/tests)
        endforeach()
    endif()
    
    message(STATUS "Google Test found. Unit tests will be built.")
else()
    message(WARNING "Google Test not found. Tests will not be built.")
    message(STATUS "To install Google Test on Ubuntu/Debian: sudo apt-get install libgtest-dev")
    message(STATUS "To install Google Test on macOS: brew install googletest")
endif()

# Documentation (optional)
find_package(Doxygen QUIET)
if(DOXYGEN_FOUND)
    set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/docs/Doxyfile.in)
    set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
    
    if(EXISTS ${DOXYGEN_IN})
        configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
        add_custom_target(docs
            COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
            COMMENT "Generating API documentation with Doxygen"
            VERBATIM)
    endif()
endif()

# Install targets
install(TARGETS lbm_core lbm_analysis
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin)

install(DIRECTORY include/ DESTINATION include)

# Print configuration summary
message(STATUS "")
message(STATUS "=== LBM Poiseuille Configuration Summary ===")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "C++ standard: ${CMAKE_CXX_STANDARD}")
message(STATUS "Compiler: ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "Install prefix: ${CMAKE_INSTALL_PREFIX}")
message(STATUS "Google Test: ${GTest_FOUND}")
message(STATUS "Doxygen: ${DOXYGEN_FOUND}")
message(STATUS "===========================================")
message(STATUS "")