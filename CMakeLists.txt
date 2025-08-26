cmake_minimum_required(VERSION 3.16)
project(rinox)

# Options
option(RINOX_ENABLE_ASAN "Enable AddressSanitizer" OFF)
option(RINOX_ENABLE_TESTS "Build unit tests" ON)
option(RINOX_ENABLE_EXPERIMENTS "Build experiments" OFF)
option(RINOX_VERBOSE_COMPILE "Print compilation-time parameters" ON)
set(RINOX_NUM_VARS_SIGN 6 CACHE STRING "Default cube size for signatures")
set(RINOX_MAX_CUTS_SIZE 6 CACHE STRING "Default maximum number of leaves in cuts")
set(RINOX_MAX_CUBE_SPFD 12 CACHE STRING "Default maximum cube size exactly represented")
set(RINOX_MAX_NUM_LEAVES 12 CACHE STRING "Default maximum number of leaves in the window")
set(RINOX_BENCHMARK_DIR "${CMAKE_SOURCE_DIR}/benchmarks" CACHE PATH "Directory containing benchmark suites")
set(RINOX_EXPERIMENTS_DIR "${CMAKE_SOURCE_DIR}/experiments" CACHE PATH "Directory containing the experiments")
set(RINOX_SOURCE_DIR "${CMAKE_SOURCE_DIR}" CACHE PATH "Source directory")


# Sanitizer
if (RINOX_ENABLE_ASAN)
  set(SANITIZE_FLAGS "-fsanitize=address -fno-omit-frame-pointer")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SANITIZE_FLAGS}")
  set(CMAKE_LINKER_FLAGS "${CMAKE_LINKER_FLAGS} ${SANITIZE_FLAGS}")
endif()

set(CMAKE_CXX_STANDARD 17)

# SIMD checks
if(UNIX)
  include(CheckCXXCompilerFlag)
  check_cxx_compiler_flag("-mavx2" HAS_MAVX2)
  check_cxx_compiler_flag("-mfma" HAS_MFMA)
  if(HAS_MAVX2 AND HAS_MFMA)
    message(STATUS "[rinox] Enabling AVX2 for SIMD operations")
    add_compile_options(-mavx2 -mfma)
  endif()
endif()

if(MSVC)
  add_compile_options(/arch:AVX2)
  add_compile_options(/EHsc /bigobj)
endif()

# Dependencies
add_subdirectory(lib/mockturtle)
add_subdirectory(lib/replxx)
add_subdirectory(lib/cli11)
add_subdirectory(inspect)

# RapidJSON (header-only) — submodule at lib/rapidjson
set(RAPIDJSON_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/lib/rapidjson/include)
if(NOT EXISTS "${RAPIDJSON_INCLUDE_DIR}/rapidjson/document.h")
  message(FATAL_ERROR "RapidJSON headers not found at ${RAPIDJSON_INCLUDE_DIR}.
Make sure you added the submodule: git submodule add https://github.com/Tencent/rapidjson.git lib/rapidjson")
endif()

# Drag (header-only) — submodule at lib/drag
set(DRAG_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/lib/drag/include)
if(NOT EXISTS "${DRAG_INCLUDE_DIR}/drag/graph.hpp")
  message(FATAL_ERROR "Drag headers not found at ${DRAG_INCLUDE_DIR}.
Make sure you added the submodule: git submodule add https://github.com/bigno78/drag.git lib/drag")
endif()

# Shared CLI logic as a library (exclude main.cpp)
file(GLOB_RECURSE MAD_CLI_SOURCES CONFIGURE_DEPENDS src/cli/*.cpp)
list(FILTER MAD_CLI_SOURCES EXCLUDE REGEX "src/cli/main.cpp$")
add_library(rho_cli ${MAD_CLI_SOURCES})

target_link_libraries(rho_cli PUBLIC mockturtle replxx CLI11::CLI11 )

# Include paths
target_include_directories(rho_cli
  PUBLIC
    ${CMAKE_SOURCE_DIR}/src          # so "cli/context.hpp" works
  SYSTEM PRIVATE
    ${RAPIDJSON_INCLUDE_DIR}         # RapidJSON for sources in rho_cli
  PUBLIC
    ${DRAG_INCLUDE_DIR}
)

target_compile_definitions(rho_cli PUBLIC
  RINOX_NUM_VARS_SIGN=${RINOX_NUM_VARS_SIGN}
  RINOX_MAX_CUTS_SIZE=${RINOX_MAX_CUTS_SIZE}
  RINOX_MAX_CUBE_SPFD=${RINOX_MAX_CUBE_SPFD}
  RINOX_MAX_NUM_LEAVES=${RINOX_MAX_NUM_LEAVES}
)

if(RINOX_VERBOSE_COMPILE)
  message(STATUS "Configuring mad-hatter with:")
  message(STATUS "  RINOX_NUM_VARS_SIGN = ${RINOX_NUM_VARS_SIGN}")
  message(STATUS "  RINOX_MAX_CUTS_SIZE = ${RINOX_MAX_CUTS_SIZE}")
  message(STATUS "  RINOX_MAX_CUBE_SPFD = ${RINOX_MAX_CUBE_SPFD}")
  message(STATUS "  RINOX_MAX_NUM_LEAVES = ${RINOX_MAX_NUM_LEAVES}")
endif()

# CLI executable
add_executable(rinox src/cli/main.cpp)
target_link_libraries(rinox PRIVATE rho_cli rinox_inspect )

# Public header-only interface library
add_library(rinox_headers INTERFACE)
add_library(rinox::rinox ALIAS rinox_headers)

target_include_directories(rinox_headers
  INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
  SYSTEM INTERFACE
    ${RAPIDJSON_INCLUDE_DIR}
)

target_compile_definitions(rinox_headers INTERFACE
  RINOX_NUM_VARS_SIGN=${RINOX_NUM_VARS_SIGN}
  RINOX_MAX_CUTS_SIZE=${RINOX_MAX_CUTS_SIZE}
  RINOX_MAX_CUBE_SPFD=${RINOX_MAX_CUBE_SPFD}
  RINOX_MAX_NUM_LEAVES=${RINOX_MAX_NUM_LEAVES}
)

# Expose to all targets that need it
add_library(rinox_build_options INTERFACE)
target_compile_definitions(rinox_build_options INTERFACE
  RINOX_BENCHMARK_DIR="$<IF:$<BOOL:${RINOX_BENCHMARK_DIR}>,${RINOX_BENCHMARK_DIR},${CMAKE_SOURCE_DIR}/benchmarks>"
  RINOX_EXPERIMENTS_DIR="$<IF:$<BOOL:${RINOX_EXPERIMENTS_DIR}>,${RINOX_EXPERIMENTS_DIR},${CMAKE_SOURCE_DIR}/experiments>"
  RINOX_SOURCE_DIR="$<IF:$<BOOL:${RINOX_SOURCE_DIR}>,${RINOX_SOURCE_DIR},${CMAKE_SOURCE_DIR}>"
)

# Make rho_cli consumers inherit headers/defines too
target_link_libraries(rho_cli PUBLIC rinox_build_options rinox_headers rinox_inspect )
# Tests
if (RINOX_ENABLE_TESTS)
  enable_testing()
  add_subdirectory(test)
endif()

if (RINOX_ENABLE_EXPERIMENTS)
  add_subdirectory(experiments)
endif()
