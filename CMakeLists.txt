cmake_minimum_required(VERSION 3.25)
project(GeomCore VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

option(GEOM_BUILD_TESTS "Build unit tests" ON)
option(GEOM_BUILD_BENCH "Build benchmarks" OFF)
option(GEOM_WERROR "Treat warnings as errors" OFF)

include(cmake/GeomWarnings.cmake)

# Coverage instrumentation comes from GCC/Clang's --coverage runtime, which
# MSVC does not provide. Reject the option loudly on MSVC instead of building
# an uninstrumented binary that silently reports 0% coverage.
option(GEOM_COVERAGE "Enable coverage instrumentation (GCC/Clang)" OFF)
if(GEOM_COVERAGE)
  if(MSVC)
    message(FATAL_ERROR "GEOM_COVERAGE requires GCC or Clang")
  endif()
  target_compile_options(geom_warnings INTERFACE --coverage -O0 -g)
  target_link_options(geom_warnings INTERFACE --coverage)
endif()

add_library(geom_core INTERFACE)
add_library(GeomCore::core ALIAS geom_core)
target_include_directories(geom_core INTERFACE
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>)
target_link_libraries(geom_core INTERFACE geom_warnings)

if(GEOM_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()

# Benchmarks are opt-in (GEOM_BUILD_BENCH=OFF by default) and are not
# registered with ctest: they measure and report, they are not pass/fail gates
# on shared CI. Run the executable manually to reproduce the numbers.
if(GEOM_BUILD_BENCH)
  add_executable(geom_bench_predicates bench/bench_predicates.cpp)
  target_link_libraries(geom_bench_predicates PRIVATE geom_core)
endif()
