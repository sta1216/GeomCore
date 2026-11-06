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

add_library(geom_core INTERFACE)
add_library(GeomCore::core ALIAS geom_core)
target_include_directories(geom_core INTERFACE
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>)
target_link_libraries(geom_core INTERFACE geom_warnings)

if(GEOM_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()
