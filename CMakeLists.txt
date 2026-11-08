cmake_minimum_required(VERSION 3.25)
project(GeomCore VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

option(GEOM_BUILD_TESTS "Build unit tests" ON)
option(GEOM_BUILD_BENCH "Build benchmarks" OFF)
option(GEOM_WERROR "Treat warnings as errors" OFF)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

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

option(GEOM_TSAN "Build with ThreadSanitizer (GCC/Clang)" OFF)
if(GEOM_TSAN)
  if(MSVC)
    message(FATAL_ERROR "GEOM_TSAN requires GCC or Clang")
  endif()
  target_compile_options(geom_warnings INTERFACE -fsanitize=thread -g)
  target_link_options(geom_warnings INTERFACE -fsanitize=thread)
endif()

add_library(geom_core INTERFACE)
add_library(GeomCore::core ALIAS geom_core)
target_include_directories(geom_core INTERFACE
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
# Export the language level so installed consumers compile the headers with
# C++20 instead of the compiler default.
target_compile_features(geom_core INTERFACE cxx_std_20)
target_link_libraries(geom_core INTERFACE geom_warnings)

# Installed package: headers plus the geom_core/geom_warnings interface targets.
# doctest is a test-only dependency and must never reach the install tree.
set_target_properties(geom_core PROPERTIES EXPORT_NAME core)
set_target_properties(geom_warnings PROPERTIES EXPORT_NAME warnings)

install(TARGETS geom_core geom_warnings EXPORT geomTargets)
install(DIRECTORY include/ DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}")
install(EXPORT geomTargets
  FILE geomTargets.cmake
  NAMESPACE GeomCore::
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/GeomCore")
configure_package_config_file(
  cmake/GeomCoreConfig.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/GeomCoreConfig.cmake"
  INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/GeomCore")
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/GeomCoreConfigVersion.cmake"
  VERSION "${PROJECT_VERSION}"
  COMPATIBILITY SameMajorVersion)
install(FILES
  "${CMAKE_CURRENT_BINARY_DIR}/GeomCoreConfig.cmake"
  "${CMAKE_CURRENT_BINARY_DIR}/GeomCoreConfigVersion.cmake"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/GeomCore")

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

  add_executable(geom_bench_m1b bench/bench_m1b.cpp)
  target_link_libraries(geom_bench_m1b PRIVATE geom_core)

  add_executable(geom_bench_m2 bench/bench_m2.cpp)
  target_link_libraries(geom_bench_m2 PRIVATE geom_core)
endif()
