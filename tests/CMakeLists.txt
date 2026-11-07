include(FetchContent)

if(EXISTS ${PROJECT_SOURCE_DIR}/third_party/doctest/doctest/doctest.h)
  add_library(doctest INTERFACE)
  add_library(doctest::doctest ALIAS doctest)
  target_include_directories(doctest INTERFACE ${PROJECT_SOURCE_DIR}/third_party/doctest)
else()
  # doctest v2.4.11 declares a CMake minimum older than 3.5, which CMake 4.x
  # rejects outright. This is CMake's documented escape hatch; it only relaxes
  # the policy check for the fetched subproject declared below.
  if(CMAKE_VERSION VERSION_GREATER_EQUAL 4.0)
    set(CMAKE_POLICY_VERSION_MINIMUM 3.5)
  endif()
  FetchContent_Declare(doctest
    GIT_REPOSITORY https://github.com/doctest/doctest.git
    GIT_TAG v2.4.11
    GIT_SHALLOW TRUE)
  # Test-only dependency: keep doctest's install() rules out of GeomCore's
  # install tree. This must be a CACHE entry with FORCE because doctest's own
  # cmake_minimum_required(VERSION 3.0) resets CMP0077, so its option() call
  # would ignore a plain variable and re-default DOCTEST_NO_INSTALL to OFF.
  set(DOCTEST_NO_INSTALL ON CACHE BOOL "Skip installing the fetched doctest" FORCE)
  FetchContent_MakeAvailable(doctest)
endif()

set(GEOM_TEST_SOURCES
  support/test_main.cpp
  unit/test_smoke.cpp
  unit/core_result_test.cpp
  unit/numeric_expansion_test.cpp
  unit/predicates2_orient_test.cpp
  unit/predicates2_incircle_test.cpp
  unit/core_handle_arena_test.cpp
  unit/primitives_test.cpp
  unit/curves_test.cpp
  unit/surfaces_test.cpp
  unit/hull2d_test.cpp
  unit/delaunay2_test.cpp
  unit/bvh_test.cpp
  unit/coverage_gaps_test.cpp
  unit/polygon2_test.cpp
  unit/kdtree_test.cpp
  unit/triangulate2_test.cpp
  unit/arrangement2_test.cpp
  unit/boolean2_test.cpp
)

list(APPEND GEOM_TEST_SOURCES integration/m0_exit_test.cpp)

add_executable(geom_tests ${GEOM_TEST_SOURCES})
target_link_libraries(geom_tests PRIVATE geom_core doctest::doctest)
# doctest forward-declares std::ostream by default; MSVC then fails to compile
# CHECK on std::string_view (its operator<< needs the real <ostream>). Ask
# doctest to include the standard headers instead of forward-declaring them.
target_compile_definitions(geom_tests PRIVATE DOCTEST_CONFIG_USE_STD_HEADERS)
add_test(NAME geom_tests COMMAND geom_tests)
