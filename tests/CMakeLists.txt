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
  FetchContent_MakeAvailable(doctest)
endif()

set(GEOM_TEST_SOURCES
  support/test_main.cpp
  unit/test_smoke.cpp
)

add_executable(geom_tests ${GEOM_TEST_SOURCES})
target_link_libraries(geom_tests PRIVATE geom_core doctest::doctest)
add_test(NAME geom_tests COMMAND geom_tests)
