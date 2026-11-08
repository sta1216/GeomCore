add_library(geom_warnings INTERFACE)

if(MSVC)
  target_compile_options(geom_warnings INTERFACE
    /W4 /permissive- /fp:precise /Zc:__cplusplus /Zc:preprocessor /utf-8)
  target_compile_definitions(geom_warnings INTERFACE _CRT_SECURE_NO_WARNINGS NOMINMAX)
else()
  target_compile_options(geom_warnings INTERFACE
    -Wall -Wextra -Wconversion -Wshadow -Wpedantic
    -fno-fast-math -ffp-contract=off)
  # Clang folds -Wsign-conversion into -Wconversion and flags size_t indexing
  # in every loop; GCC and MSVC do not. The library headers are index-clean, but
  # keeping this Clang-only noise out lets -Werror stay meaningful elsewhere.
  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    target_compile_options(geom_warnings INTERFACE -Wno-sign-conversion)
  endif()
endif()

if(GEOM_WERROR)
  if(MSVC)
    target_compile_options(geom_warnings INTERFACE /WX)
  else()
    target_compile_options(geom_warnings INTERFACE -Werror)
  endif()
endif()
