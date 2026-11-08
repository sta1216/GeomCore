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
  # GCC 13 under -fsanitize=thread reports a false positive
  # -Wmaybe-uninitialized for a std::variant moved through vector growth
  # (SolidEditor::add_curve inlined into brep_cylinder_test). Ordinary GCC
  # builds keep the warning; only the TSan configuration suppresses it.
  if(GEOM_TSAN AND CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    target_compile_options(geom_warnings INTERFACE -Wno-maybe-uninitialized)
  endif()
endif()

if(GEOM_WERROR)
  if(MSVC)
    target_compile_options(geom_warnings INTERFACE /WX)
  else()
    target_compile_options(geom_warnings INTERFACE -Werror)
  endif()
endif()
