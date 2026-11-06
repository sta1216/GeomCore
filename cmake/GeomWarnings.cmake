add_library(geom_warnings INTERFACE)

if(MSVC)
  target_compile_options(geom_warnings INTERFACE
    /W4 /permissive- /fp:precise /Zc:__cplusplus /Zc:preprocessor)
  target_compile_definitions(geom_warnings INTERFACE _CRT_SECURE_NO_WARNINGS NOMINMAX)
else()
  target_compile_options(geom_warnings INTERFACE
    -Wall -Wextra -Wconversion -Wshadow -Wpedantic
    -fno-fast-math -ffp-contract=off)
endif()

if(GEOM_WERROR)
  if(MSVC)
    target_compile_options(geom_warnings INTERFACE /WX)
  else()
    target_compile_options(geom_warnings INTERFACE -Werror)
  endif()
endif()
