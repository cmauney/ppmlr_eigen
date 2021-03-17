set(build_debug "$<CONFIG:Debug>")
set(build_release "$<CONFIG:Release>")

target_compile_features(ppmlr_cxx::flags
  INTERFACE
    cxx_std_20
)

target_compile_options(ppmlr_cxx::flags
  INTERFACE
    $<${build_debug}:
      "-g;-O0"
    >
    $<${build_release}:
     "-O3"
    >
)

target_include_directories(ppmlr_cxx::flags
  INTERFACE
    ${CMAKE_SOURCE_DIR}/include
)

