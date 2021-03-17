project(ppmlr_cxx)
cmake_minimum_required(VERSION 3.15)

enable_language(CXX)
include(CMakeDetermineCXXCompiler)
include(CheckCXXCompilerFlag)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

add_library(ppmlr_cxx::flags INTERFACE IMPORTED)

include(PPMLR_dependencies)
include(PPMLR_build_flags)

set(exe_name "ppmlr.x")
add_executable(${exe_name})
target_sources(${exe_name}
  PRIVATE
    src/main.cpp
)

target_link_libraries(${exe_name}
  PRIVATE
    ppmlr_cxx::flags
)

