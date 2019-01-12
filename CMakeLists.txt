cmake_minimum_required (VERSION 3.1)
project (ppmlr_eigen VERSION 0.1 LANGUAGES CXX)

#set(PROJECT_LIBS "")

#set(CMAKE_CXX_STANDARD 17)
#set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package (Eigen3 3.3.39 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIRS} include/ /usr/include)

#add_compile_options(-g3 -Wall -Wextra -Og -ggdb -fmax-errors=1 -fno-omit-frame-pointer -fno-inline)
#add_compile_options(-O3 -march=native -mtune=native -Winline -DNDEBUG)
set(CMAKE_VERBOSE_MAKEFILE ON)
set(ISL_OPTIONS -ftree-vectorize -floop-interchange -floop-strip-mine -floop-block -ftree-loop-distribution -ftree-loop-distribute-patterns)
set(DBG_OPTIONS -g3 -Wall -Wextra -Og -ggdb -fmax-errors=1 -fno-omit-frame-pointer -fno-inline)
set(OPT_OPTIONS -O3 -march=native -mtune=native -Winline -DNDEBUG)

add_executable (ppm1d src/main.cpp)
target_compile_features(ppm1d PRIVATE cxx_std_17)
target_compile_options(ppm1d PRIVATE ${OPT_OPTIONS} ${ISL_OPTIONS})
target_link_libraries(ppm1d m)
