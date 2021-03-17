find_package(Eigen3 3.3.9 REQUIRED)
find_package(fmt REQUIRED)

target_link_libraries(ppmlr_cxx::flags
  INTERFACE
    Eigen3::Eigen
    fmt::fmt
)
