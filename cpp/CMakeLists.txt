cmake_minimum_required(VERSION 3.0)
add_library(matrix-math INTERFACE)
install(FILES include/matrix.hpp include/so3.hpp
  DESTINATION include/matrix-math)

