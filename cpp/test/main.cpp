#include <iostream>
#include "matrix.hpp"

namespace vm = vector_math;

int main(int argc, char* argv[]) {

  double M1_data[9];
  double M2_data[9];
  double M3_data[9];
  vm::matrix<3,3> M1(M1_data);
  vm::matrix<3,3> M2(M2_data);
  vm::matrix<3,3> M3(M3_data);

  for (unsigned int i=0; i<9; i++) {
    M1_data[i] = i;
    M2_data[i] = 10 * i;
    M3_data[i] = 100 * i;
  }

  auto M = M1 - M2 + M3;

  for (unsigned int row=0; row<3; row++) {
    for (unsigned int col=0; col<3; col++) {
      std::cout << M(row, col) << ", ";
    }
    std::cout << std::endl;
  }
  
  return 0;
}
