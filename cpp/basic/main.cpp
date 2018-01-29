#include <iostream>
#include "matrix.hpp"

int main(int argc, char* argv[]) {

  float M1_data[9];
  float M2_data[9];
  float M3_data[9];
  matrix M1(M1_data);
  matrix M2(M2_data);
  matrix M3(M3_data);

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
