#include <iostream>
#include <Eigen/Dense>
#include "sumc.h"

int main(int argc, char** argv) {

Eigen::MatrixXd T(3, 4);

T << 1, 1, 1, 1,
     1, 1, 1, 1,
     1, 1, 1, 1;

std::cout << T << std::endl;
std::cout << sumL(T) << std::endl;
std::cout << sumE(T) << std::endl;
return 0;
}
