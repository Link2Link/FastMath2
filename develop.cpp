#include <iostream>
#include <cmath>

#include "FastMath.hpp"

int main() {

    using namespace FastMath;
    using namespace FastMath::Impl;
    Matrix<double, 2, 2> A;
    double theta = 1;
    A(0,0) = cos(theta);
    A(0,1) = sin(theta);
    A(1,0) = -sin(theta);
    A(1,1) = cos(theta);

    std::cout << A << std::endl;
    std::cout << A.logm() << std::endl;

    return 0;

}

