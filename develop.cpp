#include <iostream>
#include <cmath>

#include "Core/Matrix.hpp"

int main() {

    using namespace FastMath;
    using namespace FastMath::Impl;
    Matrix<double, 3, 3> X = Matrix3d::Rand();
    Matrix<double, 3, 3> Y = Matrix3d::Identity() * 0.5;

    std::cout << X << std::endl;
    std::cout << Y << std::endl;
    std::cout << constrain(X, 0.2, 0.6) << std::endl;

    return 0;

}

