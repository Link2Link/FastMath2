#include <iostream>
#include <cmath>
#include "FastMath.hpp"
#include "Algorithm/quadprog.hpp"

int main() {

    using namespace FastMath;
    using namespace FastMath::Algorithm;


    Matrix<double, 3, 3> A = Matrix<double, 3, 3>::Rand();

    Vector3d x_{1,2,3};
    auto x = Vector2Dual(x_);
    auto y = Vector2Dual(Vector3d::Zeros());
    Dual<double, 3> z;
    A.setIdentity();

    y = A * x;

    for (size_t k = 0; k < 3; ++k)
        z += y(k);

    std::cout << z << std::endl;

    return 0;

}

