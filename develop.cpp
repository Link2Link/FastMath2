#include <iostream>
#include <cmath>
#include "FastMath.hpp"
#include "Algorithm/quadprog.hpp"

int main() {

    using namespace FastMath;
    using namespace FastMath::Algorithm;

    SquareMatrix<double, 2> H = {1, -1, -1, 2};
    Vector2d f = {-2, -6};
    Matrix<double, 3, 2> A = {1,1,-1,2,2,1};
    Vector3d b = {2,2,3};

    Matrix<double, 0, 2> Ae;
    Vector<double, 0> be;

    Vector<double, 2> x;

    std::cout << quadprog<2, 3, 0>(H, f, A, b, Ae, be, x) << std::endl;
    std::cout << x << std::endl;

    return 0;

}

