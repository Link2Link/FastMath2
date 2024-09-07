#include <iostream>
#include <cmath>
#include "FastMath.hpp"
#include "Algorithm/quadprog.hpp"
#include "Algorithm/qcqp.hpp"
int main() {

    using namespace FastMath;
    using namespace FastMath::Algorithm;

//    double A[4] = {1,0,0,1};
    SquareMatrix<double, 2> A;
    A.setIdentity();
    A += 3 * Matrix2d::Identity() ;
    A(0, 1) = 2;
    A(1, 0) = 2;

//    double b[2] = {-2, -3};
    Vector2d b = {-2, -3};

//    double d[2] = {1, 1};
    Vector2d d = {1,2};
    double r = 0.5;

//    double x[2] = {};
    Vector2d x;

    QCQP(x, A, b, d, r);

    std::cout << x  << std::endl;

    std::cout << "cost : " <<  (0.5*x.transpose()*A*x + x.transpose()*b).Scalar() << std::endl;

    return 0;

}

