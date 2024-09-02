#include <iostream>
#include <cmath>
#include "FastMath.hpp"
#include "Algorithm/quadprog.hpp"

int main() {

    using namespace FastMath;
    using namespace FastMath::Algorithm;

    SquareMatrix<double, 4> T = {1,2,3,4,5,6,7,8};
    SquareMatrix<double, 3> R = T.slice<3,3>(0,0);
    std::cout << R << std::endl;

    return 0;

}

