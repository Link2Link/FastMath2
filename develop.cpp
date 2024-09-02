#include <iostream>
#include <cmath>
#include "FastMath.hpp"
#include "Algorithm/quadprog.hpp"

int main() {

    using namespace FastMath;
    using namespace FastMath::Algorithm;

    SquareMatrix<double, 4> T;
    SquareMatrix<double, 3> R;
    R = T.slice<3,3>(0,0);

    return 0;

}

