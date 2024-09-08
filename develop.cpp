#include <iostream>
#include "FastMath.hpp"
int main() {

    fm::Matrix<double, 4, 3> A = {1,2,3,4,5,6,1};

    fm::Matrix<double, 4, 4> Q;
    fm::Matrix<double, 4, 3> R;

    auto success = A.QR(Q, R);
    std::cout << success << std::endl;
    std::cout << A << std::endl;
    std::cout << Q << std::endl;
    std::cout << R << std::endl;
    std::cout << Q*R << std::endl;
    return 0;
}

