#include <iostream>
#include "FastMath.hpp"
int main() {

    fm::Matrix3d R(2);
    std::cout << R.power(3) << std::endl;
    std::cout << (R*R*R == R.power(3)) << std::endl;

    return 0;
}

