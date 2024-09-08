#include <iostream>
#include "FastMath.hpp"
int main() {

    fm::Matrix3d so3;
    so3(1, 2) = -fm::M_PI_F / 2;
    so3(2, 1) = fm::M_PI_F /2 ;

    std::cout << so3 << std::endl;
    std::cout << so3.expm() << std::endl;
    std::cout << so3.expm().logm() << std::endl;


    return 0;
}

