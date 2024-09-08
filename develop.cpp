#include <iostream>
#include "FastMath.hpp"
int main() {

    fm::Matrix3d A = {1,2,3,4,5,6};
    std::cout << A.rank() << std::endl; // 矩阵的秩
    std::cout << A.norm() << std::endl;     // 矩阵的L2范数
    std::cout << A.norm_one() << std::endl;     // 矩阵的L1范数
    std::cout << A.trace() << std::endl;    // 矩阵的迹
    std::cout << A.det() << std::endl; // 矩阵的行列式
    return 0;
}

