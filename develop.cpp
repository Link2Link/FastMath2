#include <iostream>
#include "FastMath.hpp"
int main() {

    fm::Matrix4d T;
    T.setIdentity();

    // 将矩阵左上角写为单位正乘3
    T.slice<3,3>(0,0) = fm::Matrix3d::Identity() * 3;
    std::cout << T << std::endl;

    // 将左上角3x3矩阵赋值给A
    fm::Matrix3d A = T.slice<3,3>(0,0);
    std::cout << A << std::endl;

    // slice的结果使用mat()后会生成一个新矩阵，不会改变原本的矩阵
    T.slice<3,3>(0,0).mat() += 3;
    std::cout << T << std::endl;

    // slice后直接进行+=，会改变原本的矩阵中的数据
    T.slice<3,3>(0,0) -= 3;
    std::cout << T << std::endl;

    // 若希望使用切片之后的结果做矩阵运算，例如求逆，需要使用mat()将其转为矩阵
    fm::Matrix3d U, S, V;
    auto svd_success = T.slice<3, 3>(0,0).mat().SVD(U, S, V);
    if (svd_success)
        std::cout << U*S*V.transpose() << std::endl;
    else
        std::cout << "svd failure" << std::endl;

    return 0;
}

