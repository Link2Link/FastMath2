#include <iostream>
#include "FastMath.hpp"
int main() {

    fm::Matrix4d T;
    T.setIdentity();

    // ���������Ͻ�дΪ��λ����3
    T.slice<3,3>(0,0) = fm::Matrix3d::Identity() * 3;
    std::cout << T << std::endl;

    // �����Ͻ�3x3����ֵ��A
    fm::Matrix3d A = T.slice<3,3>(0,0);
    std::cout << A << std::endl;

    // slice�Ľ��ʹ��mat()�������һ���¾��󣬲���ı�ԭ���ľ���
    T.slice<3,3>(0,0).mat() += 3;
    std::cout << T << std::endl;

    // slice��ֱ�ӽ���+=����ı�ԭ���ľ����е�����
    T.slice<3,3>(0,0) -= 3;
    std::cout << T << std::endl;

    // ��ϣ��ʹ����Ƭ֮��Ľ�����������㣬�������棬��Ҫʹ��mat()����תΪ����
    fm::Matrix3d U, S, V;
    auto svd_success = T.slice<3, 3>(0,0).mat().SVD(U, S, V);
    if (svd_success)
        std::cout << U*S*V.transpose() << std::endl;
    else
        std::cout << "svd failure" << std::endl;

    return 0;
}

