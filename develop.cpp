#include <iostream>
#include "FastMath.hpp"
int main() {

    fm::Matrix3d A = {1,2,3,4,5,6};
    std::cout << A.rank() << std::endl; // �������
    std::cout << A.norm() << std::endl;     // �����L2����
    std::cout << A.norm_one() << std::endl;     // �����L1����
    std::cout << A.trace() << std::endl;    // ����ļ�
    std::cout << A.det() << std::endl; // ���������ʽ
    return 0;
}

