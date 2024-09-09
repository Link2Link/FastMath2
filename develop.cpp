#include <iostream>
#include "FastMath.hpp"
int main() {

    double R_array[3*3] = {1,2,3,4,56,6,7,8,9};
    fm::Matrix3d R;
    R.loadFrom(R_array);    //loadFrom �������ߴ���м��
    std::cout << R << std::endl;

    R.loadFromArrayPtr(R_array);
    std::cout << R << std::endl;   // loadFromArrayPtr ֻ�贫����Ԫ�ص�ָ�룬����Գߴ���м��

    R = R * 2;

    double R_array_out[9];
    R.copyTo(R_array_out);     // ��R���ڲ����ݿ���������R_array_out��
    std::cout << fm::Impl::MatToString(&R_array_out[0], 3, 3) << std::endl; //MatToString�Ƕ���ͨ����Ĵ�ӡ֧��

    R.copyToArrayPtr(R_array_out);     //��R���ڲ����ݿ�����R_array_outָ����׵�ַ��
    std::cout << fm::Impl::MatToString(&R_array_out[0], 3, 3) << std::endl; //MatToString�Ƕ���ͨ����Ĵ�ӡ֧��



    return 0;
}

