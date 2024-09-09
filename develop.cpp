#include <iostream>
#include "FastMath.hpp"
int main() {

    double R_array[3*3] = {1,2,3,4,56,6,7,8,9};
    fm::Matrix3d R;
    R.loadFrom(R_array);    //loadFrom 会对数组尺寸进行检查
    std::cout << R << std::endl;

    R.loadFromArrayPtr(R_array);
    std::cout << R << std::endl;   // loadFromArrayPtr 只需传入受元素的指针，不会对尺寸进行检查

    R = R * 2;

    double R_array_out[9];
    R.copyTo(R_array_out);     // 将R中内部数据拷贝到数组R_array_out中
    std::cout << fm::Impl::MatToString(&R_array_out[0], 3, 3) << std::endl; //MatToString是对普通数组的打印支持

    R.copyToArrayPtr(R_array_out);     //将R中内部数据拷贝到R_array_out指向的首地址处
    std::cout << fm::Impl::MatToString(&R_array_out[0], 3, 3) << std::endl; //MatToString是对普通数组的打印支持



    return 0;
}

