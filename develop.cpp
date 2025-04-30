#include <iostream>
#include "FastMath.hpp"
#include "FastMath/Algorithm/CurveFitting.hpp"

#include "vector"
using namespace std;



int main() {
    double A[3][4] = {1,2,3,4,5,6,7,8,9,10,11,12};
    double B[3] = {1,2,3};
    auto x = fm::MAT<3,4>(A);
    std::cout << fm::MAT<3,4>(A);
    std::cout << fm::MAT<3>(B).T();

    std::cout << fm::MAT<3,4>(A).resize<1,12>();

    return 0;
}

