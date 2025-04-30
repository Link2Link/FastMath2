#include <iostream>
#include "FastMath.hpp"
#include "FastMath/Algorithm/CurveFitting.hpp"

#include "vector"
using namespace std;



int main() {
    double A[4][3] = {1,2,3,4,5,6,7,8,9,10,11,12};
    double B[3] = {1,2,3};

    std::cout << fm::Matrix<double, 4, 3>(A) * fm::Matrix<double, 3, 1>(B) << std::endl;

    return 0;
}

