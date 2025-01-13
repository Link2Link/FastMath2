#include <iostream>
#include "FastMath.hpp"
#include "FastMath/Algorithm/CurveFitting.hpp"

int main() {
    constexpr size_t N = 5;
    fm::Vector<double, N> xc{0, 0.3, 0.6, 0.9, 1.2};
    fm::Vector<double, N> yc;

    for (int k = 0; k < N; ++k) {
        yc(k) = std::sin(xc(k));
    }

    double x = 0.5;
    auto w = FastMath::Algorithm::RBFtrain(xc, yc);

    std::cout << w << std::endl;
    std::cout << FastMath::Algorithm::RBFeval(0.5, xc, w) << std::endl;

    return 0;
}

