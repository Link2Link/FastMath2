#include <iostream>
#include "FastMath.hpp"
#include "FastMath/Algorithm/CurveFitting.hpp"

#include "vector"
using namespace std;


using namespace std;



using namespace std;




int main() {
    constexpr size_t N = 11;
    fm::Vector<double, N> xc{0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
    fm::Vector<double, N> yc{1, -1, 2, -2, 3, -4, 5, -6, 7, -8, 0};

    FastMath::Algorithm::CubicSpline spline(xc, yc);
    // vector<double> test_points = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0};
    for (double t = -0.1; t < 1.1; t+=0.01) {

        std::cout << t << "," << spline.interpolate(t) << std::endl;
    }




    return 0;
}

