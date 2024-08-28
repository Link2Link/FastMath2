#include <iostream>
#include "Impl/basic/rand.hpp"
#include "Impl/linear/matrix.hpp"
#include "Impl/linear/equation.hpp"
#include "Impl/nonlinear/solver.hpp"
#include <cmath>





int main() {

    using namespace FastMath::Impl;
    using namespace FastMath::Impl::NL;
    using namespace std;


    int i, n;
    double xr[6],xi[6],eps;
    double a[7]={-30.0,10.5,-10.5,1.5,4.5,-7.5,1.5};

    i=qrrt<6>(a,xr,xi);
    if (i>0)
    { for (i=0; i<=5; i++)
            printf("%lf, J %f\n", xr[i], xi[i]);
    }
    return 0;

}

