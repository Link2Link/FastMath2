#include "doctest/doctest.h"
#include <iostream>
#include "Impl/linear/matrix.hpp"
#include "Impl/linear/equation.hpp"
#include "Impl/nonlinear/solver.hpp"
#include <cmath>
#include<limits>
using namespace FastMath;
using namespace std;
using namespace FastMath::Impl;
using namespace FastMath::Impl::NL;

double dhrtf(double x)
{
    double z;
    z=(((((x-5.0)*x+3.0)*x+1.0)*x-7.0)*x+7.0)*x-20.0;
    return(z);
}

TEST_CASE("dhrt test") {

    int i,n;
    double x[6];
    n=dhrt<6>(-2.0,5.0,0.2,x,dhrtf, 1E-10);
    for (i=0; i<=n-1; i++)
        CHECK(dhrtf(x[i]) == doctest::Approx(0).epsilon(1E-8));
}

double newtf(double x)
{
    return(x*x*(x-1.0)-1.0);
}

double newtdf(double x)
{
    return(3.0*x*x - 2.0*x);
}

TEST_CASE("newt test") {
    int k;
    double x,eps;
    double newtf(double x), newtdf(double x);
    eps=0.000001;  x=1.5;
    k=newt(&x,newtf,newtdf);

    CHECK( k > 0);
    CHECK(newtf(x) == doctest::Approx(0).epsilon(1E-8));
}


double atknf(double x)
{
    return(6.0-x*x) ;
}

TEST_CASE("aitken test") {
    int k;
    double x;
    x = 0.0;
    k = atkn(&x,atknf);
    CHECK( k > 0);
    CHECK(atknf(x) == doctest::Approx(0).epsilon(1E-8));
}

TEST_CASE("qrrt test"){
    int i;
    double xr[6],xi[6];
    double a[7]={-30.0,10.5,-10.5,1.5,4.5,-7.5,1.5};

    i=qrrt<6>(a,xr,xi);

    double x = xr[0];

    double result = a[0];
    for (int k = 0; k < 6; ++k)
    {
        result += a[k + 1] * x;
        x *= xr[0];
    }

    CHECK(result == doctest::Approx(0).epsilon(1E-8));
}