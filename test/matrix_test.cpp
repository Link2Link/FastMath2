#include "doctest/doctest.h"
#include <iostream>
#include "Impl/linear/matrix.hpp"
#include <cmath>
#include<limits>
using namespace FastMath;
using namespace std;
using namespace FastMath::Impl;
TEST_CASE("matrix product") {
    int i,j;
    double a[4][5]={ {1.0,3.0,-2.0,0.0,4.0},
                     {-2.0,-1.0,5.0,-7.0,2.0},
                     {0.0,8.0,4.0,1.0,-5.0},
                     {3.0,-3.0,2.0,-4.0,1.0}};
    double c[4][3],b[5][3]={ {4.0,5.0,-1.0},
                             {2.0,-2.0,6.0},{7.0,8.0,1.0},
                             {0.0,3.0,-5.0},{9.0,8.0,-6.0}};

    double expected[4][3] = {    32,    15,    -9,
                                43 ,   27 ,   24,
                                -1 ,  -21 ,   77,
                                29 ,   33 ,   -5};
    tmul(&a[0][0],4,5,&b[0][0],5,3,&c[0][0]);
    CHECK(isEqual<4,3>(&c[0][0], &expected[0][0]));

}

TEST_CASE("matrix inverse")
{
    int i, j;
    double a[4][4]={ {0.2368,0.2471,0.2568,1.2671},
                     {1.1161,0.1254,0.1397,0.1490},
                     {0.1582,1.1675,0.1768,0.1871},
                     {0.1968,0.2071,1.2168,0.2271}};
    double b[4][4];
    double c[4][4];
    double expected[4][4];
    MatCopy<4,4>(MatRef(a), MatRef(b));
    i=inv<4>(&b[0][0]);
    tmul(MatRef(a),4,4,MatRef(b),4,4,&c[0][0]);

    setIdentity(MatRef(expected), 4, 4);

    CHECK(isEqual<4, 4>(MatRef(expected), MatRef(c)));

}

TEST_CASE("positive definite matrix inverse")
{
    int i,j;
    double a[4][4]={ {5.0,7.0,6.0,5.0},
                     {7.0,10.0,8.0,7.0},
                     {6.0,8.0,10.0,9.0},
                     {5.0,7.0,9.0,10.0}};
    double b[4][4],c[4][4];
    double expected[4][4];
    setIdentity(MatRef(expected), 4, 4);

    MatCopy<4,4>(MatRef(a), MatRef(b));
    i=ssgj<4>(&b[0][0]);

    tmul(&a[0][0],4,4,&b[0][0],4,4,&c[0][0]);

    CHECK(isEqual<4,4>(MatRef(c), MatRef(expected)));
}

