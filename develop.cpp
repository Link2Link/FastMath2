#include <iostream>
#include "Impl/basic/rand.hpp"
#include "Impl/linear/matrix.hpp"
#include "Impl/linear/equation.hpp"
#include <cmath>


#define N  5

int main() {
    constexpr int m = 5;
    constexpr int n = 2;
    using namespace FastMath::Impl;
    using namespace std;


    int i, j;
    double a[N][N];
    double r[N], x[N], b[N] = {1.0, 0.0, 0.0, 0.0, 1.0};
    double b_[N];
    for (i = 0; i < N; i++)
        for (j = 0; j < N; j++) a[i][j] = 1.0 / (1.0 + i + j);
    cout << "A:\n";
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) cout << setw(10) << a[i][j];
        cout << endl;
    }
    cout << "b:\n";
    for (i = 0; i < N; i++) cout << b[i] << "    ";
    cout << endl;
    cout << "x:\n";
    bingt<N>(&a[0][0], b, x);
    for (i = 0; i < N; i++)
        printf("%f\n", x[i]);
    cout << "e:\n";
    for (i = 0; i < N; i++) {
        r[i] = 0;
        for (j = 0; j < N; j++) r[i] = r[i] + a[i][j] * x[j];
        r[i] = b[i] - r[i];
        cout << "r(" << i << ")=" << r[i] << endl;
    }

    tmul(&a[0][0], N,N, &x[0], N, 1, &b_[0]);

    cout << "b_:\n";
    for (i = 0; i < N; i++) cout << b[i] << "    ";
    cout << endl;
    return 0;


}

