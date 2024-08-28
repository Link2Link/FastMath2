/******************************************************************************
  文 件 名   : solver.hpp
  作    者   : Barry
  生成日期   : 2024-08-28
  最近修改   :
  功能描述   :  非线性方程、方程组求解
  函数列表   :
  修改历史   :
  1.日    期   : 2024-08-28
    作    者   : Barry
    修改内容   : 创建文件
******************************************************************************/

#ifndef FASTMATH2_SOLVER_HPP
#define FASTMATH2_SOLVER_HPP

#include <cmath>
#include <iostream>
#include "Impl/linear/matrix.hpp"
#include "Impl/linear/equation.hpp"

namespace FastMath::Impl::NL {
    /**
     * 使用对分法搜索f(x)=0在区间[a, b]内的实根
     * @tparam m    实根个数的预估值。
     * @param a     求根区间的左端点。
     * @param b     求根区间的右端点。
     * @param h     搜索求根所采用的步长。
     * @param x[m]  存放返回的实根。实根个数由函数值返回。
     * @param f     方程左端函数f(x)的函数名。
     * @param eps   控制精度要求。
     * @return      函数返回搜索到的实根个数。若此值等于m，则有可能没有搜索完。
     */

    template<size_t m>
    inline int dhrt(double a, double b, double h, double x[],
                    double (*f)(double), double eps = 1E-8) {
        int n, js;
        double z, y, z1, y1, z0, y0;
        if (a > b) {
            z = a;
            a = b;
            b = z;
        }
        n = 0;
        z = a;
        y = (*f)(z);
        while ((z <= b + h / 2.0) && (n != m)) {
            if (fabs(y) < eps) {
                n = n + 1;
                x[n - 1] = z;
                z = z + h / 2.0;
                y = (*f)(z);
            } else {
                z1 = z + h;
                y1 = (*f)(z1);
                if (fabs(y1) < eps) {
                    n = n + 1;
                    x[n - 1] = z1;
                    z = z1 + h / 2.0;
                    y = (*f)(z);
                } else if (y * y1 > 0.0) {
                    y = y1;
                    z = z1;
                } else {
                    js = 0;
                    while (js == 0) {
                        if (fabs(z1 - z) < eps) {
                            n = n + 1;
                            x[n - 1] = (z1 + z) / 2.0;
                            z = z1 + h / 2.0;
                            y = (*f)(z);
                            js = 1;
                        } else {
                            z0 = (z1 + z) / 2.0;
                            y0 = (*f)(z0);
                            if (fabs(y0) < eps) {
                                x[n] = z0;
                                n = n + 1;
                                js = 1;
                                z = z0 + h / 2.0;
                                y = (*f)(z);
                            } else if ((y * y0) < 0.0) {
                                z1 = z0;
                                y1 = y0;
                            } else {
                                z = z0;
                                y = y0;
                            }
                        }
                    }
                }
            }
        }
        return (n);
    }


    /**
     * 牛顿迭代法求解非线性方程 f(x) = 0
     * @param x     存放方程根的初值。返回迭代终值。
     * @param f     方程左端函数f(x)的函数名。
     * @param df    方程左端函数f(x)一阶导函数名。
     * @param eps   控制精度要求。
     * @return
     */

    inline int newt(double *x, double (*f)(double), double (*df)(double), double eps = 1E-8) {
        int k, interation;
        double y, dy, d, p, x0, x1;
        interation = 500;         //最大迭代次数
        k = 0;
        x0 = *x;
        y = (*f)(x0);
        dy = (*df)(x0);
        d = eps + 1.0;
        while ((d >= eps) && (k != interation)) {
            if (fabs(dy) + 1.0 == 1.0)     //出现df(x)/dx=0
            {
                return (-1);
            }
            x1 = x0 - y / dy;            //迭代
            y = (*f)(x1);
            dy = (*df)(x1);
            d = fabs(x1 - x0);
            p = fabs(y);
            if (p > d) d = p;
            x0 = x1;
            k = k + 1;
        }
        *x = x0;
        return (k);
    }


    /**
     * aitken迭代法,求解 方程 f(x) = 0
     * @param x         存放方程根的初值。返回迭代终值。
     * @param f         函数f(x)的函数名。
     * @param eps       控制精度要求。
     * @param max_it    最大迭代次数
     * @return  函数返回迭代次数。
     */
    inline int atkn(double *x, double (*f)(double), double eps = 1E-8, int max_it = 500) {
        int flag, k, interation;
        double u, v, x0;
        interation = max_it;         //最大迭代次数
        k = 0;
        x0 = *x;
        flag = 0;
        while ((flag == 0) && (k != interation)) {
            k = k + 1;
            u = (*f)(x0) + x0;
            v = (*f)(u) + u;
            if (fabs(u - v) < eps) {
                x0 = v;
                flag = 1;
            } else
                x0 = v - (v - u) * (v - u) / (v - 2.0 * u + x0);
        }
        *x = x0;
        return (k);
    }


    /**
     * 多项式方程求根QR方法。
     * 将多项式构造为一上H矩阵的特征方程，对此矩阵求特征值进而得到多项式的所有根
     * @tparam n    多项式次数。
     * @param a     存放n次多项式的n+1个系数。 P(x) = a0 + a1 x + a2 x^2 + ... + an x^n
     * @param xr    返回n个根的实部。
     * @param xi    返回n个根的虚部。
     * @param eps   控制精度要求。
     * @return  函数返回在求上H矩阵特征值时返回的标志值。若>0则正常。
     */
    template<size_t n>
    inline int qrrt(double a[], double xr[], double xi[], double eps = 1E-12) {
        int i, j;
        double q[n * n];
        for (j = 0; j <= n - 1; j++)
            q[j] = -a[n - j - 1] / a[n];
        for (j = n; j <= n * n - 1; j++)
            q[j] = 0.0;
        for (i = 0; i <= n - 2; i++)
            q[(i + 1) * n + i] = 1.0;
        i = hhqr(q, n, xr, xi, eps);
        return (i);
    }


}


#endif //FASTMATH2_SOLVER_HPP
