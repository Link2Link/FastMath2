/******************************************************************************
  文 件 名   : optimization.hpp
  作    者   : Barry
  生成日期   : 2024-08-28
  最近修改   :
  功能描述   : 
  函数列表   :
  修改历史   :
  1.日    期   : 2024-08-28
    作    者   : Barry
    修改内容   : 创建文件
******************************************************************************/

#ifndef FASTMATH2_OPTIMIZATION_HPP
#define FASTMATH2_OPTIMIZATION_HPP

#include <cmath>
#include <iostream>
#include "Impl/Common.hpp"
#include "Impl/basic/rand.hpp"
#include "Impl/linear/matrix.hpp"
#include "Impl/linear/equation.hpp"

namespace FastMath::Impl::NL {


    /**
     * 一维极值连分式法求f(x)的极值, 可用于单变量无约束优化问题
     * @param x     存放极值点初值。返回极值点。
     * @param f     计算目标函数f(x)值的函数名。
     * @param df    计算目标函数一阶导数值的函数名。
     * @param eps   控制精度要求。
     * @return      函数返回标志值。若>0为极大值点；若<0为极小值点;若=0则不是极值点。
     */

    inline int max1(double *x, double (*f)(double), double (*df)(double), double eps = 1E-10, int max_it = 20) {
        int i, j, m, jt, flag, k;
        double xx, h, h1, h2, dx, y[10], b[10], z;
        flag = max_it;           //最大迭代次数
        k = 0;
        jt = 1;
        h2 = 0.0;
        while (jt == 1) {
            j = 0;
            while (j <= 7) {
                if (j <= 2)
                    xx = *x + 0.01 * j;
                else
                    xx = h2;
                z = (*df)(xx);
                if (fabs(z) < eps) {
                    jt = 0;
                    j = 10;
                } else {
                    h1 = z;
                    h2 = xx;
                    if (j == 0) {
                        y[0] = h1;
                        b[0] = h2;
                    } else {
                        y[j] = h1;
                        m = 0;
                        i = 0;
                        while ((m == 0) && (i <= j - 1)) {
                            if (fabs(h2 - b[i]) + 1.0 == 1.0)
                                m = 1;
                            else
                                h2 = (h1 - y[i]) / (h2 - b[i]);
                            i = i + 1;
                        }
                        b[j] = h2;
                        if (m != 0)
                            b[j] = 1.0e+35;
                        h2 = 0.0;
                        for (i = j - 1; i >= 0; i--)
                            h2 = -y[i] / (b[i + 1] + h2);
                        h2 = h2 + b[0];
                    }
                    j = j + 1;
                }
            }
            *x = h2;
            k = k + 1;
            if (k == flag)
                jt = 0; // 达到最大迭代次数
        }
        xx = *x;
        h = (*f)(xx);
        if (fabs(xx) <= 1.0)
            dx = 1.0e-05;
        else
            dx = fabs(xx * 1.0e-05);
        xx = *x - dx;
        h1 = (*f)(xx);
        xx = *x + dx;
        h2 = (*f)(xx);
        if ((h1 + h2 - 2.0 * h) > 0.0)
            k = -1;
        else if ((h1 + h2 - 2.0 * h) < 0.0)
            k = 1;
        else
            k = 0;  // 达到最大迭代次数
        return (k);
    }

    /**
     * n维极值连分式法, 可用于求解n维无约束优化问题
     * @param n         自变量个数。
     * @param x[n]      存放极值点初值。返回极值点。
     * @param f         计算目标函数值与各偏导数值的函数名。
     *                  函数原型为 double (*f)(double [] x, int n, int j)
     *                  n为自变量个数，j=0时返回目标函数值，j=1，2，...返回目标函数对第1,2...个变量的偏导数值
     * @param eps       控制精度要求。
     * @param max_it    最大迭代次数
     * @return  函数返回标志值。若>0为极大值点；若<0为极小值点;若=0为鞍点。
     */

    inline int maxn(int n, double x[], double (*f)(double [], int, int), double eps = 1E-10, int max_it = 20) {
        int i, j, m, kk, jt, il, k;
        double y[10], b[10], p, z, t, h1, h2, ff, dx;
        k = 0;
        jt = max_it;
        h2 = 0.0;
        while (jt != 0) {
            t = 0.0;
            for (i = 1; i <= n; i++) {
                ff = (*f)(x, n, i);
                t = t + fabs(ff);
            }
            if (t < eps) jt = 0;
            else {
                for (i = 0; i <= n - 1; i++) {
                    il = 5;
                    while (il != 0) {
                        j = 0;
                        t = x[i];
                        il = il - 1;
                        while (j <= 7) {
                            if (j <= 2) z = t + j * 0.01;
                            else z = h2;
                            x[i] = z;
                            ff = (*f)(x, n, i + 1);
                            if (fabs(ff) + 1.0 == 1.0) {
                                j = 10;
                                il = 0;
                            } else {
                                h1 = ff;
                                h2 = z;
                                if (j == 0) {
                                    y[0] = h1;
                                    b[0] = h2;
                                } else {
                                    y[j] = h1;
                                    m = 0;
                                    kk = 0;
                                    while ((m == 0) && (kk <= j - 1)) {
                                        p = h2 - b[kk];
                                        if (fabs(p) + 1.0 == 1.0) m = 1;
                                        else h2 = (h1 - y[kk]) / p;
                                        kk = kk + 1;
                                    }
                                    b[j] = h2;
                                    if (m != 0) b[j] = 1.0e+35;
                                    h2 = 0.0;
                                    for (kk = j - 1; kk >= 0; kk--) h2 = -y[kk] / (b[kk + 1] + h2);
                                    h2 = h2 + b[0];
                                }
                                j = j + 1;
                            }
                        }
                        x[i] = h2;
                    }
                    x[i] = z;
                }
                jt = jt - 1;
            }
        }
        k = 1;
        ff = (*f)(x, n, 0);
        x[n] = ff;
        dx = 0.00001;
        t = x[0];
        x[0] = t + dx;
        h1 = (*f)(x, n, 0);
        x[0] = t - dx;
        h2 = (*f)(x, n, 0);
        x[0] = t;
        t = h1 + h2 - 2.0 * ff;
        if (t > 0.0) k = -1;
        j = 1;
        jt = 1;
        while (jt == 1) {
            j = j + 1;
            dx = 0.00001;
            jt = 0;
            t = x[j - 1];
            x[j - 1] = t + dx;
            h2 = (*f)(x, n, 0);
            x[j - 1] = t - dx;
            h1 = (*f)(x, n, 0);
            x[j - 1] = t;
            t = h1 + h2 - 2.0 * ff;
            if ((t * k < 0.0) && (j < n)) jt = 1;
        }
        if (t * k > 0.0) k = 0;
        return (k);
    }


    /**
     * 不等式约束线性规划问题
     * @tparam m    不等式约束条件个数。
     * @tparam n    变量个数。
     * @param a[m][m+n]     左边n列存放不等式约束条件左端的系数矩阵，右边m列为单位矩阵。
     * @param b[m]          存放不等式约束条件右端项值。
     * @param c[m+n]        存放目标函数中的系数，其中后m个分量为0。
     * @param x[m+n]        前n个分量返回目标函数f的极小值点的n个坐标；第n+1个分量返回目标函数f的极小值。其余用于内部使用
     * @return
     */
    template<size_t m, size_t n>
    int lplq(double a[], double b[], double c[], double x[]) {
        int i, mn, k, j;
        double s, z, dd, y;
        int js[m];
        double p[m * m];
        double d[m * (m + n)];
        for (i = 0; i <= m - 1; i++) js[i] = n + i;
        mn = m + n;
        s = 0.0;
        while (1 == 1) {
            for (i = 0; i <= m - 1; i++)
                for (j = 0; j <= m - 1; j++) p[i * m + j] = a[i * mn + js[j]];
            k = inv<m>(p);
            if (k == 0) {
                x[n] = s;
                return (k);
            }
            tmul(p, m, m, a, m, mn, d);
            for (i = 0; i <= mn - 1; i++) x[i] = 0.0;
            for (i = 0; i <= m - 1; i++) {
                s = 0.0;
                for (j = 0; j <= m - 1; j++) s = s + p[i * m + j] * b[j];
                x[js[i]] = s;
            }
            k = -1;
            dd = 1.0e-35;
            for (j = 0; j <= mn - 1; j++) {
                z = 0.0;
                for (i = 0; i <= m - 1; i++) z = z + c[js[i]] * d[i * mn + j];
                z = z - c[j];
                if (z > dd) {
                    dd = z;
                    k = j;
                }
            }
            if (k == -1) {
                s = 0.0;
                for (j = 0; j <= n - 1; j++) s = s + c[j] * x[j];
                x[n] = s;
                return (1);
            }
            j = -1;
            dd = 1.0e+20;
            for (i = 0; i <= m - 1; i++)
                if (d[i * mn + k] >= 1.0e-20) {
                    y = x[js[i]] / d[i * mn + k];
                    if (y < dd) {
                        dd = y;
                        j = i;
                    }
                }
            if (j == -1) {
                x[n] = s;
                return (-1);
            }
            js[j] = k;
        }
        return 0;
    }


}


#endif //FASTMATH2_OPTIMIZATION_HPP
