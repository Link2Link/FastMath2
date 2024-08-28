/******************************************************************************
  文 件 名   : equation.hpp
  作    者   : Barry
  生成日期   : 2024-08-27
  最近修改   :
  功能描述   :  线性方程组求解
  函数列表   :
                gauss           全选主元高斯消去法求解n阶线性代数方程组Ax=b
                gauss_jordan    全选主元高斯-约当消去法求解n阶线性代数方程组Ax=b
                ldle            分解法求解系数矩阵为对称且右端有m组常数向量的线性代数方程AX=C
                chlk            对称正定方程组的平方根法 用Cholesky分解求解系数矩阵为对称正定，且右端有m组常数向量的线性代数方程AX=C
                seidel          高斯-赛德尔迭代法求解系数矩阵对角占优的方程Ax=b
                grad            共轭梯度法求解n阶对称正定方程组Ax=b,仅适用于对称正定矩阵A，内部不检查是否对称正定
                gmqr            线性最小二乘问题的Householder法
                bingt           病态方程组求解，Ax=b
  修改历史   :
  1.日    期   : 2024-08-27
    作    者   : Barry
    修改内容   : 创建文件
******************************************************************************/

#ifndef FASTMATH_EQUATION_HPP
#define FASTMATH_EQUATION_HPP

#include "Impl/linear/matrix.hpp"

namespace FastMath::Impl {

    /**
     * 全选主元高斯消去法求解n阶线性代数方程组Ax=b
     * @tparam n n         方程组阶数
     * @tparam T 模板类型，支持complex或double
     * @param a a[n][n]   系数矩阵。返回时被破坏。
     * @param b b[n]      常数向量，返回解向量。若系数矩阵奇异，则返回0向量。
     * @return
     */
    template<size_t n, class T>
    //模板声明T为类型参数支持复数计算
    //若系数矩阵奇异，则程序显示错误信息，函数返回0标志值。
    inline int gauss(T *a, T *b) {
        int k, i, j, is, p, q;
        int js[n];
        double d, t;
        T s;

        for (k = 0; k <= n - 2; k++)             //消元过程
        {
            d = 0.0;                       //全选主元
            for (i = k; i <= n - 1; i++)
                for (j = k; j <= n - 1; j++) {
                    t = ffabs(a[i * n + j]);
                    if (t > d) {
                        d = t;
                        js[k] = j;
                        is = i;
                    }

                }
            if (d + 1.0 == 1.0)               //系数矩阵奇异，求解失败！
            {
                for (i = 0; i < n; i++) b[i] = init(s);
                return 0;   // 系数矩阵奇异，求解失败！
            }
            if (js[k] != k)                 //列交换
            {
                for (i = 0; i <= n - 1; i++) {
                    p = i * n + k;
                    q = i * n + js[k];
                    s = a[p];
                    a[p] = a[q];
                    a[q] = s;
                }
            }
            if (is != k)                    //行交换
            {
                for (j = k; j <= n - 1; j++) {
                    p = k * n + j;
                    q = is * n + j;
                    s = a[p];
                    a[p] = a[q];
                    a[q] = s;
                }
                s = b[k];
                b[k] = b[is];
                b[is] = s;
            }
            s = a[k * n + k];
            for (j = k + 1; j <= n - 1; j++)        //归一化
            {
                p = k * n + j;
                a[p] = a[p] / s;
            }
            b[k] = b[k] / s;
            for (i = k + 1; i <= n - 1; i++)        //消元
            {
                for (j = k + 1; j <= n - 1; j++) {
                    p = i * n + j;
                    a[p] = a[p] - a[i * n + k] * a[k * n + j];
                }
                b[i] = b[i] - a[i * n + k] * b[k];
            }
        }
        s = a[(n - 1) * n + n - 1];
        if (ffabs(s) + 1.0 == 1.0)        //系数矩阵奇异，求解失败！
        {
            for (i = 0; i < n; i++) b[i] = init(s);
            return 0;   // 系数矩阵奇异，求解失败！
        }
        b[n - 1] = b[n - 1] / s;                //回代过程
        for (i = n - 2; i >= 0; i--) {
            s = init(s);
            for (j = i + 1; j <= n - 1; j++) s = s + a[i * n + j] * b[j];
            b[i] = b[i] - s;
        }
        js[n - 1] = n - 1;
        for (k = n - 1; k >= 0; k--)       //恢复
            if (js[k] != k) {
                s = b[k];
                b[k] = b[js[k]];
                b[js[k]] = s;
            }
        return 1;
    }


    /**
     * 全选主元高斯-约当消去法求解n阶线性代数方程组Ax=b
     * @tparam n n         方程组阶数
     * @tparam T 模板类型，支持complex或double
     * @param a a[n][n]   系数矩阵。返回时被破坏。
     * @param b b[n]      常数向量，返回解向量。若系数矩阵奇异，则返回0向量。
     * @return
     */
    template<size_t n, typename T>
    int gauss_jordan(T *a, T *b) {
        int k, i, j, is, p, q;
        int js[n];
        double d, t;
        T s;
        for (k = 0; k <= n - 1; k++)             //消去过程
        {
            d = 0.0;                       //全选主元
            for (i = k; i <= n - 1; i++)
                for (j = k; j <= n - 1; j++) {
                    t = ffabs(a[i * n + j]);
                    if (t > d) {
                        d = t;
                        js[k] = j;
                        is = i;
                    }

                }
            if (d + 1.0 == 1.0)               //系数矩阵奇异，求解失败！
            {
                for (i = 0; i < n; i++) b[i] = init(s);
                return 0; //系数矩阵奇异，求解失败！
            }
            if (js[k] != k)                 //列交换
            {
                for (i = 0; i <= n - 1; i++) {
                    p = i * n + k;
                    q = i * n + js[k];
                    s = a[p];
                    a[p] = a[q];
                    a[q] = s;
                }
            }
            if (is != k)                    //行交换
            {
                for (j = k; j <= n - 1; j++) {
                    p = k * n + j;
                    q = is * n + j;
                    s = a[p];
                    a[p] = a[q];
                    a[q] = s;
                }
                s = b[k];
                b[k] = b[is];
                b[is] = s;
            }
            s = a[k * n + k];
            for (j = k + 1; j <= n - 1; j++)        //归一化
            {
                p = k * n + j;
                a[p] = a[p] / s;
            }
            b[k] = b[k] / s;
            for (i = 0; i <= n - 1; i++)        //消元
            {
                if (i != k) {
                    for (j = k + 1; j <= n - 1; j++) {
                        p = i * n + j;
                        a[p] = a[p] - a[i * n + k] * a[k * n + j];
                    }
                    b[i] = b[i] - a[i * n + k] * b[k];
                }
            }
        }
        for (k = n - 1; k >= 0; k--)       //恢复
            if (js[k] != k) {
                s = b[k];
                b[k] = b[js[k]];
                b[js[k]] = s;
            }
        return 1;
    }

    /**
     * 分解法求解系数矩阵为对称且右端有m组常数向量的线性代数方程AX=C
     * @param a a[n][n]        存放系数矩阵，返回时将被破坏。 a为对称矩阵
     * @param n n              方程组的阶数。
     * @param m m              方程组右端常数向量的组数。
     * @param c c[n][m]        存放方程组右端m组常数向量。返回m组解向量。
     * @return 函数返回标志值。若<0则表示系数矩阵非对称；若=0则表示失败；若>0则表示正常。
     */


    int ldle(double a[], int n, int m, double c[]) {
        int i, j, l, k, u, v, w, k1, k2, k3;
        double p;
        for (i = 0; i < n; i++)
            for (j = 0; j < i - 1; j++)
                if (a[i * n + j] != a[j * n + i]) {
                    return -1;  //矩阵不对称
                }
        if (fabs(a[0]) + 1.0 == 1.0) return 0;  // 求解失败
        for (i = 1; i <= n - 1; i++) {
            u = i * n;
            a[u] = a[u] / a[0];
        }
        for (i = 1; i <= n - 2; i++) {
            u = i * n + i;
            for (j = 1; j <= i; j++) {
                v = i * n + j - 1;
                l = (j - 1) * n + j - 1;
                a[u] = a[u] - a[v] * a[v] * a[l];
            }
            p = a[u];
            if (fabs(p) + 1.0 == 1.0) return 0; // 求解失败
            for (k = i + 1; k <= n - 1; k++) {
                u = k * n + i;
                for (j = 1; j <= i; j++) {
                    v = k * n + j - 1;
                    l = i * n + j - 1;
                    w = (j - 1) * n + j - 1;
                    a[u] = a[u] - a[v] * a[l] * a[w];
                }
                a[u] = a[u] / p;
            }
        }
        u = n * n - 1;
        for (j = 1; j <= n - 1; j++) {
            v = (n - 1) * n + j - 1;
            w = (j - 1) * n + j - 1;
            a[u] = a[u] - a[v] * a[v] * a[w];
        }
        p = a[u];
        if (fabs(p) + 1.0 == 1.0) return 0; // 求解失败
        for (j = 0; j <= m - 1; j++)
            for (i = 1; i <= n - 1; i++) {
                u = i * m + j;
                for (k = 1; k <= i; k++) {
                    v = i * n + k - 1;
                    w = (k - 1) * m + j;
                    c[u] = c[u] - a[v] * c[w];
                }
            }
        for (i = 1; i <= n - 1; i++) {
            u = (i - 1) * n + i - 1;
            for (j = i; j <= n - 1; j++) {
                v = (i - 1) * n + j;
                w = j * n + i - 1;
                a[v] = a[u] * a[w];
            }
        }
        for (j = 0; j <= m - 1; j++) {
            u = (n - 1) * m + j;
            c[u] = c[u] / p;
            for (k = 1; k <= n - 1; k++) {
                k1 = n - k;
                k3 = k1 - 1;
                u = k3 * m + j;
                for (k2 = k1; k2 <= n - 1; k2++) {
                    v = k3 * n + k2;
                    w = k2 * m + j;
                    c[u] = c[u] - a[v] * c[w];
                }
                c[u] = c[u] / a[k3 * n + k3];
            }
        }
        return (1);
    }

    /**
     * 对称正定方程组的平方根法
     * 用Cholesky分解求解系数矩阵为对称正定，且右端有m组常数向量的线性代数方程AX=C
     * @param a a[n][n]        存放系数矩阵，返回时将被破坏。 a为对称正定
     * @param n n              方程组的阶数。
     * @param m m              方程组右端常数向量的组数。
     * @param d c[n][m]        存放方程组右端m组常数向量。返回m组解向量。
     * @return 函数返回标志值。若<0则表示系数矩阵非对称；若=0则表示失败；若>0则表示正常。
     */
    int chlk(double a[], int n, int m, double d[]) {
        int i, j, k, u, v;
        if ((a[0] + 1.0 == 1.0) || (a[0] < 0.0)) return (0);
        a[0] = sqrt(a[0]);
        for (j = 1; j <= n - 1; j++) a[j] = a[j] / a[0];
        for (i = 1; i <= n - 1; i++) {
            u = i * n + i;
            for (j = 1; j <= i; j++) {
                v = (j - 1) * n + i;
                a[u] = a[u] - a[v] * a[v];
            }
            if ((a[u] + 1.0 == 1.0) || (a[u] < 0.0)) return (0);
            a[u] = sqrt(a[u]);
            if (i != (n - 1)) {
                for (j = i + 1; j <= n - 1; j++) {
                    v = i * n + j;
                    for (k = 1; k <= i; k++)
                        a[v] = a[v] - a[(k - 1) * n + i] * a[(k - 1) * n + j];
                    a[v] = a[v] / a[u];
                }
            }
        }
        for (j = 0; j <= m - 1; j++) {
            d[j] = d[j] / a[0];
            for (i = 1; i <= n - 1; i++) {
                u = i * n + i;
                v = i * m + j;
                for (k = 1; k <= i; k++)
                    d[v] = d[v] - a[(k - 1) * n + i] * d[(k - 1) * m + j];
                d[v] = d[v] / a[u];
            }
        }
        for (j = 0; j <= m - 1; j++) {
            u = (n - 1) * m + j;
            d[u] = d[u] / a[n * n - 1];
            for (k = n - 1; k >= 1; k--) {
                u = (k - 1) * m + j;
                for (i = k; i <= n - 1; i++) {
                    v = (k - 1) * n + i;
                    d[u] = d[u] - a[v] * d[i * m + j];
                }
                v = (k - 1) * n + k - 1;
                d[u] = d[u] / a[v];
            }
        }
        return (1);
    }

    /**
     * 高斯-赛德尔迭代法求解系数矩阵对角占优的方程Ax=b
     * @param a a[n][n]     系数矩阵
     * @param b b[n]        常数向量
     * @param n n           方程组阶数
     * @param x x[n]        返回满足精度要求的解向量。若系数矩阵非对角优势，返回解向量0。
     * @param eps eps         精度要求
     * @return 若系数矩阵非对角优势，返回0标志值。否则返回非0标志值,表示迭代次数
     */

    int seidel(double *a, double *b, int n, double *x, double eps = 1E-8) {
        int i, j, u, v;
        double p, t, s, q;
        int it = 0;
        for (i = 0; i <= n - 1; i++) {
            u = i * n + i;
            p = 0.0;
            x[i] = 0.0;           //置解向量初值
            for (j = 0; j <= n - 1; j++)
                if (i != j) {
                    v = i * n + j;
                    p = p + fabs(a[v]);
                }
            if (p >= fabs(a[u]))     //检查系数矩阵是否对角优势
            {
                return 0;   //系数矩阵非对角优势
            }
        }
        p = eps + 1.0;
        while (p >= eps) {
            it++;
            p = 0.0;
            for (i = 0; i <= n - 1; i++) {
                t = x[i];
                s = 0.0;
                for (j = 0; j <= n - 1; j++)
                    if (j != i) s = s + a[i * n + j] * x[j];
                x[i] = (b[i] - s) / a[i * n + i];
                q = fabs(x[i] - t) / (1.0 + fabs(x[i]));
                if (q > p) p = q;
            }
        }
        return it;
    }

    /**
     * 使用共轭梯度法求解n阶对称正定方程组Ax=b,仅适用于对称正定矩阵A，内部不检查是否对称正定
     * @tparam n    方程组阶数。
     * @param a     a[n][n]      存放对称正定矩阵。
     * @param b     b[n]         存放方程组右端常数向量。
     * @param x     x[n]         返回方程组解向量。
     * @param eps   控制精度要求。
     */
    template<size_t n>
    void grad(double a[], double b[], double x[], double eps = 1E-8) {
        int i, k;
        double alpha, beta, d, e;
        double p[n], r[n], s[n], q[n];
        for (i = 0; i <= n - 1; i++) {
            x[i] = 0.0;
            p[i] = b[i];
            r[i] = b[i];
        }
        i = 0;
        while (i <= n - 1) {
            tmul(a, n, n, p, n, 1, s);
            d = 0.0;
            e = 0.0;
            for (k = 0; k <= n - 1; k++) {
                d = d + p[k] * b[k];
                e = e + p[k] * s[k];
            }
            alpha = d / e;
            for (k = 0; k <= n - 1; k++)
                x[k] = x[k] + alpha * p[k];
            tmul(a, n, n, x, n, 1, q);
            d = 0.0;
            for (k = 0; k <= n - 1; k++) {
                r[k] = b[k] - q[k];
                d = d + r[k] * s[k];
            }
            beta = d / e;
            d = 0.0;
            for (k = 0; k <= n - 1; k++) d = d + r[k] * r[k];
            d = sqrt(d);
            if (d < eps) {
                return;
            }
            for (k = 0; k <= n - 1; k++)
                p[k] = r[k] - beta * p[k];
            i = i + 1;
        }
        return;
    }


    /**
     * 线性最小二乘问题的Householder法
     * @tparam m    m              方程个数，也是系数矩阵行数。
     * @tparam n    n              未知数个数，也是系数矩阵的列数。要求 m>=n 。
     * @param a     a[m][n]        超定方程组的系数矩阵，返回时存放QR分解式中的R矩阵。
     * @param b     b[m]           存放方程组右端的常数向量。返回时前n个分量存放方程组的最小二乘解。
     * @param q     q[m][m]        返回QR分解式中的正交矩阵Q。
     * @return      函数返回标志值。若=0则表示失败；否则表示正常。
     */
    template<size_t m, size_t n>
    int gmqr(double a[], double b[], double q[]) {
        int i, j;
        double d;
        double c[n];
        i = maqr(a, m, n, q);
        if (i == 0) {
            return (0);
        }
        for (i = 0; i <= n - 1; i++) {
            d = 0.0;
            for (j = 0; j <= m - 1; j++) d = d + q[j * m + i] * b[j];
            c[i] = d;
        }
        b[n - 1] = c[n - 1] / a[n * n - 1];
        for (i = n - 2; i >= 0; i--) {
            d = 0.0;
            for (j = i + 1; j <= n - 1; j++) d = d + a[i * n + j] * b[j];
            b[i] = (c[i] - d) / a[i * n + i];
        }
        return (1);
    }

    /**
     * 病态方程组求解，Ax=b
     * @tparam n    n           方程组阶数
     * @param a     a[n][n]     系数矩阵
     * @param b     b[n]        常数向量
     * @param x     x[n]        返回解向量。若系数矩阵奇异，返回0向量。
     * @param eps   eps         精度要求
     * @param mat_it    最大迭代次数
     * @return      若系数矩阵奇异或校正达到最大迭代次数次还不满足精度要求，返回0标志值。正常则返回迭代次数
     */
    template<size_t n>
    int bingt(double a[], double b[], double x[], double eps=1E-8, int mat_it=10) {
        int i, j, k;
        double q, qq;
        double p[n * n], r[n], e[n];
        k = 0;
        for (i = 0; i <= n - 1; i++)
            for (j = 0; j <= n - 1; j++) p[i * n + j] = a[i * n + j];
        for (i = 0; i <= n - 1; i++) x[i] = b[i];
        i = gauss<n>(p, x);
        if (i == 0) {
            return 0;
        }
        q = 1.0 + eps;
        while (q >= eps) {
            if (k == mat_it) {
                return 0;   //校正达到最大次数
            }
            k = k + 1;
            for (i = 0; i <= n - 1; i++) {
                e[i] = 0;
                for (j = 0; j <= n - 1; j++) e[i] = e[i] + a[i * n + j] * x[j];
            }
            for (i = 0; i <= n - 1; i++) r[i] = b[i] - e[i];
            for (i = 0; i <= n - 1; i++)
                for (j = 0; j <= n - 1; j++) p[i * n + j] = a[i * n + j];
            i = gauss<n>(p, r);
            if (i == 0) {
                return 0;
            }
            q = 0.0;
            for (i = 0; i <= n - 1; i++) {
                qq = fabs(r[i]) / (1.0 + fabs(x[i] + r[i]));
                if (qq > q) q = qq;
            }
            for (i = 0; i <= n - 1; i++) x[i] = x[i] + r[i];
        }
        return k;
    }


}


#endif //FASTMATH_EQUATION_HPP
