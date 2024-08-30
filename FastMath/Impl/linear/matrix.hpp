/******************************************************************************
  文 件 名   : matrix.hpp
  作    者   : Barry
  生成日期   : 2024-08-26
  最近修改   :
  功能描述   : 矩阵运算、分解
  函数列表   :
                tmul        矩阵相乘
                inv         矩阵求逆
                ssgj        对称正定矩阵求逆
                sdet        计算行列式
                rank        计算秩
                chol        对称正定矩阵的Cholesky分解
                lluu        方阵LU分解
                maqr         QR分解，要求行数m大于列数n
                muav         SVD分解
                ginv         SVD伪逆
                strq         约化对称矩阵为对角阵
                sstq         求对角阵特征值
                hhbg         约化一般矩阵为上H矩阵
                hhqr         求上H矩阵特征值
                eig_jacobi   特征值分解 雅克比法
                eig_jcbj     特征值分解 雅克比过关法
                eig_top1     最大特征值

  修改历史   :
  1.日    期   : 2024-08-26
    作    者   : Barry
    修改内容   : 创建文件
******************************************************************************/

#ifndef FASTMATH_MATRIX_HPP
#define FASTMATH_MATRIX_HPP

#include <cmath>
#include <cstring>
#include <sstream>
#include <iostream>

#include "Impl/basic/complex.hpp"

#define MatRef(name)    &name[0][0]
#define FastMatRef(name)    &name._data[0][0]

namespace FastMath::Impl {
    inline double init(double p)       //实数初始化
    {
        p = 0.0;
        return (p);
    }

    inline complex init(complex p)     //复数初始化
    {
        p = complex(0.0, 0.0);
        return (p);
    }

    inline double ffabs(double p)        //计算实数的绝对值
    {
        double q;
        q = fabs(p);
        return (q);
    }

    inline double ffabs(complex p)      //计算复数的模
    {
        double q;
        q = p.cfabs();
        return (q);
    }

    inline double ff(double p)         //计算1.0/p
    {
        double q;
        q = 1.0 / p;
        return (q);
    }

    inline complex ff(const complex &p)       //计算(1+j0)/p
    {
        complex q;
        q = complex(1.0, 0.0) / p;
        return (q);
    }

    /**
     * 矩阵数组格式化为string类型
     * @param a a[m][n] 矩阵数组
     * @param m m 行数
     * @param n n 列数
     * @param w w 格式化时每个数据的占位宽度
     * @return 返回string类型的字符串
     */
    inline std::string MatToString(double a[], int m, int n, int w = 10) {
        std::stringstream ss;
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++)
                ss << std::setw(w) << a[i * n + j] << "    ";
            ss << std::endl;
        }
        return ss.str();
    }

    /**
     * 检查两个矩阵是否相等
     * @tparam m 矩阵行数
     * @tparam n 矩阵列数
     * @param a a[m][n]
     * @param b b[m][n]
     * @param eps 检查精度
     * @return 返回是否相等
     */

    template<size_t m, size_t n>
    inline bool isEqual(double a[], double b[], double eps = 1E-12)
    {
        bool eq = true;
        for (int i = 0; i < m; ++i)
            for (int j = 0; j < n; ++j)
                eq &= std::abs(a[i*n + j] - b[i*n + j]) <= eps;
        return eq;
    }

    /**
     * 矩阵数组拷贝, 将a拷贝到b
     * @tparam m 行数
     * @tparam n 列数
     * @tparam T 元素类型 double或complex
     * @param a a[m][n]
     * @param b b[m][n]
     */
    template<size_t m, size_t n, typename T>
    inline void MatCopy(T a[], T b[])
    {
        std::memcpy(b, a, sizeof(T)*m*n);
    }

    /**
     * 将矩阵设置为全0
     * @param a a[m][n]
     * @param m 行数
     * @param n 列数
     */
    inline void setZero(double a[], int m, int n)
    {
        std::memset(a, 0, sizeof(double)*m*n);
    }

    /**
     * 将矩阵对角线设为1
     * @param a a[m][n]
     * @param m 行数
     * @param n 列数
     */
    inline void setIdentity(double a[], int m, int n)
    {
        int num_d = m > n ? n : m;
        std::memset(a, 0, sizeof(double)*m*n);
        for (int i = 0; i < num_d; ++i)
                a[i*n + i] = 1.0;
    }

    /*  矩阵相乘
     *  a, ma, na       矩阵A[ma][na]
     *  b, mb, nb       矩阵B[mb][nb]
     *  c, ma, nb       乘积矩阵C[ma][nb] = A[ma][na] * B[mb][nb]
     * */
    template<typename T>
    //模板声明T为类型参数，用来兼容实数和复数计算
    inline void tmul(T a[], int ma, int na, T b[], int mb, int nb, T c[]) {
        int i, j, k, u;
        assert(na == mb);
        for (i = 0; i <= ma - 1; i++)
            for (j = 0; j <= nb - 1; j++) {
                u = i * nb + j;
                c[u] = init(c[u]);            //乘积矩阵元素初始化
                for (k = 0; k <= mb - 1; k++)
                    c[u] = c[u] + a[i * na + k] * b[k * nb + j];
            }
    }

    /*  矩阵求逆
     *  a   原矩阵。返回逆矩阵。
     *  N   方阵维度
     * */
    template<size_t N, class T>
    //模板声明T为类型参数
    inline int inv(T a[])       //若矩阵奇异，则返回标志值0，否则返回标志值非0。
    {
        constexpr int n = N;
        int i, j, k, l, u, v;
        int is[n], js[n];
        double d, q;
        T p;
        for (k = 0; k <= n - 1; k++) {
            d = 0.0;
            for (i = k; i <= n - 1; i++)      //选主元
                for (j = k; j <= n - 1; j++) {
                    l = i * n + j;
                    q = ffabs(a[l]);        //计算元素绝对值（模）
                    if (q > d) {
                        d = q;
                        is[k] = i;
                        js[k] = j;
                    }
                }
            if (d + 1.0 == 1.0)            //矩阵奇异
            {
                return (0);             //返回奇异标志值
            }
            if (is[k] != k)
                for (j = 0; j <= n - 1; j++)     //行交换
                {
                    u = k * n + j;
                    v = is[k] * n + j;
                    p = a[u];
                    a[u] = a[v];
                    a[v] = p;
                }
            if (js[k] != k)
                for (i = 0; i <= n - 1; i++)     //列交换
                {
                    u = i * n + k;
                    v = i * n + js[k];
                    p = a[u];
                    a[u] = a[v];
                    a[v] = p;
                }
            l = k * n + k;
            a[l] = ff(a[l]);            //计算1/a[l]
            for (j = 0; j <= n - 1; j++)      //归一化
                if (j != k) {
                    u = k * n + j;
                    a[u] = a[u] * a[l];
                }
            for (i = 0; i <= n - 1; i++)      //消元计算
                if (i != k)
                    for (j = 0; j <= n - 1; j++)
                        if (j != k) {
                            u = i * n + j;
                            a[u] = a[u] - a[i * n + k] * a[k * n + j];
                        }
            for (i = 0; i <= n - 1; i++)
                if (i != k) {
                    u = i * n + k;
                    a[u] = (a[u] - a[u] - a[u]) * a[l];
                }
        }
        for (k = n - 1; k >= 0; k--)      //恢复行列交换
        {
            if (js[k] != k)
                for (j = 0; j <= n - 1; j++) {
                    u = k * n + j;
                    v = js[k] * n + j;
                    p = a[u];
                    a[u] = a[v];
                    a[v] = p;
                }
            if (is[k] != k)
                for (i = 0; i <= n - 1; i++) {
                    u = i * n + k;
                    v = i * n + is[k];
                    p = a[u];
                    a[u] = a[v];
                    a[v] = p;
                }
        }
        return (1);
    }


    /*
     * 对称正定矩阵求逆
     *  a   原矩阵。返回逆矩阵。
     *  N   方阵维度
     *
     */
    template<size_t N>
    inline int ssgj(double a[]) {
        constexpr int n = N;
        int i, j, k, m;
        double w, g;
        double b[n];
        for (k = 0; k <= n - 1; k++) {
            w = a[0];
            if (fabs(w) + 1.0 == 1.0) {
                return (0);
            }
            m = n - k - 1;
            for (i = 1; i <= n - 1; i++) {
                g = a[i * n];
                b[i] = g / w;
                if (i <= m) b[i] = -b[i];
                for (j = 1; j <= i; j++)
                    a[(i - 1) * n + j - 1] = a[i * n + j] + g * b[j];
            }
            a[n * n - 1] = 1.0 / w;
            for (i = 1; i <= n - 1; i++)
                a[(n - 1) * n + i - 1] = b[i];
        }
        for (i = 0; i <= n - 2; i++)
            for (j = i + 1; j <= n - 1; j++)
                a[i * n + j] = a[j * n + i];
        return (1);
    }


    /*  行列式求值，高斯消取法变换得到上三角矩阵，然后计算主对角线乘积
     *  a[n][n]       存放方阵A的元素。返回时被破坏。
     *
     * */
    inline double sdet(double a[], int n) {
        int i, j, k, is, js, l, u, v;
        double f, det, q, d;
        f = 1.0;
        det = 1.0;
        for (k = 0; k <= n - 2; k++) {
            q = 0.0;
            for (i = k; i <= n - 1; i++)
                for (j = k; j <= n - 1; j++) {
                    l = i * n + j;
                    d = fabs(a[l]);
                    if (d > q) {
                        q = d;
                        is = i;
                        js = j;
                    }
                }
            if (q + 1.0 == 1.0) {
                det = 0.0;
                return (det);
            }
            if (is != k) {
                f = -f;
                for (j = k; j <= n - 1; j++) {
                    u = k * n + j;
                    v = is * n + j;
                    d = a[u];
                    a[u] = a[v];
                    a[v] = d;
                }
            }
            if (js != k) {
                f = -f;
                for (i = k; i <= n - 1; i++) {
                    u = i * n + js;
                    v = i * n + k;
                    d = a[u];
                    a[u] = a[v];
                    a[v] = d;
                }
            }
            l = k * n + k;
            det = det * a[l];
            for (i = k + 1; i <= n - 1; i++) {
                d = a[i * n + k] / a[l];
                for (j = k + 1; j <= n - 1; j++) {
                    u = i * n + j;
                    a[u] = a[u] - d * a[k * n + j];
                }
            }
        }
        det = f * det * a[n * n - 1];
        return (det);
    }


    /* 矩阵求秩，使用全选主元高斯消取法
     * a[m][n]      存放m×n阶矩阵A的元素。返回时将被破坏。
     * 函数返回A的秩。
     * */
    inline int rank(double a[], int m, int n) {
        int i, j, k, nn, is, js, l, ll, u, v;
        double q, d;
        nn = m;
        if (m >= n) nn = n;
        k = 0;
        for (l = 0; l <= nn - 1; l++) {
            q = 0.0;
            for (i = l; i <= m - 1; i++)
                for (j = l; j <= n - 1; j++) {
                    ll = i * n + j;
                    d = fabs(a[ll]);
                    if (d > q) {
                        q = d;
                        is = i;
                        js = j;
                    }
                }
            if (q + 1.0 == 1.0) return (k);
            k = k + 1;
            if (is != l) {
                for (j = l; j <= n - 1; j++) {
                    u = l * n + j;
                    v = is * n + j;
                    d = a[u];
                    a[u] = a[v];
                    a[v] = d;
                }
            }
            if (js != l) {
                for (i = l; i <= m - 1; i++) {
                    u = i * n + js;
                    v = i * n + l;
                    d = a[u];
                    a[u] = a[v];
                    a[v] = d;
                }
            }
            ll = l * n + l;
            for (i = l + 1; i <= n - 1; i++) {
                d = a[i * n + l] / a[ll];
                for (j = l + 1; j <= n - 1; j++) {
                    u = i * n + j;
                    a[u] = a[u] - d * a[l * n + j];
                }
            }
        }
        return (k);
    }


    /* 正定矩阵的Cholesky分解
     *  a[n][n]       存放对称正定矩阵A,会被破坏
     *  返回时其下三角部分存放分解得到的下三角阵L，其余元素均为0。
     *  函数返回标志值。若等于0，则表示失败；若大于0，则表示正常。
     *
     * */
    inline int chol(double a[], int n) {
        int i, j, k, u, l;
        if ((a[0] + 1.0 == 1.0) || (a[0] < 0.0)) {
            return (0);
        }
        a[0] = sqrt(a[0]);
        for (i = 1; i <= n - 1; i++) {
            u = i * n;
            a[u] = a[u] / a[0];
        }
        for (j = 1; j <= n - 1; j++) {
            l = j * n + j;
            for (k = 0; k <= j - 1; k++) {
                u = j * n + k;
                a[l] = a[l] - a[u] * a[u];
            }
            if ((a[l] + 1.0 == 1.0) || (a[l] < 0.0)) {
                return (0);
            }
            a[l] = sqrt(a[l]);
            for (i = j + 1; i <= n - 1; i++) {
                u = i * n + j;
                for (k = 0; k <= j - 1; k++)
                    a[u] = a[u] - a[i * n + k] * a[j * n + k];
                a[u] = a[u] / a[l];
            }
        }
        for (i = 0; i <= n - 2; i++)
            for (j = i + 1; j <= n - 1; j++) a[i * n + j] = 0.0;
        return (1);
    }

    /*
     *  矩阵的LU分解
     *  a[n][n]        存放n阶矩阵A。返回时存放Q矩阵。Q = L + U - I
     *  l[n][n]        返回时存放下三角矩阵L。
     *  u[n][n]        返回时存放上三角矩阵U。
     *  函数返回标志值。若为0，则表示失败；若不为0，则表示正常。
     *
     *  计算量 O(n^3), 本方法没有进行主元选取，因此数值不稳定
     * */

    inline int lluu(double a[], int n, double l[], double u[]) {
        int i, j, k, w, v, ll;
        for (k = 0; k <= n - 2; k++) {
            ll = k * n + k;
            if (fabs(a[ll]) + 1.0 == 1.0)   //右下角子阵中a[k][k]=0
            {
                return 0;
            }
            for (i = k + 1; i <= n - 1; i++) {
                w = i * n + k;
                a[w] = a[w] / a[ll];
            }
            for (i = k + 1; i <= n - 1; i++) {
                w = i * n + k;
                for (j = k + 1; j <= n - 1; j++) {
                    v = i * n + j;
                    a[v] = a[v] - a[w] * a[k * n + j];
                }
            }
        }
        for (i = 0; i <= n - 1; i++) {
            for (j = 0; j < i; j++)      //L-I
            {
                w = i * n + j;
                l[w] = a[w];
                u[w] = 0.0;
            }
            w = i * n + i;
            l[w] = 1.0;
            u[w] = a[w];
            for (j = i + 1; j <= n - 1; j++)    //U
            {
                w = i * n + j;
                l[w] = 0.0;
                u[w] = a[w];
            }
        }
        return 1;
    }


    /*
     *  QR分解。
     *  a[m][n]   存放m×n的实矩阵A。要求 m>=n。
     *  返回时其右上三角部分存放QR分解中的上三角阵R。
     *  q[m][m]   返回QR分解中的正交矩阵Q。
     *  函数返回标志值。若为0，则表示失败；若不为0，则表示正常。
     *
     *  使用Householder变换进行QR分解
     * */

    inline int maqr(double a[], int m, int n, double q[]) {
        int i, j, k, l, nn, p, jj;
        double u, alpha, w, t;
        if (m < n) {
            return 0;
        }
        for (i = 0; i <= m - 1; i++)
            for (j = 0; j <= m - 1; j++) {
                l = i * m + j;
                q[l] = 0.0;
                if (i == j) q[l] = 1.0;
            }
        nn = n;
        if (m == n) nn = m - 1;
        for (k = 0; k <= nn - 1; k++) {
            u = 0.0;
            l = k * n + k;
            for (i = k; i <= m - 1; i++) {
                w = fabs(a[i * n + k]);
                if (w > u) u = w;
            }
            alpha = 0.0;
            for (i = k; i <= m - 1; i++) {
                t = a[i * n + k] / u;
                alpha = alpha + t * t;
            }
            if (a[l] > 0.0) u = -u;
            alpha = u * sqrt(alpha);
            if (fabs(alpha) + 1.0 == 1.0) {
                return 0;
            }
            u = sqrt(2.0 * alpha * (alpha - a[l]));
            if ((u + 1.0) != 1.0) {
                a[l] = (a[l] - alpha) / u;
                for (i = k + 1; i <= m - 1; i++) {
                    p = i * n + k;
                    a[p] = a[p] / u;
                }
                for (j = 0; j <= m - 1; j++) {
                    t = 0.0;
                    for (jj = k; jj <= m - 1; jj++) t = t + a[jj * n + k] * q[jj * m + j];
                    for (i = k; i <= m - 1; i++) {
                        p = i * m + j;
                        q[p] = q[p] - 2.0 * t * a[i * n + k];
                    }
                }
                for (j = k + 1; j <= n - 1; j++) {
                    t = 0.0;
                    for (jj = k; jj <= m - 1; jj++) t = t + a[jj * n + k] * a[jj * n + j];
                    for (i = k; i <= m - 1; i++) {
                        p = i * n + j;
                        a[p] = a[p] - 2.0 * t * a[i * n + k];
                    }
                }
                a[l] = alpha;
                for (i = k + 1; i <= m - 1; i++) a[i * n + k] = 0.0;
            }
        }
        for (i = 0; i <= m - 2; i++)
            for (j = i + 1; j <= m - 1; j++) {
                p = i * m + j;
                l = j * m + i;
                t = q[p];
                q[p] = q[l];
                q[l] = t;
            }
        return 1;
    }


    /*
     * 奇异值分解, 使用Householder变换进行SVD分解, A = U * a * VT
     *
     * a[m][n]     存放m×n的实矩阵A。返回时其对角线给出奇异值（以非递增次序排列），其余元素均为0。
     * u[m][m]     返回左奇异向量U。
     * vt[n][n]     返回右奇异向量VT 。
     * eps         给定的精度要求。
     * ka          其值为max(m，n)＋1。
     * */

    inline void ppp(double a[], double e[], double s[], double vt[], int m, int n);

    inline void sss(double fg[2], double cs[2]);

    template<size_t m, size_t n>
    inline int muav(double a[], double u[], double vt[], double eps = 1E-8) {
        constexpr int ka = (m > n ? m : n) + 1;

        int i, j, k, l, it, ll, kk, ix, iy, mm, nn, iz, m1, ks;
        double d, dd, t, sm, sm1, em1, sk, ek, b, c, shh, fg[2], cs[2];

        double s[ka], e[ka], w[ka];
        it = 60;
        k = n;
        if (m - 1 < n) k = m - 1;
        l = m;
        if (n - 2 < m) l = n - 2;
        if (l < 0) l = 0;
        ll = k;
        if (l > k) ll = l;
        if (ll >= 1) {
            for (kk = 1; kk <= ll; kk++) {
                if (kk <= k) {
                    d = 0.0;
                    for (i = kk; i <= m; i++) {
                        ix = (i - 1) * n + kk - 1;
                        d = d + a[ix] * a[ix];
                    }
                    s[kk - 1] = sqrt(d);
                    if (s[kk - 1] != 0.0) {
                        ix = (kk - 1) * n + kk - 1;
                        if (a[ix] != 0.0) {
                            s[kk - 1] = fabs(s[kk - 1]);
                            if (a[ix] < 0.0) s[kk - 1] = -s[kk - 1];
                        }
                        for (i = kk; i <= m; i++) {
                            iy = (i - 1) * n + kk - 1;
                            a[iy] = a[iy] / s[kk - 1];
                        }
                        a[ix] = 1.0 + a[ix];
                    }
                    s[kk - 1] = -s[kk - 1];
                }
                if (n >= kk + 1) {
                    for (j = kk + 1; j <= n; j++) {
                        if ((kk <= k) && (s[kk - 1] != 0.0)) {
                            d = 0.0;
                            for (i = kk; i <= m; i++) {
                                ix = (i - 1) * n + kk - 1;
                                iy = (i - 1) * n + j - 1;
                                d = d + a[ix] * a[iy];
                            }
                            d = -d / a[(kk - 1) * n + kk - 1];
                            for (i = kk; i <= m; i++) {
                                ix = (i - 1) * n + j - 1;
                                iy = (i - 1) * n + kk - 1;
                                a[ix] = a[ix] + d * a[iy];
                            }
                        }
                        e[j - 1] = a[(kk - 1) * n + j - 1];
                    }
                }
                if (kk <= k) {
                    for (i = kk; i <= m; i++) {
                        ix = (i - 1) * m + kk - 1;
                        iy = (i - 1) * n + kk - 1;
                        u[ix] = a[iy];
                    }
                }
                if (kk <= l) {
                    d = 0.0;
                    for (i = kk + 1; i <= n; i++)
                        d = d + e[i - 1] * e[i - 1];
                    e[kk - 1] = sqrt(d);
                    if (e[kk - 1] != 0.0) {
                        if (e[kk] != 0.0) {
                            e[kk - 1] = fabs(e[kk - 1]);
                            if (e[kk] < 0.0) e[kk - 1] = -e[kk - 1];
                        }
                        for (i = kk + 1; i <= n; i++)
                            e[i - 1] = e[i - 1] / e[kk - 1];
                        e[kk] = 1.0 + e[kk];
                    }
                    e[kk - 1] = -e[kk - 1];
                    if ((kk + 1 <= m) && (e[kk - 1] != 0.0)) {
                        for (i = kk + 1; i <= m; i++) w[i - 1] = 0.0;
                        for (j = kk + 1; j <= n; j++)
                            for (i = kk + 1; i <= m; i++)
                                w[i - 1] = w[i - 1] + e[j - 1] * a[(i - 1) * n + j - 1];
                        for (j = kk + 1; j <= n; j++)
                            for (i = kk + 1; i <= m; i++) {
                                ix = (i - 1) * n + j - 1;
                                a[ix] = a[ix] - w[i - 1] * e[j - 1] / e[kk];
                            }
                    }
                    for (i = kk + 1; i <= n; i++)
                        vt[(i - 1) * n + kk - 1] = e[i - 1];
                }
            }
        }

        mm = n;
        if (m + 1 < n) mm = m + 1;
        if (k < n) s[k] = a[k * n + k];
        if (m < mm) s[mm - 1] = 0.0;
        if (l + 1 < mm) e[l] = a[l * n + mm - 1];
        e[mm - 1] = 0.0;
        nn = m;
        if (m > n) nn = n;
        if (nn >= k + 1) {
            for (j = k + 1; j <= nn; j++) {
                for (i = 1; i <= m; i++)
                    u[(i - 1) * m + j - 1] = 0.0;
                u[(j - 1) * m + j - 1] = 1.0;
            }
        }
        if (k >= 1) {
            for (ll = 1; ll <= k; ll++) {
                kk = k - ll + 1;
                iz = (kk - 1) * m + kk - 1;
                if (s[kk - 1] != 0.0) {
                    if (nn >= kk + 1)
                        for (j = kk + 1; j <= nn; j++) {
                            d = 0.0;
                            for (i = kk; i <= m; i++) {
                                ix = (i - 1) * m + kk - 1;
                                iy = (i - 1) * m + j - 1;
                                d = d + u[ix] * u[iy] / u[iz];
                            }
                            d = -d;
                            for (i = kk; i <= m; i++) {
                                ix = (i - 1) * m + j - 1;
                                iy = (i - 1) * m + kk - 1;
                                u[ix] = u[ix] + d * u[iy];
                            }
                        }
                    for (i = kk; i <= m; i++) {
                        ix = (i - 1) * m + kk - 1;
                        u[ix] = -u[ix];
                    }
                    u[iz] = 1.0 + u[iz];
                    if (kk - 1 >= 1)
                        for (i = 1; i <= kk - 1; i++)
                            u[(i - 1) * m + kk - 1] = 0.0;
                } else {
                    for (i = 1; i <= m; i++)
                        u[(i - 1) * m + kk - 1] = 0.0;
                    u[(kk - 1) * m + kk - 1] = 1.0;
                }
            }
        }
        for (ll = 1; ll <= n; ll++) {
            kk = n - ll + 1;
            iz = kk * n + kk - 1;
            if ((kk <= l) && (e[kk - 1] != 0.0)) {
                for (j = kk + 1; j <= n; j++) {
                    d = 0.0;
                    for (i = kk + 1; i <= n; i++) {
                        ix = (i - 1) * n + kk - 1;
                        iy = (i - 1) * n + j - 1;
                        d = d + vt[ix] * vt[iy] / vt[iz];
                    }
                    d = -d;
                    for (i = kk + 1; i <= n; i++) {
                        ix = (i - 1) * n + j - 1;
                        iy = (i - 1) * n + kk - 1;
                        vt[ix] = vt[ix] + d * vt[iy];
                    }
                }
            }
            for (i = 1; i <= n; i++)
                vt[(i - 1) * n + kk - 1] = 0.0;
            vt[iz - n] = 1.0;
        }
        for (i = 1; i <= m; i++)
            for (j = 1; j <= n; j++)
                a[(i - 1) * n + j - 1] = 0.0;
        m1 = mm;
        it = 60;
        while (1 == 1) {
            if (mm == 0) {
                ppp(a, e, s, vt, m, n);
                return (1);
            }
            if (it == 0) {
                ppp(a, e, s, vt, m, n);
                return (-1);
            }
            kk = mm - 1;
            while ((kk != 0) && (fabs(e[kk - 1]) != 0.0)) {
                d = fabs(s[kk - 1]) + fabs(s[kk]);
                dd = fabs(e[kk - 1]);
                if (dd > eps * d) kk = kk - 1;
                else e[kk - 1] = 0.0;
            }
            if (kk == mm - 1) {
                kk = kk + 1;
                if (s[kk - 1] < 0.0) {
                    s[kk - 1] = -s[kk - 1];
                    for (i = 1; i <= n; i++) {
                        ix = (i - 1) * n + kk - 1;
                        vt[ix] = -vt[ix];
                    }
                }
                while ((kk != m1) && (s[kk - 1] < s[kk])) {
                    d = s[kk - 1];
                    s[kk - 1] = s[kk];
                    s[kk] = d;
                    if (kk < n)
                        for (i = 1; i <= n; i++) {
                            ix = (i - 1) * n + kk - 1;
                            iy = (i - 1) * n + kk;
                            d = vt[ix];
                            vt[ix] = vt[iy];
                            vt[iy] = d;
                        }
                    if (kk < m)
                        for (i = 1; i <= m; i++) {
                            ix = (i - 1) * m + kk - 1;
                            iy = (i - 1) * m + kk;
                            d = u[ix];
                            u[ix] = u[iy];
                            u[iy] = d;
                        }
                    kk = kk + 1;
                }
                it = 60;
                mm = mm - 1;
            } else {
                ks = mm;
                while ((ks > kk) && (fabs(s[ks - 1]) != 0.0)) {
                    d = 0.0;
                    if (ks != mm) d = d + fabs(e[ks - 1]);
                    if (ks != kk + 1) d = d + fabs(e[ks - 2]);
                    dd = fabs(s[ks - 1]);
                    if (dd > eps * d) ks = ks - 1;
                    else s[ks - 1] = 0.0;
                }
                if (ks == kk) {
                    kk = kk + 1;
                    d = fabs(s[mm - 1]);
                    t = fabs(s[mm - 2]);
                    if (t > d) d = t;
                    t = fabs(e[mm - 2]);
                    if (t > d) d = t;
                    t = fabs(s[kk - 1]);
                    if (t > d) d = t;
                    t = fabs(e[kk - 1]);
                    if (t > d) d = t;
                    sm = s[mm - 1] / d;
                    sm1 = s[mm - 2] / d;
                    em1 = e[mm - 2] / d;
                    sk = s[kk - 1] / d;
                    ek = e[kk - 1] / d;
                    b = ((sm1 + sm) * (sm1 - sm) + em1 * em1) / 2.0;
                    c = sm * em1;
                    c = c * c;
                    shh = 0.0;
                    if ((b != 0.0) || (c != 0.0)) {
                        shh = sqrt(b * b + c);
                        if (b < 0.0) shh = -shh;
                        shh = c / (b + shh);
                    }
                    fg[0] = (sk + sm) * (sk - sm) - shh;
                    fg[1] = sk * ek;
                    for (i = kk; i <= mm - 1; i++) {
                        sss(fg, cs);
                        if (i != kk) e[i - 2] = fg[0];
                        fg[0] = cs[0] * s[i - 1] + cs[1] * e[i - 1];
                        e[i - 1] = cs[0] * e[i - 1] - cs[1] * s[i - 1];
                        fg[1] = cs[1] * s[i];
                        s[i] = cs[0] * s[i];
                        if ((cs[0] != 1.0) || (cs[1] != 0.0))
                            for (j = 1; j <= n; j++) {
                                ix = (j - 1) * n + i - 1;
                                iy = (j - 1) * n + i;
                                d = cs[0] * vt[ix] + cs[1] * vt[iy];
                                vt[iy] = -cs[1] * vt[ix] + cs[0] * vt[iy];
                                vt[ix] = d;
                            }
                        sss(fg, cs);
                        s[i - 1] = fg[0];
                        fg[0] = cs[0] * e[i - 1] + cs[1] * s[i];
                        s[i] = -cs[1] * e[i - 1] + cs[0] * s[i];
                        fg[1] = cs[1] * e[i];
                        e[i] = cs[0] * e[i];
                        if (i < m)
                            if ((cs[0] != 1.0) || (cs[1] != 0.0))
                                for (j = 1; j <= m; j++) {
                                    ix = (j - 1) * m + i - 1;
                                    iy = (j - 1) * m + i;
                                    d = cs[0] * u[ix] + cs[1] * u[iy];
                                    u[iy] = -cs[1] * u[ix] + cs[0] * u[iy];
                                    u[ix] = d;
                                }
                    }
                    e[mm - 2] = fg[0];
                    it = it - 1;
                } else {
                    if (ks == mm) {
                        kk = kk + 1;
                        fg[1] = e[mm - 2];
                        e[mm - 2] = 0.0;
                        for (ll = kk; ll <= mm - 1; ll++) {
                            i = mm + kk - ll - 1;
                            fg[0] = s[i - 1];
                            sss(fg, cs);
                            s[i - 1] = fg[0];
                            if (i != kk) {
                                fg[1] = -cs[1] * e[i - 2];
                                e[i - 2] = cs[0] * e[i - 2];
                            }
                            if ((cs[0] != 1.0) || (cs[1] != 0.0))
                                for (j = 1; j <= n; j++) {
                                    ix = (j - 1) * n + i - 1;
                                    iy = (j - 1) * n + mm - 1;
                                    d = cs[0] * vt[ix] + cs[1] * vt[iy];
                                    vt[iy] = -cs[1] * vt[ix] + cs[0] * vt[iy];
                                    vt[ix] = d;
                                }
                        }
                    } else {
                        kk = ks + 1;
                        fg[1] = e[kk - 2];
                        e[kk - 2] = 0.0;
                        for (i = kk; i <= mm; i++) {
                            fg[0] = s[i - 1];
                            sss(fg, cs);
                            s[i - 1] = fg[0];
                            fg[1] = -cs[1] * e[i - 1];
                            e[i - 1] = cs[0] * e[i - 1];
                            if ((cs[0] != 1.0) || (cs[1] != 0.0))
                                for (j = 1; j <= m; j++) {
                                    ix = (j - 1) * m + i - 1;
                                    iy = (j - 1) * m + kk - 2;
                                    d = cs[0] * u[ix] + cs[1] * u[iy];
                                    u[iy] = -cs[1] * u[ix] + cs[0] * u[iy];
                                    u[ix] = d;
                                }
                        }
                    }
                }
            }
        }
        return (1);
    }

    void ppp(double a[], double e[], double s[], double vt[], int m, int n) {
        int i, j, p, q;
        double d;
        if (m >= n) i = n;
        else i = m;
        for (j = 1; j <= i - 1; j++) {
            a[(j - 1) * n + j - 1] = s[j - 1];
            a[(j - 1) * n + j] = e[j - 1];
        }
        a[(i - 1) * n + i - 1] = s[i - 1];
        if (m < n) a[(i - 1) * n + i] = e[i - 1];
        for (i = 1; i <= n - 1; i++)
            for (j = i + 1; j <= n; j++) {
                p = (i - 1) * n + j - 1;
                q = (j - 1) * n + i - 1;
                d = vt[p];
                vt[p] = vt[q];
                vt[q] = d;
            }
        return;
    }

    void sss(double fg[2], double cs[2]) {
        double r, d;
        if ((fabs(fg[0]) + fabs(fg[1])) == 0.0) {
            cs[0] = 1.0;
            cs[1] = 0.0;
            d = 0.0;
        } else {
            d = sqrt(fg[0] * fg[0] + fg[1] * fg[1]);
            if (fabs(fg[0]) > fabs(fg[1])) {
                d = fabs(d);
                if (fg[0] < 0.0) d = -d;
            }
            if (fabs(fg[1]) >= fabs(fg[0])) {
                d = fabs(d);
                if (fg[1] < 0.0) d = -d;
            }
            cs[0] = fg[0] / d;
            cs[1] = fg[1] / d;
        }
        r = 1.0;
        if (fabs(fg[0]) > fabs(fg[1])) r = cs[1];
        else if (cs[0] != 0.0) r = 1.0 / cs[0];
        fg[0] = d;
        fg[1] = r;
        return;
    }

    /*
     *  矩阵广义逆, 使用SVD分解计算广义逆
     *  a[m]n]     存放m×n的实矩阵A。
     *              返回时其对角线给出奇异值（以非递增次序排列），其余元素均为0。
     *  aa[n][m]   返回A的广义逆  。
     *  eps        给定的精度要求。
     *  u[m][m]    返回左奇异向量U。
     *  vt[n][n]    返回右奇异向量VT。
     *  函数返回标志值。若小于0，则表示失败；若大于0，则表示正常。
     * */

    template<size_t m, size_t n>
    inline int ginv(double a[], double aa[], double u[], double vt[], double eps = 1E-10) {
        int i, j, k, l, t, p, q, f;
        i = muav<m, n>(a, u, vt, eps);
        if (i < 0) return (-1);
        j = n;
        if (m < n) j = m;
        j = j - 1;
        k = 0;
        while ((k <= j) && (a[k * n + k] != 0.0)) k = k + 1;
        k = k - 1;
        for (i = 0; i <= n - 1; i++)
            for (j = 0; j <= m - 1; j++) {
                t = i * m + j;
                aa[t] = 0.0;
                for (l = 0; l <= k; l++) {
                    f = l * n + i;
                    p = j * m + l;
                    q = l * n + l;
                    aa[t] = aa[t] + vt[f] * u[p] / a[q];
                }
            }
        return (1);
    }


    /*
     * 约化对称矩阵为对角阵, A =  Q * T * Q'
     * 用Householder变换将n阶是对称矩阵约化为对称三角阵
     *
     * a[n][n]        存放n阶实对称矩阵A
     * q[n][n]        返回豪斯荷尔德变换的乘积矩阵Q。 Q为正交矩阵
     * b[n]           返回对称三角阵T中的主对角线元素。
     * c[n]           前n-1个元素返回对称三角阵T中的次对角线元素。
     * 若矩阵非对称，则显示错误信息，并返回0标志值。否则返回非0标志值。
     *
     * */
    inline int strq(double a[], int n, double q[], double b[], double c[]) {
        int i, j, k, u;
        double h, f, g, h2;
        for (i = 0; i < n; i++)
            for (j = 0; j < i - 1; j++)
                if (a[i * n + j] != a[j * n + i]) {
                    return 0;
                }
        for (i = 0; i <= n - 1; i++)
            for (j = 0; j <= n - 1; j++) {
                u = i * n + j;
                q[u] = a[u];
            }
        for (i = n - 1; i >= 1; i--) {
            h = 0.0;
            if (i > 1)
                for (k = 0; k <= i - 1; k++) {
                    u = i * n + k;
                    h = h + q[u] * q[u];
                }
            if (h + 1.0 == 1.0) {
                c[i] = 0.0;
                if (i == 1) c[i] = q[i * n + i - 1];
                b[i] = 0.0;
            } else {
                c[i] = sqrt(h);
                u = i * n + i - 1;
                if (q[u] > 0.0) c[i] = -c[i];
                h = h - q[u] * c[i];
                q[u] = q[u] - c[i];
                f = 0.0;
                for (j = 0; j <= i - 1; j++) {
                    q[j * n + i] = q[i * n + j] / h;
                    g = 0.0;
                    for (k = 0; k <= j; k++) g = g + q[j * n + k] * q[i * n + k];
                    if (j + 1 <= i - 1)
                        for (k = j + 1; k <= i - 1; k++) g = g + q[k * n + j] * q[i * n + k];
                    c[j] = g / h;
                    f = f + g * q[j * n + i];
                }
                h2 = f / (h + h);
                for (j = 0; j <= i - 1; j++) {
                    f = q[i * n + j];
                    g = c[j] - h2 * f;
                    c[j] = g;
                    for (k = 0; k <= j; k++) {
                        u = j * n + k;
                        q[u] = q[u] - f * c[k] - g * q[i * n + k];
                    }
                }
                b[i] = h;
            }
        }
        for (i = 0; i <= n - 2; i++) c[i] = c[i + 1];
        c[n - 1] = 0.0;
        b[0] = 0.0;
        for (i = 0; i <= n - 1; i++) {
            if ((b[i] != 0.0) && (i - 1 >= 0))
                for (j = 0; j <= i - 1; j++) {
                    g = 0.0;
                    for (k = 0; k <= i - 1; k++) g = g + q[i * n + k] * q[k * n + j];
                    for (k = 0; k <= i - 1; k++) {
                        u = k * n + j;
                        q[u] = q[u] - g * q[k * n + i];
                    }
                }
            u = i * n + i;
            b[i] = q[u];
            q[u] = 1.0;
            if (i - 1 >= 0)
                for (j = 0; j <= i - 1; j++) {
                    q[i * n + j] = 0.0;
                    q[j * n + i] = 0.0;
                }
        }
        return 1;
    }


    /*
     *  求对称三对角阵的特征值
     *  b[n]       存放n阶实对称三角阵的主对角线上的元素。返回时存放全部特征值。
     *  c[n]       前n-1个元素存放实对称三角阵的次对角线上的元素。
     *  q[n][n]    若存放n阶单位矩阵，则返回实对称三对角阵T的特征向量组；
     *             若存放函数strq()所返回的一般实对称矩阵A的豪斯荷尔德变换的乘积矩阵Q，则返回实对称矩阵A的特征向量组。
     *             其中q中的第j列为与数组b中第j个特征值对应的特征向量。
     *
     *  max_it      为最大迭代次数
     * 返回的标志值为迭代次数。
     * */

    inline int sstq(int n, double b[], double c[], double q[], double eps = 1E-8, int max_it = 100) {
        int i, j, k, m, it, u, v;
        double d, f, h, g, p, r, e, s;
        c[n - 1] = 0.0;
        d = 0.0;
        f = 0.0;
        for (j = 0; j <= n - 1; j++) {
            it = 0;
            h = eps * (fabs(b[j]) + fabs(c[j]));
            if (h > d) d = h;
            m = j;
            while ((m <= n - 1) && (fabs(c[m]) > d)) m = m + 1;
            if (m != j) {
                do {
                    if (it == max_it) {
                        return (it);
                    }
                    it = it + 1;
                    g = b[j];
                    p = (b[j + 1] - g) / (2.0 * c[j]);
                    r = sqrt(p * p + 1.0);
                    if (p >= 0.0) b[j] = c[j] / (p + r);
                    else b[j] = c[j] / (p - r);
                    h = g - b[j];
                    for (i = j + 1; i <= n - 1; i++) b[i] = b[i] - h;
                    f = f + h;
                    p = b[m];
                    e = 1.0;
                    s = 0.0;
                    for (i = m - 1; i >= j; i--) {
                        g = e * c[i];
                        h = e * p;
                        if (fabs(p) >= fabs(c[i])) {
                            e = c[i] / p;
                            r = sqrt(e * e + 1.0);
                            c[i + 1] = s * p * r;
                            s = e / r;
                            e = 1.0 / r;
                        } else {
                            e = p / c[i];
                            r = sqrt(e * e + 1.0);
                            c[i + 1] = s * c[i] * r;
                            s = 1.0 / r;
                            e = e / r;
                        }
                        p = e * b[i] - s * g;
                        b[i + 1] = h + s * (e * g + s * b[i]);
                        for (k = 0; k <= n - 1; k++) {
                            u = k * n + i + 1;
                            v = u - 1;
                            h = q[u];
                            q[u] = s * q[v] + e * h;
                            q[v] = e * q[v] - s * h;
                        }
                    }
                    c[j] = s * p;
                    b[j] = e * p;
                } while (fabs(c[j]) > d);
            }
            b[j] = b[j] + f;
        }
        for (i = 0; i <= n - 1; i++) {
            k = i;
            p = b[i];
            if (i + 1 <= n - 1) {
                j = i + 1;
                while ((j <= n - 1) && (b[j] <= p)) {
                    k = j;
                    p = b[j];
                    j = j + 1;
                }
            }
            if (k != i) {
                b[k] = b[i];
                b[i] = p;
                for (j = 0; j <= n - 1; j++) {
                    u = j * n + i;
                    v = j * n + k;
                    p = q[u];
                    q[u] = q[v];
                    q[v] = p;
                }
            }
        }
        return (it);
    }


    /*
     * 约化一般实矩阵为上H阵, 用初等变换将一般实矩阵约化为上H矩阵，即Hessenberg矩阵
     *
     * a[n][n]    一般实矩阵。返回上H矩阵。
     * */

    inline void hhbg(double a[], int n) {
        int i, j, k, u, v;
        double d, t;
        for (k = 1; k <= n - 2; k++) {
            d = 0.0;
            for (j = k; j <= n - 1; j++) {
                u = j * n + k - 1;
                t = a[u];
                if (fabs(t) > fabs(d)) {
                    d = t;
                    i = j;
                }
            }
            if (fabs(d) + 1.0 != 1.0) {
                if (i != k) {
                    for (j = k - 1; j <= n - 1; j++) {
                        u = i * n + j;
                        v = k * n + j;
                        t = a[u];
                        a[u] = a[v];
                        a[v] = t;
                    }
                    for (j = 0; j <= n - 1; j++) {
                        u = j * n + i;
                        v = j * n + k;
                        t = a[u];
                        a[u] = a[v];
                        a[v] = t;
                    }
                }
                for (i = k + 1; i <= n - 1; i++) {
                    u = i * n + k - 1;
                    t = a[u] / d;
                    a[u] = 0.0;
                    for (j = k; j <= n - 1; j++) {
                        v = i * n + j;
                        a[v] = a[v] - t * a[k * n + j];
                    }
                    for (j = 0; j <= n - 1; j++) {
                        v = j * n + k;
                        a[v] = a[v] + t * a[j * n + i];
                    }
                }
            }
        }
        return;
    }


    /*
     * 求上H矩阵特征值的QR方法,
     * a[n][n]      上H矩阵
     * u[n]         返回n个特征值的实部
     * v[n]         返回n个特征值的虚部
     * eps          精度要求
     *
     * */

    inline int hhqr(double a[], int n, double u[], double v[], double eps = 1E-8, int max_it = 100) {
        int jt = max_it;
        int m, it, i, j, k, l, ii, jj, kk, ll;
        double b, c, w, g, xy, p, q, r, x, s, e, f, z, y;
        it = 0;
        m = n;
        while (m != 0) {
            l = m - 1;
            while ((l > 0) && (fabs(a[l * n + l - 1]) > eps *
                                                        (fabs(a[(l - 1) * n + l - 1]) + fabs(a[l * n + l]))))
                l = l - 1;
            ii = (m - 1) * n + m - 1;
            jj = (m - 1) * n + m - 2;
            kk = (m - 2) * n + m - 1;
            ll = (m - 2) * n + m - 2;
            if (l == m - 1) {
                u[m - 1] = a[(m - 1) * n + m - 1];
                v[m - 1] = 0.0;
                m = m - 1;
                it = 0;
            } else if (l == m - 2) {
                b = -(a[ii] + a[ll]);
                c = a[ii] * a[ll] - a[jj] * a[kk];
                w = b * b - 4.0 * c;
                y = sqrt(fabs(w));
                if (w > 0.0)     //计算两个实特征值
                {
                    xy = 1.0;
                    if (b < 0.0) xy = -1.0;
                    u[m - 1] = (-b - xy * y) / 2.0;
                    u[m - 2] = c / u[m - 1];
                    v[m - 1] = 0.0;
                    v[m - 2] = 0.0;
                } else           //计算复特征值
                {
                    u[m - 1] = -b / 2.0;
                    u[m - 2] = u[m - 1];
                    v[m - 1] = y / 2.0;
                    v[m - 2] = -v[m - 1];
                }
                m = m - 2;
                it = 0;
            } else {
                if (it >= jt)               //超过最大迭代次数
                {
                    return 0;
                }
                it = it + 1;
                for (j = l + 2; j <= m - 1; j++) a[j * n + j - 2] = 0.0;
                for (j = l + 3; j <= m - 1; j++) a[j * n + j - 3] = 0.0;
                for (k = l; k <= m - 2; k++) {
                    if (k != l) {
                        p = a[k * n + k - 1];
                        q = a[(k + 1) * n + k - 1];
                        r = 0.0;
                        if (k != m - 2) r = a[(k + 2) * n + k - 1];
                    } else {
                        x = a[ii] + a[ll];
                        y = a[ll] * a[ii] - a[kk] * a[jj];
                        ii = l * n + l;
                        jj = l * n + l + 1;
                        kk = (l + 1) * n + l;
                        ll = (l + 1) * n + l + 1;
                        p = a[ii] * (a[ii] - x) + a[jj] * a[kk] + y;
                        q = a[kk] * (a[ii] + a[ll] - x);
                        r = a[kk] * a[(l + 2) * n + l + 1];
                    }
                    if ((fabs(p) + fabs(q) + fabs(r)) != 0.0) {
                        xy = 1.0;
                        if (p < 0.0) xy = -1.0;
                        s = xy * sqrt(p * p + q * q + r * r);
                        if (k != l) a[k * n + k - 1] = -s;
                        e = -q / s;
                        f = -r / s;
                        x = -p / s;
                        y = -x - f * r / (p + s);
                        g = e * r / (p + s);
                        z = -x - e * q / (p + s);
                        for (j = k; j <= m - 1; j++) {
                            ii = k * n + j;
                            jj = (k + 1) * n + j;
                            p = x * a[ii] + e * a[jj];
                            q = e * a[ii] + y * a[jj];
                            r = f * a[ii] + g * a[jj];
                            if (k != m - 2) {
                                kk = (k + 2) * n + j;
                                p = p + f * a[kk];
                                q = q + g * a[kk];
                                r = r + z * a[kk];
                                a[kk] = r;
                            }
                            a[jj] = q;
                            a[ii] = p;
                        }
                        j = k + 3;
                        if (j >= m - 1) j = m - 1;
                        for (i = l; i <= j; i++) {
                            ii = i * n + k;
                            jj = i * n + k + 1;
                            p = x * a[ii] + e * a[jj];
                            q = e * a[ii] + y * a[jj];
                            r = f * a[ii] + g * a[jj];
                            if (k != m - 2) {
                                kk = i * n + k + 2;
                                p = p + f * a[kk];
                                q = q + g * a[kk];
                                r = r + z * a[kk];
                                a[kk] = r;
                            }
                            a[jj] = q;
                            a[ii] = p;
                        }
                    }
                }
            }
        }
        return 1;
    }

    /*
     * jacobi法求实对称矩阵的特征值分解
     * a[n][n]        实对称矩阵。对角线元素返回特征值。
     * v[n][n]        返回特征向量
     * eps            精度要求
     * mat_it         最大迭代次数
     *
     * */

    template<size_t n>
    inline int eig_jacobi(double a[], double v[], double eps = 1E-8, int mat_it = 200) {
        int i, j, p, q, u, w, t, s, count;
        double fm, cn, sn, omega, x, y, d;
        for (i = 0; i < n; i++)
            for (j = i + 1; j < n; j++)
                if (a[i * n + j] != a[j * n + i]) {
                    // 矩阵不对称
                    return 0;
                }
        for (i = 0; i <= n - 1; i++) {
            v[i * n + i] = 1.0;
            for (j = 0; j <= n - 1; j++)
                if (i != j) v[i * n + j] = 0.0;
        }
        count = 1;
        while (count <= mat_it) {
            fm = 0.0;
            for (i = 1; i <= n - 1; i++)
                for (j = 0; j <= i - 1; j++) {
                    d = fabs(a[i * n + j]);
                    if ((i != j) && (d > fm)) {
                        fm = d;
                        p = i;
                        q = j;
                    }
                }
            if (fm < eps) return (count);
            count = count + 1;
            u = p * n + q;
            w = p * n + p;
            t = q * n + p;
            s = q * n + q;
            x = -a[u];
            y = (a[s] - a[w]) / 2.0;
            omega = x / sqrt(x * x + y * y);
            if (y < 0.0) omega = -omega;
            sn = 1.0 + sqrt(1.0 - omega * omega);
            sn = omega / sqrt(2.0 * sn);
            cn = sqrt(1.0 - sn * sn);
            fm = a[w];
            a[w] = fm * cn * cn + a[s] * sn * sn + a[u] * omega;
            a[s] = fm * sn * sn + a[s] * cn * cn - a[u] * omega;
            a[u] = 0.0;
            a[t] = 0.0;
            for (j = 0; j <= n - 1; j++)
                if ((j != p) && (j != q)) {
                    u = p * n + j;
                    w = q * n + j;
                    fm = a[u];
                    a[u] = fm * cn + a[w] * sn;
                    a[w] = -fm * sn + a[w] * cn;
                }
            for (i = 0; i <= n - 1; i++)
                if ((i != p) && (i != q)) {
                    u = i * n + p;
                    w = i * n + q;
                    fm = a[u];
                    a[u] = fm * cn + a[w] * sn;
                    a[w] = -fm * sn + a[w] * cn;
                }
            for (i = 0; i <= n - 1; i++) {
                u = i * n + p;
                w = i * n + q;
                fm = v[u];
                v[u] = fm * cn + v[w] * sn;
                v[w] = -fm * sn + v[w] * cn;
            }
        }
        return (count);
    }

    /*
     *  jacobi过关法求是对称矩阵的特征值分解
     *  a[n][n]     实对称矩阵。对角线元素返回特征值。
     *  v[n][n]     返回特征向量
     *  eps         精度要求
     *  若矩阵不对称，则显示错误信息，并返回0标志值。
     *
     *
     * */
    template<size_t n>
    inline int eig_jcbj(double a[], double v[], double eps = 1E-8) {
        int i, j, p, q, u, w, t, s;
        double ff, fm, cn, sn, omega, x, y, d;
        for (i = 0; i < n; i++)
            for (j = i + 1; j < n; j++)
                if (a[i * n + j] != a[j * n + i]) {
                    // 矩阵不对称
                    return 0;
                }
        for (i = 0; i <= n - 1; i++)      //特征向量初始化
        {
            v[i * n + i] = 1.0;
            for (j = 0; j <= n - 1; j++)
                if (i != j) v[i * n + j] = 0.0;
        }
        ff = 0.0;
        for (i = 1; i <= n - 1; i++)
            for (j = 0; j <= i - 1; j++) {
                d = a[i * n + j];
                ff = ff + d * d;
            }
        ff = sqrt(2.0 * ff);
        ff = ff / (1.0 * n);
        while (ff >= eps) {
            d = 0.0;
            for (i = 1; (i <= n - 1) && (d <= ff); i++)
                for (j = 0; (j <= i - 1) && (d <= ff); j++) {
                    d = fabs(a[i * n + j]);
                    p = i;
                    q = j;
                }
            if (d <= ff) ff = ff / (1.0 * n);
            else {
                u = p * n + q;
                w = p * n + p;
                t = q * n + p;
                s = q * n + q;
                x = -a[u];
                y = (a[s] - a[w]) / 2.0;
                omega = x / sqrt(x * x + y * y);
                if (y < 0.0) omega = -omega;
                sn = 1.0 + sqrt(1.0 - omega * omega);
                sn = omega / sqrt(2.0 * sn);
                cn = sqrt(1.0 - sn * sn);
                fm = a[w];
                a[w] = fm * cn * cn + a[s] * sn * sn + a[u] * omega;
                a[s] = fm * sn * sn + a[s] * cn * cn - a[u] * omega;
                a[u] = 0.0;
                a[t] = 0.0;
                for (j = 0; j <= n - 1; j++)
                    if ((j != p) && (j != q)) {
                        u = p * n + j;
                        w = q * n + j;
                        fm = a[u];
                        a[u] = fm * cn + a[w] * sn;
                        a[w] = -fm * sn + a[w] * cn;
                    }
                for (i = 0; i <= n - 1; i++)
                    if ((i != p) && (i != q)) {
                        u = i * n + p;
                        w = i * n + q;
                        fm = a[u];
                        a[u] = fm * cn + a[w] * sn;
                        a[w] = -fm * sn + a[w] * cn;
                    }
                for (i = 0; i <= n - 1; i++) {
                    u = i * n + p;
                    w = i * n + q;
                    fm = v[u];
                    v[u] = fm * cn + v[w] * sn;
                    v[w] = -fm * sn + v[w] * cn;
                }
            }
        }
        return 1;
    }


    /* 乘幂法计算最大特征值即特征向量
     * a[n][n]        实矩阵
     * v[n]           特征向量
     * eps            精度要求
     * 函数返回绝对值最大的特征值。
     *
     * */
    template<size_t n>
    double eig_top1(double a[], double v[], double eps = 1E-8, int mat_it = 1000) {
        int i, j, k, flag = 1, iteration;
        double lambda, sum, z, err, t, d, f;
        double u[n];
        iteration = 0;
        do {
            iteration++;
            for (i = 0; i < n; i++)       //计算u=Av
            {
                sum = 0.0;
                for (j = 0; j < n; j++) { sum = sum + a[i * n + j] * v[j]; }
                u[i] = sum;
            }
            d = 0.0;            //计算向量的范数
            for (k = 0; k < n; k++) d = d + u[k] * u[k];
            d = sqrt(d);
            for (i = 0; i < n; i++) { v[i] = u[i] / d; }
            if (iteration > 1) {
                err = fabs((d - t) / d);
                f = 1;
                if (v[0] * z < 0) f = -1;
                if (err < eps) { flag = 0; }
            }
            if (flag == 1) {
                t = d;
                z = v[0];
            }
            if (iteration >= mat_it) flag = 0;
        } while (flag == 1);
        lambda = f * d;
        return (lambda);
    }

}

#endif //FASTMATH_MATRIX_HPP
