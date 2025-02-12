/******************************************************************************
  文 件 名   : filter.hpp
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

#ifndef FASTMATH2_FILTER_HPP
#define FASTMATH2_FILTER_HPP

#include <cmath>
#include "Impl/Common.hpp"

namespace FastMath::Impl::util {

    /**
     * Fourier级数逼近
     * @param n 等距点数为2n+1。
     * @param f[2n+1] 存放区间[0,2*3.1415926]内的2n+1个等距点处的函数值。
     * @param a[n+1] 返回Fourier级数系数。
     * @param b[n+1] 返回Fourier级数系数。
     *
     * 生成的傅里叶级数为 f(x) = 1/2*a_0 + sum_k( a_k cos kx + b_k sin kx)
     */

    inline void four(int n, double f[], double a[], double b[]) {
        int i, j;
        double t, c, s, c1, s1, u1, u2, u0;
        t = M_TWOPI_F / (2.0 * n + 1.0);
        c = cos(t);
        s = sin(t);
        t = 2.0 / (2.0 * n + 1.0);
        c1 = 1.0;
        s1 = 0.0;
        for (i = 0; i <= n; i++) {
            u1 = 0.0;
            u2 = 0.0;
            for (j = 2 * n; j >= 1; j--) {
                u0 = f[j] + 2.0 * c1 * u1 - u2;
                u2 = u1;
                u1 = u0;
            }
            a[i] = t * (f[0] + u1 * c1 - u2);
            b[i] = t * u1 * s1;
            u0 = c * c1 - s * s1;
            s1 = c * s1 + s * c1;
            c1 = u0;
        }
        return;
    }


    /**
     * 快速Fourier变换
     * @param n     采样点数 n = 2^k
     * @param k     满足 n = 2^k
     * @param pr[n]    flag = 0 : 存放采样输入的实部，返回变换的模;    flag = 1 : 存放傅里叶变换的n个实部， 返回傅里叶变换的模
     * @param pi[n]    flag = 0 : 存放采样输入的虚部，返回变换的幅角;    flag = 1 : 存放傅里叶变换的n个虚部， 返回傅里叶变换的幅角
     * @param fr[n]    flag = 0 : 返回变换(或逆变换)的实部; flag = 1 : 返回傅里叶逆变换的实部
     * @param fi[n]    flag = 0 : 返回变换(或逆变换)的虚部; flag = 1 : 返回傅里叶逆变换的虚部
     * @param flag  存放标志。flag=0表示作变换；flag=1表示作逆变换。
     */

    inline void kfft(int n, int k, double pr[], double pi[],
                     double fr[], double fi[], int flag) {
        int it, m, is, i, j, nv, kk;
        double p, q, s, vr, vi, poddr, poddi;
        for (it = 0; it <= n - 1; it++) {
            m = it;
            is = 0;
            for (i = 0; i <= k - 1; i++) {
                j = m / 2;
                is = 2 * is + (m - 2 * j);
                m = j;
            }
            fr[it] = pr[is];
            fi[it] = pi[is];
        }
        pr[0] = 1.0;
        pi[0] = 0.0;
        p = M_TWOPI_F / (1.0 * n);
        pr[1] = cos(p);
        pi[1] = -sin(p);
        if (flag != 0) pi[1] = -pi[1];      //逆变换
        for (i = 2; i <= n - 1; i++) {
            p = pr[i - 1] * pr[1];
            q = pi[i - 1] * pi[1];
            s = (pr[i - 1] + pi[i - 1]) * (pr[1] + pi[1]);
            pr[i] = p - q;
            pi[i] = s - p - q;
        }
        for (it = 0; it <= n - 2; it = it + 2) {
            vr = fr[it];
            vi = fi[it];
            fr[it] = vr + fr[it + 1];
            fi[it] = vi + fi[it + 1];
            fr[it + 1] = vr - fr[it + 1];
            fi[it + 1] = vi - fi[it + 1];
        }
        m = n / 2;
        nv = 2;
        for (kk = k - 2; kk >= 0; kk--) {
            m = m / 2;
            nv = 2 * nv;
            for (it = 0; it <= (m - 1) * nv; it = it + nv)
                for (j = 0; j <= (nv / 2) - 1; j++) {
                    p = pr[m * j] * fr[it + j + nv / 2];
                    q = pi[m * j] * fi[it + j + nv / 2];
                    s = pr[m * j] + pi[m * j];
                    s = s * (fr[it + j + nv / 2] + fi[it + j + nv / 2]);
                    poddr = p - q;
                    poddi = s - p - q;
                    fr[it + j + nv / 2] = fr[it + j] - poddr;
                    fi[it + j + nv / 2] = fi[it + j] - poddi;
                    fr[it + j] = fr[it + j] + poddr;
                    fi[it + j] = fi[it + j] + poddi;
                }
        }
        if (flag != 0)      //逆变换
            for (i = 0; i <= n - 1; i++) {
                fr[i] = fr[i] / (1.0 * n);
                fi[i] = fi[i] / (1.0 * n);
            }
        for (i = 0; i <= n - 1; i++)     //计算变换的模与幅角
        {
            pr[i] = std::sqrt(fr[i] * fr[i] + fi[i] * fi[i]);
            if (fabs(fr[i]) < 1E-8 * fabs(fi[i])) {
                if ((fi[i] * fr[i]) > 0) pi[i] = 90.0;
                else pi[i] = -90.0;
            } else pi[i] = atan2(fi[i], fr[i]) * 360.0 / M_TWOPI_F;
        }
        return;
    }


    /**
     * 快速Walsh变换
     * @param n         输入序列的长度。
     * @param k         满足 n = 2^k
     * @param p[n]      存放长度为n的给定输入序列。
     * @param x[n]     返回给定输入序列的Walsh变换序列。
     */
    inline void kfwt(int n, int k, double p[], double x[]) {
        int m, l, it, ii, i, j, is;
        double q;
        m = 1;
        l = n;
        it = 2;
        x[0] = 1;
        ii = n / 2;
        x[ii] = 2;
        for (i = 1; i <= k - 1; i++) {
            m = m + m;
            l = l / 2;
            it = it + it;
            for (j = 0; j <= m - 1; j++) x[j * l + l / 2] = it + 1 - x[j * l];
        }
        for (i = 0; i <= n - 1; i++) {
            ii = (int) (x[i] - 1);
            x[i] = p[ii];
        }
        l = 1;
        for (i = 1; i <= k; i++) {
            m = n / (2 * l) - 1;
            for (j = 0; j <= m; j++) {
                it = 2 * l * j;
                for (is = 0; is <= l - 1; is++) {
                    q = x[it + is] + x[it + is + l];
                    x[it + is + l] = x[it + is] - x[it + is + l];
                    x[it + is] = q;
                }
            }
            l = 2 * l;
        }
        return;
    }


    /**
     * Kalman滤波
     * @tparam n        动态系统的维数。
     * @tparam m        观测系统的维数
     * @tparam k        观测序列长度。
     * @param f[n][n]    系统状态转移矩阵
     * @param q[n][n]    模型噪声W的协方差阵。
     * @param r[m][m]    观测噪声V的协方差阵。
     * @param h[m][n]    观测矩阵
     * @param y[k][m]    观测向量序列。
     * @param x[k][n]    x[0][j]存放初值。其余各行返回状态向量估值序列。
     * @param p[n][n]    存放初值。返回最后时刻的估计误差协方差阵。
     * @param g[n][m]    返回最后时刻的稳定增益矩阵。
     * @return
     */
    template<size_t n, size_t m, size_t k>
    inline int kalman(double f[], double q[], double r[],
                      double h[], double y[], double x[], double p[], double g[]) {
        int i, j, kk, ii, jj, js;
        double e[m * m];
        constexpr int l = m > n ? m : n;
        double a[l * l];
        double b[l * l];
        for (i = 0; i <= n - 1; i++)
            for (j = 0; j <= n - 1; j++) {
                ii = i * l + j;
                a[ii] = 0.0;
                for (kk = 0; kk <= n - 1; kk++)
                    a[ii] = a[ii] + p[i * n + kk] * f[j * n + kk];
            }
        for (i = 0; i <= n - 1; i++)
            for (j = 0; j <= n - 1; j++) {
                ii = i * n + j;
                p[ii] = q[ii];
                for (kk = 0; kk <= n - 1; kk++)
                    p[ii] = p[ii] + f[i * n + kk] * a[kk * l + j];
            }
        for (ii = 2; ii <= k; ii++) {
            for (i = 0; i <= n - 1; i++)
                for (j = 0; j <= m - 1; j++) {
                    jj = i * l + j;
                    a[jj] = 0.0;
                    for (kk = 0; kk <= n - 1; kk++)
                        a[jj] = a[jj] + p[i * n + kk] * h[j * n + kk];
                }
            for (i = 0; i <= m - 1; i++)
                for (j = 0; j <= m - 1; j++) {
                    jj = i * m + j;
                    e[jj] = r[jj];
                    for (kk = 0; kk <= n - 1; kk++)
                        e[jj] = e[jj] + h[i * n + kk] * a[kk * l + j];
                }
            js = inv<m>(e);
            if (js == 0) {
                return (js);
            }
            for (i = 0; i <= n - 1; i++)
                for (j = 0; j <= m - 1; j++) {
                    jj = i * m + j;
                    g[jj] = 0.0;
                    for (kk = 0; kk <= m - 1; kk++)
                        g[jj] = g[jj] + a[i * l + kk] * e[j * m + kk];
                }
            for (i = 0; i <= n - 1; i++) {
                jj = (ii - 1) * n + i;
                x[jj] = 0.0;
                for (j = 0; j <= n - 1; j++)
                    x[jj] = x[jj] + f[i * n + j] * x[(ii - 2) * n + j];
            }
            for (i = 0; i <= m - 1; i++) {
                jj = i * l;
                b[jj] = y[(ii - 1) * m + i];
                for (j = 0; j <= n - 1; j++)
                    b[jj] = b[jj] - h[i * n + j] * x[(ii - 1) * n + j];
            }
            for (i = 0; i <= n - 1; i++) {
                jj = (ii - 1) * n + i;
                for (j = 0; j <= m - 1; j++)
                    x[jj] = x[jj] + g[i * m + j] * b[j * l];
            }
            if (ii < k) {
                for (i = 0; i <= n - 1; i++)
                    for (j = 0; j <= n - 1; j++) {
                        jj = i * l + j;
                        a[jj] = 0.0;
                        for (kk = 0; kk <= m - 1; kk++)
                            a[jj] = a[jj] - g[i * m + kk] * h[kk * n + j];
                        if (i == j) a[jj] = 1.0 + a[jj];
                    }
                for (i = 0; i <= n - 1; i++)
                    for (j = 0; j <= n - 1; j++) {
                        jj = i * l + j;
                        b[jj] = 0.0;
                        for (kk = 0; kk <= n - 1; kk++)
                            b[jj] = b[jj] + a[i * l + kk] * p[kk * n + j];
                    }
                for (i = 0; i <= n - 1; i++)
                    for (j = 0; j <= n - 1; j++) {
                        jj = i * l + j;
                        a[jj] = 0.0;
                        for (kk = 0; kk <= n - 1; kk++)
                            a[jj] = a[jj] + b[i * l + kk] * f[j * n + kk];
                    }
                for (i = 0; i <= n - 1; i++)
                    for (j = 0; j <= n - 1; j++) {
                        jj = i * n + j;
                        p[jj] = q[jj];
                        for (kk = 0; kk <= n - 1; kk++)
                            p[jj] = p[jj] + f[i * n + kk] * a[j * l + kk];
                    }
            }
        }
        return (js);
    }


    /**
     * alpha_beta_gemma滤波器
     * @tparam n    量测数据点数。
     * @param x[n]  n个等间隔点上的量测值。
     * @param t     采样周期。
     * @param a     滤波器结构参数Alpha。
     * @param b     滤波器结构参数Beta。
     * @param c     滤波器结构参数Gamma。
     * @param y     返回n个等间隔点上的滤波估值。
     * @param ss    滤波器初始位置
     * @param vv    滤波器初始速度
     * @param aa    滤波器初始加速度
     * @param se    指针返回结束时的位置
     * @param ve    指针返回结束时的速度
     * @param ae    指针返回结束时的加速度
     */

    template<size_t n>
    inline void kabg(double x[], double t, double a, double b, double c, double y[],
                     double ss = 0.0, double vv = 0.0, double aa = 0.0,
                     double* se = nullptr, double* ve = nullptr, double* ae= nullptr) {
        int i;
        double s1, v1, a1;
        for (i = 0; i <= n - 1; i++) {
            s1 = ss + t * vv + t * t * aa / 2.0;
            v1 = vv + t * aa;
            a1 = aa;
            ss = s1 + a * (x[i] - s1);
            y[i] = ss;
            vv = v1 + b * (x[i] - s1);
            aa = a1 + 2.0 * c * (x[i] - s1) / (t * t);
        }
        if (se != nullptr)
            *se = ss;
        if (ve != nullptr)
            *ve = vv;
        if (ae != nullptr)
            *ae = aa;
        return;
    }

}

#endif //FASTMATH2_FILTER_HPP
