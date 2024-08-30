/******************************************************************************
  文 件 名   : function.hpp
  作    者   : Barry
  生成日期   : 2024-08-29
  最近修改   :
  功能描述   : 特殊函数计算
  函数列表   :
  修改历史   :
  1.日    期   : 2024-08-29
    作    者   : Barry
    修改内容   : 创建文件
******************************************************************************/

#ifndef FASTMATH2_FUNCTION_HPP
#define FASTMATH2_FUNCTION_HPP

#include <cmath>
#include "Impl/Common.hpp"

namespace FastMath::Impl::util {
    /**
     * Gamma函数
     * @param x 自变量值。要求x>0。
     * @return  函数返回Gamma函数值, 若x<=0, 返回-1
     */
    inline double gamma(double x) {
        int i;
        double y, t, s, u;
        double a[11] = {0.0000677106, -0.0003442342,
                        0.0015397681, -0.0024467480, 0.0109736958,
                        -0.0002109075, 0.0742379071, 0.0815782188,
                        0.4118402518, 0.4227843370, 1.0};
        if (x <= 0.0) {
            return (-1.0);
        }
        y = x;
        if (y <= 1.0) {
            t = 1.0 / (y * (y + 1.0));
            y = y + 2.0;
        } else if (y <= 2.0) {
            t = 1.0 / y;
            y = y + 1.0;
        } else if (y <= 3.0) t = 1.0;
        else {
            t = 1.0;
            while (y > 3.0) {
                y = y - 1.0;
                t = t * y;
            }
        }
        s = a[0];
        u = y - 2.0;
        for (i = 1; i <= 10; i++) s = s * u + a[i];
        s = s * t;
        return (s);
    }

    /**
     * Beta函数， 利用Gamma函数计算Beta函数
     * @param x
     * @param y
     * @return
     */

    inline double beta(double x, double y) {
        return gamma(x) * gamma(y) / gamma(x + y);
    }


    /**
     * 不完全Gamma函数
     * @param a 自变量值。要求 a>0。若a<=0 返回-1
     * @param x 自变量值。要求 x>=0。若a<0 返回-1
     * @return  函数返回不完全Gamma函数值。
     */
    inline double ingamma(double a, double x) {
        int n, flag;
        double p, q, d, s, s1, p0, q0, p1, q1, qq;
        if ((a <= 0.0) || (x < 0.0)) {
            return (-1.0);
        }
        if (x + 1.0 == 1.0) return (0.0);
        if (x > 1.0e+35) return (1.0);
        q = log(x);
        q = a * q;
        qq = exp(q);
        if (x < 1.0 + a) {
            p = a;
            d = 1.0 / a;
            s = d;
            n = 0;
            do {
                n = n + 1;
                p = p + 1.0;
                d = d * x / p;
                s = s + d;
                flag = (fabs(d) >= fabs(s) * 1.0e-10);
            } while ((n <= 100) && (flag));
            if (!flag) {
                s = s * exp(-x) * qq / gamma(a);
                return (s);
            }
        } else {
            s = 1.0 / x;
            p0 = 0.0;
            p1 = 1.0;
            q0 = 1.0;
            q1 = x;
            for (n = 1; n <= 100; n++) {
                p0 = p1 + (n - a) * p0;
                q0 = q1 + (n - a) * q0;
                p = x * p0 + n * p1;
                q = x * q0 + n * q1;
                if (fabs(q) + 1.0 != 1.0) {
                    s1 = p / q;
                    p1 = p;
                    q1 = q;
                    if (fabs((s1 - s) / s1) < 1.0e-10) {
                        s = s1 * exp(-x) * qq / gamma(a);
                        return (1.0 - s);
                    }
                    s = s1;
                }
                p1 = p;
                q1 = q;
            }
        }
        s = 1.0 - s * exp(-x) * qq / gamma(a);
        return (s);
    }


    /**
     * 误差函数
     * @param x 自变量值
     * @return  函数返回误差函数值
     */
    inline double errf(double x) {
        double y;
        if (x >= 0.0) y = ingamma(0.5, x * x);
        else y = -ingamma(0.5, x * x);
        return (y);
    }

    /**
     *  正态分布函数, 表示概率，也就是高斯函数的积分
     * @param a 数学期望值。
     * @param d 标准差 d*d为方差值。要求d>0。
     * @param x 随机变量值。
     * @return  函数返回正态分布函数值（0~1之间，表示概率）。
     */
    inline double gass(double a, double d, double x) {
        double y;
        if (d <= 0.0) d = 1.0e-10;
        y = 0.5 + 0.5 * errf((x - a) / (sqrt(2.0) * d));
        return (y);
    }


    /**
     * 正弦积分,即sin(t)/t从0到x的积分
     * @param x 自变量值。
     * @return  函数返回正弦积分值。
     */
    inline double sinn(double x) {
        int m, i, j;
        double s, p, ep, h, aa, bb, w, xx, g;
        double t[5] = {-0.9061798459, -0.5384693101, 0.0,
                       0.5384693101, 0.9061798459};
        double c[5] = {0.2369268851, 0.4786286705, 0.5688888889,
                       0.4786286705, 0.2369268851};
        m = 1;
        if (x == 0) return (0.0);
        h = fabs(x);
        s = fabs(0.0001 * h);
        p = 1.0e+35;
        ep = 0.000001;
        g = 0.0;
        while ((ep >= 0.0000001) && (fabs(h) > s)) {
            g = 0.0;
            for (i = 1; i <= m; i++) {
                aa = (i - 1.0) * h;
                bb = i * h;
                w = 0.0;
                for (j = 0; j <= 4; j++) {
                    xx = ((bb - aa) * t[j] + (bb + aa)) / 2.0;
                    w = w + sin(xx) / xx * c[j];
                }
                g = g + w;
            }
            g = g * h / 2.0;
            ep = fabs(g - p) / (1.0 + fabs(g));
            p = g;
            m = m + 1;
            h = fabs(x) / m;
        }
        return (g);
    }

    /**
     * 余弦积分, 也就是 -cos(t)/t 从x到正无穷的积分
     * @param x 自变量值。
     * @return  函数返回余弦积分值。
     */
    inline double coss(double x) {
        int m, i, j;
        double s, p, ep, h, aa, bb, w, xx, g, r, q;
        double t[5] = {-0.9061798459, -0.5384693101, 0.0,
                       0.5384693101, 0.9061798459};
        double c[5] = {0.2369268851, 0.4786286705, 0.5688888889,
                       0.4786286705, 0.2369268851};
        m = 1;
        if (x == 0) x = 1.0e-35;
        if (x < 0.0) x = -x;
        r = 0.57721566490153286060651;
        q = r + log(x);
        h = x;
        s = fabs(0.0001 * h);
        p = 1.0e+35;
        ep = 0.000001;
        g = 0.0;
        while ((ep >= 0.0000001) && (fabs(h) > s)) {
            g = 0.0;
            for (i = 1; i <= m; i++) {
                aa = (i - 1.0) * h;
                bb = i * h;
                w = 0.0;
                for (j = 0; j <= 4; j++) {
                    xx = ((bb - aa) * t[j] + (bb + aa)) / 2.0;
                    w = w + (1.0 - cos(xx)) / xx * c[j];
                }
                g = g + w;
            }
            g = g * h / 2.0;
            ep = fabs(g - p) / (1.0 + fabs(g));
            p = g;
            m = m + 1;
            h = x / m;
        }
        g = q - g;
        return (g);
    }


    /**
     * 指数积分, 要求x>0。  -e^(-t)/t从x到正无穷的积分
     * @param x 自变量值。
     * @return  函数返回指数积分值。
     */
    inline double expp(double x) {
        int m, i, j;
        double s, p, ep, h, aa, bb, w, xx, g, r, q;
        double t[5] = {-0.9061798459, -0.5384693101, 0.0,
                       0.5384693101, 0.9061798459};
        double c[5] = {0.2369268851, 0.4786286705, 0.5688888889,
                       0.4786286705, 0.2369268851};
        m = 1;
        if (x == 0) x = 1.0e-10;
        if (x < 0.0) x = -x;
        r = 0.57721566490153286060651;
        q = r + log(x);
        h = x;
        s = fabs(0.0001 * h);
        p = 1.0e+35;
        ep = 0.000001;
        g = 0.0;
        while ((ep >= 0.0000001) && (fabs(h) > s)) {
            g = 0.0;
            for (i = 1; i <= m; i++) {
                aa = (i - 1.0) * h;
                bb = i * h;
                w = 0.0;
                for (j = 0; j <= 4; j++) {
                    xx = ((bb - aa) * t[j] + (bb + aa)) / 2.0;
                    w = w + (exp(-xx) - 1.0) / xx * c[j];
                }
                g = g + w;
            }
            g = g * h / 2.0;
            ep = fabs(g - p) / (1.0 + fabs(g));
            p = g;
            m = m + 1;
            h = x / m;
        }
        g = q + g;
        return (g);
    }


}


#endif //FASTMATH2_FUNCTION_HPP
