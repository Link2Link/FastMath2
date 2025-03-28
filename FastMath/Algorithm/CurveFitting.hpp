/******************************************************************************
文 件 名   : CurveFitting.hpp
  作    者   : Barry
  生成日期   : 2025-01-13
  最近修改   :
  功能描述   : 各种曲线拟合
  函数列表   :
  修改历史   :
  1.日    期   : 2025-01-13
    作    者   : Barry
    修改内容   : RBF径向基拟合
  2.日    期   : 2025-03-28
    作    者   : Barry
    修改内容   : 三次样条插值 CubicSpline
******************************************************************************/

#ifndef FastMath2_CURVEFITTING_HPP
#define FastMath2_CURVEFITTING_HPP


#include "FastMath.hpp"
#include <vector>
#include <stdexcept>
#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <algorithm>
#include <iostream>

namespace FastMath::Algorithm {
    template<size_t M, size_t N>
    FastMath::Matrix<double, M, N> rbf(FastMath::Vector<double, M> c, FastMath::Vector<double, N> x, double delta) {
        FastMath::Matrix<double, M, N> r;
        double delta2 = delta * delta;
        for (size_t i = 0; i < M; i++) {
            for (size_t j = 0; j < N; j++) {
                r(i, j) = std::exp(-(c(i) - x(j)) * (c(i) - x(j)) / delta2 / 2);
            }
        }
        return r;
    }

    template<size_t N>
    FastMath::Vector<double, N> RBFtrain(FastMath::Vector<double, N> xc, FastMath::Vector<double, N> yc) {
        double d = (xc.max() - xc.min()) / N;


        auto Phi = rbf(xc, xc, d);
        auto w = Phi.pinv() * yc;

        return w;
    }

    template<size_t N>
    double RBFeval(double x, FastMath::Vector<double, N> xc, FastMath::Vector<double, N> w) {
        double d = (xc.max() - xc.min()) / N;
        auto base = rbf(xc, FastMath::Vector<double, 1>(x), d); // 计算预测点的 RBF 基函数值
        auto y = base.transpose() * w;
        return y(0);
    }


    /**
     *  CubicSpline的构造函数中进行三次样条计算，其中使用了vector的push_back，
     *  在实时系统中调用此构造函数可能会触发内存分配进而影响实时系统的实时性。
     *
     */
    template<size_t n>
    class CubicSpline {
    private:
        fm::Vector<double, n> x, y;
        fm::Vector<double, n - 1> a, b, c, d;

    public:
        CubicSpline(const fm::Vector<double, n> &x, const fm::Vector<double, n> &y) {
            if (n < 2) {
                throw std::invalid_argument("At least two points are required");
            }
            static_assert(n >= 2, "At least 2 points are required");
            for (int i = 0; i < n - 1; ++i) {
                if (x[i] >= x[i + 1]) {
                    throw std::invalid_argument("x must be strictly increasing");
                }
            }
            this->x = x;
            this->y = y;

            std::vector<double> h(n - 1);
            for (int i = 0; i < n - 1; ++i) {
                h[i] = x[i + 1] - x[i];
            }

            if (n == 2) {
                a[0] = y[0];
                b[0] = (y[1] - y[0]) / h[0];
                c[0] = 0.0;
                d[0] = 0.0;
                return;
            }

            constexpr int m = n - 2;
            std::vector<double> rhs(m);
            for (int i = 1; i <= m; ++i) {
                double delta_i = (y[i + 1] - y[i]) / h[i];
                double delta_im1 = (y[i] - y[i - 1]) / h[i - 1];
                rhs[i - 1] = 6 * (delta_i - delta_im1);
            }

            std::vector<double> a_diag;
            std::vector<double> b_diag(m);
            std::vector<double> c_diag;

            for (int k = 0; k < m; ++k) {
                b_diag[k] = 2 * (h[k] + h[k + 1]);
                if (k < m - 1) {
                    c_diag.push_back(h[k + 1]);
                }
                if (k > 0) {
                    a_diag.push_back(h[k]);
                }
            }

            std::vector<double> M_part(m);
            if (!solveThomas(a_diag, b_diag, c_diag, rhs, M_part)) {
                throw std::runtime_error("Failed to solve the tridiagonal system");
            }

            std::vector<double> M(n, 0.0);
            for (int k = 0; k < m; ++k) {
                M[k + 1] = M_part[k];
            }

            for (int i = 0; i < n - 1; ++i) {
                double h_i = h[i];
                double M_i = M[i];
                double M_ip1 = M[i + 1];
                a[i] = y[i];
                c[i] = M_i / 2.0;
                d[i] = (M_ip1 - M_i) / (6.0 * h_i);
                b[i] = (y[i + 1] - y[i]) / h_i - h_i * (2 * M_i + M_ip1) / 6.0;
            }
        }

        double interpolate(double t) const {
            if (t < x[0]) return y[0];
            if (t > x[n - 1]) return y[n - 1];

            int i = 0;
            for (int k = 1; k < n; ++k) {
                if (t < x[k]) {
                    i = k - 1; // 找到所在区间
                    break;
                }
            }

            if (i < 0) i = 0;
            else if (i >= n - 1) i = n - 2;
            double dx = t - x[i];
            return a[i] + b[i] * dx + c[i] * dx * dx + d[i] * dx * dx * dx;
        }

    private:
        bool solveThomas(const std::vector<double> &a, const std::vector<double> &b, const std::vector<double> &c,
                         const std::vector<double> &d, std::vector<double> &x) {
            int m = b.size();
            if (a.size() != m - 1 || c.size() != m - 1 || d.size() != m) return false;
            x.resize(m);
            if (m == 0) return true;

            std::vector<double> c_prime(m), d_prime(m);
            c_prime[0] = c[0] / b[0];
            d_prime[0] = d[0] / b[0];

            for (int i = 1; i < m; ++i) {
                double denominator = b[i] - a[i - 1] * c_prime[i - 1];
                if (denominator == 0) return false;
                c_prime[i] = (i < m - 1) ? c[i] / denominator : 0;
                d_prime[i] = (d[i] - a[i - 1] * d_prime[i - 1]) / denominator;
            }

            x[m - 1] = d_prime[m - 1];
            for (int i = m - 2; i >= 0; --i)
                x[i] = d_prime[i] - c_prime[i] * x[i + 1];

            return true;
        }
    };
}


#endif //FastMath2_CURVEFITTING_HPP
