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
******************************************************************************/

#ifndef FastMath2_CURVEFITTING_HPP
#define FastMath2_CURVEFITTING_HPP


#include "FastMath.hpp"

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
        auto base = rbf(xc, FastMath::Vector<double, 1>(x), d);  // 计算预测点的 RBF 基函数值
        auto y = base.transpose() * w;
        return y(0);
    }
}


#endif //FastMath2_CURVEFITTING_HPP
