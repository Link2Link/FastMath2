/******************************************************************************
  文 件 名   : complex.hpp
  作    者   : Barry
  生成日期   : 2024-08-26
  最近修改   :
  功能描述   : 
  函数列表   :
  修改历史   :
  1.日    期   : 2024-08-26
    作    者   : Barry
    修改内容   : 创建文件
******************************************************************************/

#ifndef FASTMATH_COMPLEX_HPP
#define FASTMATH_COMPLEX_HPP

#include <cmath>
#include <cstring>
#include <sstream>
#include "Impl/Common.hpp"

namespace FastMath::Impl {
    class complex {
    private:
        double R;
        double I;
    public:
        explicit complex(double real = 0, double imag = 0) : R{real}, I(imag) {}

        //输出形式为:(实部, 虚部)
        [[nodiscard]] std::string toString() const {
            std::stringstream ss;
            ss << "(" << R << ", " << I << ")";
            return ss.str();
        }

        //复数模
        [[nodiscard]] double cfabs() const {
            return std::sqrt(R * R + I * I);
        }

        //复数幅角
        [[nodiscard]] double angle() const {
            return std::atan2(I, R);
        }

        complex operator+(const complex &c2) const         //复数加法
        {
            complex c;
            c.R = R + c2.R;
            c.I = I + c2.I;
            return c;
        }

        complex &operator=(const complex &a) = default;

        complex &operator+=(const complex &c2) {
            R += c2.R;
            I += c2.I;
            return *this;
        }

        complex &operator-=(const complex &c2) {
            R -= c2.R;
            I -= c2.I;
            return *this;
        }


        complex operator-(const complex &c2) const         //复数减法
        {
            complex c;
            c.R = R - c2.R;
            c.I = I - c2.I;
            return c;
        }

        complex operator*(const complex &c2) const            //复数乘法
        {
            complex c;
            double p, q, s;
            p = R * c2.R;
            q = I * c2.I;
            s = (R + I) * (c2.R + c2.I);
            c.R = p - q;
            c.I = s - p - q;
            return c;
        }

        complex operator/(const complex &c2) const           //复数除法
        {
            complex c;
            double p, q, s, w;
            p = R * c2.R;
            q = -I * c2.I;
            s = (R + I) * (c2.R - c2.I);
            w = (c2.R) * (c2.R) + (c2.I) * (c2.I);
            if (w + 1.0 != 1.0) {
                c.R = (p - q) / w;
                c.I = (s - p - q) / w;
            } else {
                c.R = 1e+300;
                c.I = 1e+300;
            }
            return c;
        }

        [[nodiscard]] complex cpower(int n) const              //复数乘幂
        {
            complex c;
            double r, q;
            q = std::atan2(I, R);
            r = std::sqrt(R * R + I * I);
            if (r + 1.0 != 1.0) {
                r = n * std::log(r);
                r = std::exp(r);
            }
            c.R = r * std::cos(n * q);
            c.I = r * std::sin(n * q);
            return c;
        }

        void croot(int n, complex *p) const                 //复数的n次方根
        {
            complex c;
            int k;
            double r, q, t;
            if (n < 1) return;
            q = std::atan2(I, R);
            r = std::sqrt(R * R + I * I);
            if (r + 1.0 != 1.0) {
                r = (1.0 / n) * std::log(r);
                r = std::exp(r);
            }
            for (k = 0; k < n; k++) {
                t = (2.0 * k * M_PI_F + q) / n;
                c.R = r * std::cos(t);
                c.I = r * std::sin(t);
                p[k] = c;
            }
        }

        [[nodiscard]] complex cexp() const                      //复数指数
        {
            complex c;
            double p;
            p = std::exp(R);
            c.R = p * std::cos(I);
            c.I = p * std::sin(I);
            return c;
        }

        [[nodiscard]] complex clog() const                       //复数对数
        {
            complex c;
            double p;
            p = R * R + I * I;
            p = std::log(std::sqrt(p));
            c.R = p;
            c.I = std::atan2(I, R);
            return c;
        }

        [[nodiscard]] complex csin() const                       //复数正弦
        {
            complex c;
            double p, q;
            p = std::exp(I);
            q = std::exp(-I);
            c.R = std::sin(R) * (p + q) / 2;
            c.I = std::cos(R) * (p - q) / 2;
            return c;
        }

        [[nodiscard]] complex ccos() const                       //复数余弦
        {
            complex c;
            double p, q;
            p = std::exp(I);
            q = std::exp(-I);
            c.R = std::cos(R) * (p + q) / 2;
            c.I = -std::sin(R) * (p - q) / 2;
            return c;
        }


    };
}


#endif //FASTMATH_COMPLEX_HPP
