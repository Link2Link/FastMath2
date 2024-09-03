/******************************************************************************
  文 件 名   : Dual.hpp
  作    者   : Barry
  生成日期   : 2024-09-02
  最近修改   :
  功能描述   : 
  函数列表   :
  修改历史   :
  1.日    期   : 2024-09-02
    作    者   : Barry
    修改内容   : 创建文件
******************************************************************************/

#ifndef FASTMATH2_DUAL_HPP
#define FASTMATH2_DUAL_HPP

#include "Core/Matrix.hpp"

namespace FastMath{
    template <typename Scalar, size_t N>
    struct Dual
    {
        static constexpr size_t WIDTH = N;
        Scalar value {};
        Vector<Scalar, N> derivative;

        Dual() = default;

        explicit Dual(Scalar v, size_t inputDimension = 128)
        {
            value = v;
            if (inputDimension < N) {
                derivative(inputDimension) = Scalar(1);
            }
        }

        explicit Dual(Scalar v, const Vector<Scalar, N>& d) :
                value(v), derivative(d)
        {}

        Dual<Scalar, N>& operator=(const Scalar& a)
        {
            derivative.setZero();
            value = a;
            return *this;
        }

        Dual<Scalar, N>& operator +=(const Dual<Scalar, N>& a)
        {
            return (*this = *this + a);
        }

        Dual<Scalar, N>& operator *=(const Dual<Scalar, N>& a)
        {
            return (*this = *this * a);
        }

        Dual<Scalar, N>& operator -=(const Dual<Scalar, N>& a)
        {
            return (*this = *this - a);
        }

        Dual<Scalar, N>& operator /=(const Dual<Scalar, N>& a)
        {
            return (*this = *this / a);
        }

        Dual<Scalar, N>& operator +=(Scalar a)
        {
            return (*this = *this + a);
        }

        Dual<Scalar, N>& operator -=(Scalar a)
        {
            return (*this = *this - a);
        }

        Dual<Scalar, N>& operator *=(Scalar a)
        {
            return (*this = *this * a);
        }

        Dual<Scalar, N>& operator /=(Scalar a)
        {
            return (*this = *this / a);
        }

    };

    // operators

    template <typename Scalar, size_t N>
    Dual<Scalar, N> operator+(const Dual<Scalar, N>& a)
    {
        return a;
    }

    template <typename Scalar, size_t N>
    Dual<Scalar, N> operator-(const Dual<Scalar, N>& a)
    {
        return Dual<Scalar, N>(-a.value, -a.derivative);
    }


    template <typename Scalar, size_t N>
    Dual<Scalar, N> operator+(const Dual<Scalar, N>& a, const Dual<Scalar, N>& b)
    {
        return Dual<Scalar, N>(a.value + b.value, a.derivative + b.derivative);
    }

    template <typename Scalar, size_t N>
    Dual<Scalar, N> operator-(const Dual<Scalar, N>& a, const Dual<Scalar, N>& b)
    {
        return a + (-b);
    }

    template <typename Scalar, size_t N>
    Dual<Scalar, N> operator+(const Dual<Scalar, N>& a, Scalar b)
    {
        return Dual<Scalar, N>(a.value + b, a.derivative);
    }

    template <typename Scalar, size_t N>
    Dual<Scalar, N> operator-(const Dual<Scalar, N>& a, Scalar b)
    {
        return a + (-b);
    }

    template <typename Scalar, size_t N>
    Dual<Scalar, N> operator+(Scalar a, const Dual<Scalar, N>& b)
    {
        return Dual<Scalar, N>(a + b.value, b.derivative);
    }

    template <typename Scalar, size_t N>
    Dual<Scalar, N> operator-(Scalar a, const Dual<Scalar, N>& b)
    {
        return a + (-b);
    }

    template <typename Scalar, size_t N>
    Dual<Scalar, N> operator*(const Dual<Scalar, N>& a, const Dual<Scalar, N>& b)
    {
        return Dual<Scalar, N>(a.value * b.value, a.value * b.derivative + b.value * a.derivative);
    }

    template <typename Scalar, size_t M, size_t N>
    FastMath::Vector<Dual<Scalar, N>, M> operator*(const FastMath::Matrix<Scalar, M, N> A, const FastMath::Vector<Dual<Scalar, N>, N>& x)
    {
        FastMath::Vector<Dual<Scalar, N>, M> res;
        for (size_t i = 0; i < M; ++i){
            res(i).value = 0;
            for (size_t j = 0; j < N; ++j)
            {
                res(i) += (x(j) * A(i, j));
            }

        }
        return res;
    }

    template <typename Scalar1, typename Scalar2, size_t N>
    Dual<Scalar1, N> operator*(const Dual<Scalar1, N>& a, Scalar2 b)
    {
        return Dual<Scalar1, N>(a.value * Scalar1(b), a.derivative * Scalar1(b));
    }

    template <typename Scalar1, typename Scalar2, size_t N>
    Dual<Scalar2, N> operator*(Scalar1 a, const Dual<Scalar2, N>& b)
    {
        return b * Scalar2(a);
    }

    template <typename Scalar, size_t N>
    Dual<Scalar, N> operator/(const Dual<Scalar, N>& a, const Dual<Scalar, N>& b)
    {
        const Scalar inv_b_real = Scalar(1) / b.value;
        return Dual<Scalar, N>(a.value * inv_b_real, a.derivative * inv_b_real -
                                                     a.value * b.derivative * inv_b_real * inv_b_real);
    }
    template <typename Scalar, size_t N>
    Dual<Scalar, N> operator/(const Dual<Scalar, N>& a, Scalar b)
    {
        return a * (Scalar(1) / b);
    }

    template <typename Scalar, size_t N>
    Dual<Scalar, N> operator/(Scalar a, const Dual<Scalar, N>& b)
    {
        const Scalar inv_b_real = Scalar(1) / b.value;
        return Dual<Scalar, N>(a * inv_b_real, (-inv_b_real * a * inv_b_real) * b.derivative);
    }

    // basic math
    // sqrt
    template <typename Scalar, size_t N>
    Dual<Scalar, N> sqrt(Dual<Scalar, N>& a)
    {
        Scalar real = std::sqrt(a.value);
        return Dual<Scalar, N>(real, a.derivative * (Scalar(1) / (Scalar(2) * real)));
    }

    // abs
    template <typename Scalar, size_t N>
    Dual<Scalar, N> abs(Dual<Scalar, N>& a)
    {
        return a.value >= Scalar(0) ? a : -a;
    }

    // ceil
    template <typename Scalar, size_t N>
    Dual<Scalar, N> ceil(Dual<Scalar, N>& a)
    {
        return Dual<Scalar, N>(std::ceil(a.value));
    }

    // floor
    template <typename Scalar, size_t N>
    Dual<Scalar, N> floor(Dual<Scalar, N>& a)
    {
        return Dual<Scalar, N>(std::floor(a.value));
    }

    // fmod
    template <typename Scalar, size_t N>
    Dual<Scalar, N> fmod(Dual<Scalar, N>& a, Scalar mod)
    {
        return Dual<Scalar, N>(a.value - Scalar(size_t(a.value / mod)) * mod, a.derivative);
    }

    // max
    template <typename Scalar, size_t N>
    Dual<Scalar, N> max(Dual<Scalar, N>& a, Dual<Scalar, N>& b)
    {
        return a.value >= b.value ? a : b;
    }

    // min
    template <typename Scalar, size_t N>
    Dual<Scalar, N> min(Dual<Scalar, N>& a, Dual<Scalar, N>& b)
    {
        return a.value < b.value ? a : b;
    }

    // isnan
    template <typename Scalar>
    bool IsNan(Scalar a)
    {
        return std::isnan(a);
    }

    template <typename Scalar, size_t N>
    bool IsNan(const Dual<Scalar, N>& a)
    {
        return IsNan(a.value);
    }

    // isfinite
    template <typename Scalar>
    bool IsFinite(Scalar a)
    {
        return isfinite(a);
    }

    template <typename Scalar, size_t N>
    bool IsFinite(Dual<Scalar, N>& a)
    {
        return IsFinite(a.value);
    }

    // isinf
    template <typename Scalar>
    bool IsInf(Scalar a)
    {
        return isinf(a);
    }

    template <typename Scalar, size_t N>
    bool IsInf(const Dual<Scalar, N>& a)
    {
        return IsInf(a.value);
    }

    // trig
    // sin
    template <typename Scalar, size_t N>
    Dual<Scalar, N> sin(Dual<Scalar, N>& a)
    {
        return Dual<Scalar, N>(std::sin(a.value), std::cos(a.value) * a.derivative);
    }

    // cos
    template <typename Scalar, size_t N>
    Dual<Scalar, N> cos(Dual<Scalar, N>& a)
    {
        return Dual<Scalar, N>(std::cos(a.value), -std::sin(a.value) * a.derivative);
    }

    // tan
    template <typename Scalar, size_t N>
    Dual<Scalar, N> tan(Dual<Scalar, N>& a)
    {
        Scalar real = std::tan(a.value);
        return Dual<Scalar, N>(real, (Scalar(1) + real * real) * a.derivative);
    }

    // asin
    template <typename Scalar, size_t N>
    Dual<Scalar, N> asin(Dual<Scalar, N>& a)
    {
        Scalar asin_d = Scalar(1) / std::sqrt(Scalar(1) - a.value * a.value);
        return Dual<Scalar, N>(std::asin(a.value), asin_d * a.derivative);
    }

    // acos
    template <typename Scalar, size_t N>
    Dual<Scalar, N> acos(Dual<Scalar, N>& a)
    {
        Scalar acos_d = -Scalar(1) / std::sqrt(Scalar(1) - a.value * a.value);
        return Dual<Scalar, N>(std::acos(a.value), acos_d * a.derivative);
    }

    // atan
    template <typename Scalar, size_t N>
    Dual<Scalar, N> atan(Dual<Scalar, N>& a)
    {
        Scalar atan_d = Scalar(1) / (Scalar(1) + a.value * a.value);
        return Dual<Scalar, N>(std::atan(a.value), atan_d * a.derivative);
    }

    // atan2
    template <typename Scalar, size_t N>
    Dual<Scalar, N> atan2(Dual<Scalar, N>& a, Dual<Scalar, N>& b)
    {
        // derivative is equal to that of atan(a/b), so substitute a/b into atan and simplify
        Scalar atan_d = Scalar(1) / (a.value * a.value + b.value * b.value);
        return Dual<Scalar, N>(std::atan2(a.value, b.value), (a.derivative * b.value - a.value * b.derivative) * atan_d);
    }

    template<typename Scaler, size_t N>
    FastMath::Vector<Dual<Scaler, N>, N> Vector2Dual(const FastMath::Vector<Scaler, N>& vec)
    {
        FastMath::Vector<Dual<Scaler, N>, N> out;

        for (size_t k = 0; k < N; ++k)
        {
            out(k).value = vec(k);
            out(k).derivative(k) = 1;
        }
        return out;
    }

    // retrieve the derivative elements of a vector of Duals into a matrix
    template <typename Scalar, size_t M, size_t N>
    Matrix<Scalar, M, N> collectDerivatives(Matrix<Dual<Scalar, N>, M, 1>& input)
    {
        Matrix<Scalar, M, N> jac;
        for (size_t i = 0; i < M; i++) {
            jac.row(i) = input(i, 0).derivative;
        }
        return jac;
    }



    // retrieve the real (non-derivative) elements of a matrix of Duals into an equally sized matrix
    template <typename Scalar, size_t M, size_t N, size_t D>
    Matrix<Scalar, M, N> collectReals(Matrix<Dual<Scalar, D>, M, N>& input)
    {
        Matrix<Scalar, M, N> r;
        for (size_t i = 0; i < M; i++) {
            for (size_t j = 0; j < N; j++) {
                r(i,j) = input(i,j).value;
            }
        }
        return r;
    }

    template<typename Type, size_t N>
    std::ostream& operator<<(std::ostream& os,
                             const FastMath::Dual<Type, N>& dual)
    {
        os << "[";
        os << std::setw(10) << dual.value << ";";
        for (size_t j = 0; j < N; ++j) {
            os << "\t";
            os << std::setw(10) << static_cast<double>(dual.derivative(j));
        }
        os << "]";
        return os;
    }









}

#endif //FASTMATH2_DUAL_HPP
