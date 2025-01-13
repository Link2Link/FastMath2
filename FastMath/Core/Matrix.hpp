/******************************************************************************
  文 件 名   : Matrix.hpp
  作    者   : Barry
  生成日期   : 2024-08-29
  最近修改   :
  功能描述   : 
  函数列表   :
  修改历史   :
  1.日    期   : 2024-08-29
    作    者   : Barry
    修改内容   : 创建文件
******************************************************************************/

#ifndef FASTMATH2_MATRIX_HPP
#define FASTMATH2_MATRIX_HPP

#include <cmath>
#include <limits>
#include <array>
#include <cassert>
#include <cstring>
#include <type_traits>
#include <iomanip>
#include <chrono>
#include <iostream>
#include <cstdio>

#include "Impl/impl.hpp"

#define MeasureTime_START auto start = std::chrono::high_resolution_clock::now();
#define MeasureTime_END     auto end = std::chrono::high_resolution_clock::now();   \
    std::chrono::duration<double> elapsed = end - start;                            \
    std::cout << "Elapsed time: " << elapsed.count() << "s\n";
#define LOOP(N) for (size_t __k = 0; __k<N; ++__k)

namespace FastMath
{
    using namespace Impl;
    template<typename Type, size_t M, size_t N>
    class Matrix;

    template<typename Type, size_t N>
    using SquareMatrix = Matrix<Type, N, N>;

    template<typename Type, size_t N>
    using Vector = Matrix<Type, N, 1>;

#define FAST_EXTRA_MAT_TYPEDEF(N) \
    using Matrix##N##f = SquareMatrix<float, N>;         \
    using Matrix##N##d = SquareMatrix<double, N>;        \
    using Vector##N##f = Vector<double, N>;           \
    using Vector##N##d = Vector<double, N>;

    FAST_EXTRA_MAT_TYPEDEF(2)
    FAST_EXTRA_MAT_TYPEDEF(3)
    FAST_EXTRA_MAT_TYPEDEF(4)
    FAST_EXTRA_MAT_TYPEDEF(5)
    FAST_EXTRA_MAT_TYPEDEF(6)
    FAST_EXTRA_MAT_TYPEDEF(7)
    FAST_EXTRA_MAT_TYPEDEF(8)
    FAST_EXTRA_MAT_TYPEDEF(9)
#undef FAST_EXTRA_MAT_TYPEDEF

    template <typename Type, size_t P, size_t Q, size_t M, size_t N>
    class Slice;

    template<size_t M, size_t N>
    struct is_matrix {
        static constexpr bool value = (M > 1) && (N > 1);
    };

    template<size_t M, size_t N>
    struct is_square {
        static constexpr bool value = (M == N);
    };

    template<size_t M, size_t N>
    struct is_row_vector {
        static constexpr bool value = (M == 1) && (N > 1);
    };

    template<size_t M, size_t N>
    struct is_col_vector {
        static constexpr bool value = (N == 1) && (M > 1);
    };

    template<size_t M, size_t N>
    struct is_vector {
        static constexpr bool value = (is_row_vector<M,N>::value) || (is_col_vector<M,N>::value);
    };

    template<size_t M, size_t N>
    struct is_scaler {
        static constexpr bool value = (N == 1) && (M == 1);
    };

    template<size_t M, size_t N>
    struct is_less {
        static constexpr bool value = M < N;
    };

    template<size_t M, size_t N>
    struct is_less_or_equal {
        static constexpr bool value = M <= N;
    };

    template<size_t M, size_t N>
    struct is_greater {
        static constexpr bool value = M > N;
    };

    template<size_t M, size_t N>
    struct is_greater_or_equal {
        static constexpr bool value = M >= N;
    };

#define CLAIM_MATRIX static_assert(is_matrix<M,N>::value, "Not Matrix")
#define CLAIM_SQUARE_MATRIX static_assert(is_square<M,N>::value, "Not Square Matrix")
#define CLAIM_ROW_VECTOR static_assert(is_row_vector<M,N>::value, "Not Row Vector")
#define CLAIM_COL_VECTOR static_assert(is_col_vector<M,N>::value, "Not Col Vector")
#define CLAIM_VECTOR static_assert(is_vector<M,N>::value, "Not Vector")


    template<typename Type, size_t M, size_t N>
    class Matrix
    {
    public:
        Type _data[M][N] {};    // 内部数组

    public:

        // Constructors
        Matrix() = default;

        /**
         * 从单个元素构建矩阵
         * @param val
         */
        Matrix(const Type val)
        {
            for (size_t i = 0; i < M; i++) {
                for (size_t j = 0; j < N; j++) {
                    _data[i][j] = val;
                }
            }
        }

        /***
         * 从一维数组构建矩阵
         * @param dat_data_
         */
        Matrix(const Type dat_data_[M*N])
        {
            for (size_t i = 0; i < M; i++) {
                for (size_t j = 0; j < N; j++) {
                    _data[i][j] = dat_data_[N*i + j];
                }
            }
        }

        /**
         * 从二维数组构建矩阵
         * @param dat_data_
         */
        Matrix(const Type dat_data_[M][N])
        {
            for (size_t i = 0; i < M; i++) {
                for (size_t j = 0; j < N; j++) {
                    _data[i][j] = dat_data_[i][j];
                }
            }
        }

        /**
         * 从array构建矩阵，支持列表初始化
         * @param dat_data_
         */
        Matrix(const std::array<Type, M*N> dat_data_)
        {
            for (size_t i = 0; i < M; i++) {
                for (size_t j = 0; j < N; j++) {
                    _data[i][j] = dat_data_[N*i + j];
                }
            }
        }

        /**
         * 统一初始化
         * @param values
         */
        Matrix(const std::initializer_list<Type> values)
        {
            memset(_data, 0, sizeof(double)*M*N);
            auto num = values.size();
            int k = 0;
            for (auto& value : values)
            {
                *(_data[0] + k) = value;
                k++;
                if (k == M*N)
                    break;
            }
        }

        /**
         * 从其他对象初始化
         * @param other
         */
        Matrix(const Matrix &other)
        {
            for (size_t i = 0; i < M; i++) {
                for (size_t j = 0; j < N; j++) {
                    _data[i][j] = other(i, j);
                }
            }
        }


        /***
         * 从切片初始化
         * @tparam P
         * @tparam Q
         * @param in_slice
         */
        template<size_t P, size_t Q>
        Matrix(const Slice<Type, M, N, P, Q>& in_slice)
        {
            Matrix<Type, M, N>& self = *this;
            for (size_t i = 0; i < M; i++) {
                for (size_t j = 0; j < N; j++) {
                    self(i, j) = in_slice(i, j);
                }
            }
        }

        /**
         * Accessors/ Assignment etc.
         */


        inline const Type &operator()(size_t i, size_t j) const
        {
            assert(i < M);
            assert(j < N);

            return _data[i][j];
        }

        inline Type &operator()(size_t i, size_t j)
        {
            assert(i < M);
            assert(j < N);

            return _data[i][j];
        }

        inline const Type &operator()(size_t i) const {
            assert(i < M*N);

            return *(_data[0] + i);
        }

        inline Type &operator()(size_t i) {
            assert(i < M*N);

            return *(_data[0] + i);
        }


        Matrix<Type, M, N> & operator=(const Matrix<Type, M, N> &other)
        {
            if (this != &other) {
                Matrix<Type, M, N> &self = *this;
                for (size_t i = 0; i < M; i++) {
                    for (size_t j = 0; j < N; j++) {
                        self(i, j) = other(i, j);
                    }
                }
            }
            return (*this);
        }


        template<size_t P, size_t Q>
        Matrix<Type, M, N> & operator=(const Slice<Type, M, N, P, Q>& other)
        {
            Matrix<Type, M, N> &self = *this;
            for (size_t i = 0; i < M; i++) {
                for (size_t j = 0; j < N; j++) {
                    self(i, j) = other(i, j);
                }
            }
            return (*this);
        }


        void copyTo(Type dst[M*N]) const
        {
            const Matrix<Type, M, N> &self = *this;
            for (size_t i = 0; i < M; i++) {
                for (size_t j = 0; j < N; j++) {
                    dst[N*i + j] = self(i, j);
                }
            }
        }

        void loadFrom(Type src[M*N])
        {
            for (size_t i = 0; i < M; i++) {
                for (size_t j = 0; j < N; j++) {
                    _data[i][j] = src[N*i + j];
                }
            }
        }

        void copyToColumnMajor(Type dst[M*N]) const
        {
            const Matrix<Type, M, N> &self = *this;

            for (size_t i = 0; i < M; i++) {
                for (size_t j = 0; j < N; j++) {
                    dst[i+(j*M)] = self(i, j);
                }
            }
        }

        void loadFromColumnMajor(Type src[M*N]) const
        {
            const Matrix<Type, M, N> &self = *this;

            for (size_t i = 0; i < M; i++) {
                for (size_t j = 0; j < N; j++) {
                    self(i, j) = src[i+(j*M)];
                }
            }
        }

        // !warning : this function does not check size inside
        inline void copyToArrayPtr(Type* dataPtr)
        {
            memcpy(dataPtr, this->_data, sizeof(Type)*M*N);
        }

        // !warning : this function does not check size inside
        inline void loadFromArrayPtr(Type* dataPtr)
        {
            memcpy(this->_data, dataPtr, sizeof(Type)*M*N);
        }


        /**
         * Matrix Operations
         */
        // this might use a lot of programming memory
        // since it instantiates a class for every
        // required mult pair, but it provides
        // compile time size_t checking
        template<size_t P>
        Matrix<Type, M, P> operator*(const Matrix<Type, N, P> &other) const
        {

            const Matrix<Type, M, N> &self = *this;
            Matrix<Type, M, P> res{};

            for (size_t i = 0; i < M; i++) {
                for (size_t k = 0; k < P; k++) {
                    for (size_t j = 0; j < N; j++) {
                        res(i, k) += self(i, j) * other(j, k);
                    }
                }
            }

            return res;
        }


        Matrix<Type, M, N> emult(const Matrix<Type, M, N> &other) const
        {
            Matrix<Type, M, N> res;
            const Matrix<Type, M, N> &self = *this;

            for (size_t i = 0; i < M; i++) {
                for (size_t j = 0; j < N; j++) {
                    res(i, j) = self(i, j)*other(i, j);
                }
            }

            return res;
        }

        Matrix<Type, M, N> edivide(const Matrix<Type, M, N> &other) const
        {
            Matrix<Type, M, N> res;
            const Matrix<Type, M, N> &self = *this;

            for (size_t i = 0; i < M; i++) {
                for (size_t j = 0; j < N; j++) {
                    res(i, j) = self(i, j)/other(i, j);
                }
            }

            return res;
        }

        Matrix<Type, M, N> operator+(const Matrix<Type, M, N> &other) const
        {
            Matrix<Type, M, N> res;
            const Matrix<Type, M, N> &self = *this;

            for (size_t i = 0; i < M; i++) {
                for (size_t j = 0; j < N; j++) {
                    res(i, j) = self(i, j) + other(i, j);
                }
            }

            return res;
        }

        Matrix<Type, M, N> operator-(const Matrix<Type, M, N> &other) const
        {
            Matrix<Type, M, N> res;
            const Matrix<Type, M, N> &self = *this;

            for (size_t i = 0; i < M; i++) {
                for (size_t j = 0; j < N; j++) {
                    res(i, j) = self(i, j) - other(i, j);
                }
            }
            return res;
        }


        // unary minus
        Matrix<Type, M, N> operator-() const
        {
            Matrix<Type, M, N> res;
            const Matrix<Type, M, N> &self = *this;

            for (size_t i = 0; i < M; i++) {
                for (size_t j = 0; j < N; j++) {
                    res(i, j) = -self(i, j);
                }
            }

            return res;
        }


        void operator+=(const Matrix<Type, M, N> &other)
        {
            Matrix<Type, M, N> &self = *this;
            for (size_t i = 0; i < M; i++) {
                for (size_t j = 0; j < N; j++) {
                    self(i, j) += other(i, j);
                }
            }
        }

        void operator-=(const Matrix<Type, M, N> &other)
        {
            Matrix<Type, M, N> &self = *this;
            for (size_t i = 0; i < M; i++) {
                for (size_t j = 0; j < N; j++) {
                    self(i, j) -= other(i, j);
                }
            }
        }

        template<size_t P>
        void operator*=(const Matrix<Type, N, P> &other)
        {
            Matrix<Type, M, N> &self = *this;
            self = self * other;
        }

        /**
         * Scalar Operations
         */

        Matrix<Type, M, N> operator*(Type scalar) const
        {
            Matrix<Type, M, N> res;
            const Matrix<Type, M, N> &self = *this;

            for (size_t i = 0; i < M; i++) {
                for (size_t j = 0; j < N; j++) {
                    res(i, j) = self(i, j) * scalar;
                }
            }

            return res;
        }

        inline Matrix<Type, M, N> operator/(Type scalar) const
        {
            return (*this)*(1/scalar);
        }


        Matrix<Type, M, N> operator+(Type scalar) const
        {
            Matrix<Type, M, N> res;
            const Matrix<Type, M, N> &self = *this;

            for (size_t i = 0; i < M; i++) {
                for (size_t j = 0; j < N; j++) {
                    res(i, j) = self(i, j) + scalar;
                }
            }

            return res;
        }

        inline Matrix<Type, M, N> operator-(Type scalar) const
        {
            return (*this) + (-1*scalar);
        }

        void operator*=(Type scalar)
        {
            Matrix<Type, M, N> &self = *this;

            for (size_t i = 0; i < M; i++) {
                for (size_t j = 0; j < N; j++) {
                    self(i, j) *= scalar;
                }
            }
        }

        void operator/=(Type scalar)
        {
            Matrix<Type, M, N> &self = *this;
            self *= (Type(1) / scalar);
        }

        inline void operator+=(Type scalar)
        {
            Matrix<Type, M, N> &self = *this;
            for (size_t i = 0; i < M; i++) {
                for (size_t j = 0; j < N; j++) {
                    self(i, j) += scalar;
                }
            }
        }

        inline void operator-=(Type scalar)
        {
            Matrix<Type, M, N> &self = *this;
            self += (-scalar);
        }

        bool operator==(const Matrix<Type, M, N> &other) const
        {
            for (size_t i = 0; i < M; i++) {
                for (size_t j = 0; j < N; j++) {
                    if (!Impl::isEqualF(_data[i][j], other(i,j), Impl::M_EPS)) {
                        return false;
                    }
                }
            }
            return true;
        }

        bool operator!=(const Matrix<Type, M, N> &other) const
        {
            const Matrix<Type, M, N> &self = *this;
            return self != other;
        }


        void write_string(char * buf, size_t n) const
        {
            buf[0] = '\0'; // make an empty string to begin with (we need the '\0' for strlen to work)
            const Matrix<Type, M, N> &self = *this;
            for (size_t i = 0; i < M; i++) {
                for (size_t j = 0; j < N; j++) {
                    snprintf(buf + strlen(buf), n - strlen(buf), "\t%8.8g", double(self(i, j))); // directly append to the string buffer
                }
                snprintf(buf + strlen(buf), n - strlen(buf), "\n");
            }
        }

        std::string format()
        {
            // element: tab, point, 8 digits, 4 scientific notation chars; row: newline; string: \0 end
            constexpr size_t n = 15*N*M + M + 1;
            char buf[n];
            write_string(buf, n);
            return buf;
        }

        void print() const
        {
            // element: tab, point, 8 digits, 4 scientific notation chars; row: newline; string: \0 end
            constexpr size_t n = 15*N*M + M + 1;
            char buf[n];
            write_string(buf, n);
            printf("%s\n", buf);
        }

        void print(char* name) const
        {
            // element: tab, point, 8 digits, 4 scientific notation chars; row: newline; string: \0 end
            constexpr size_t n = 15*N*M + M + 1;
            char buf[n];
            write_string(buf, n);
            printf("%s\n%s\n", name, buf);
        }



        Matrix<Type, N, M> transpose() const
        {
            Matrix<Type, N, M> res;
            const Matrix<Type, M, N> &self = *this;

            for (size_t i = 0; i < M; i++) {
                for (size_t j = 0; j < N; j++) {
                    res(j, i) = self(i, j);
                }
            }
            return res;
        }

        // tranpose alias
        inline Matrix<Type, N, M> T() const
        {
            return transpose();
        }

        template<size_t P, size_t Q>
        Slice<Type, P, Q, M, N> slice(size_t x0, size_t y0) const
        {
            return Slice<Type, P, Q, M, N>(x0, y0, this);
        }


        template<size_t P, size_t Q>
        Slice<Type, P, Q, M, N> slice(size_t x0, size_t y0)
        {
            return Slice<Type, P, Q, M, N>(x0, y0, this);
        }


        const Slice<Type, 1, N, M, N> row(size_t i) const
        {
            return slice<1, N>(i,0);
        }

        Slice<Type, 1, N, M, N> row(size_t i)
        {
            return slice<1, N>(i,0);
        }

        const Slice<Type, M, 1, M, N> col(size_t j) const
        {
            return slice<M, 1>(0,j);
        }

        Slice<Type, M, 1, M, N> col(size_t j)
        {
            return slice<M, 1>(0,j);
        }

        void setRow(size_t i, const Matrix<Type, N, 1> &row_in)
        {
            slice<1,N>(i,0) = row_in.transpose();
        }

        void setRow(size_t i, Type val)
        {
            slice<1,N>(i,0) = val;
        }

        void setCol(size_t j, const Matrix<Type, M, 1> &column)
        {
            slice<M,1>(0,j) = column;
        }

        void setCol(size_t j, Type val)
        {
            slice<M,1>(0,j) = val;
        }

        void setZero()
        {
            memset(_data, 0, sizeof(_data));
        }

        inline void zero()
        {
            setZero();
        }


        void setAll(Type val)
        {
            Matrix<Type, M, N> &self = *this;

            for (size_t i = 0; i < M; i++) {
                for (size_t j = 0; j < N; j++) {
                    self(i, j) = val;
                }
            }
        }

        size_t cols()
        {
            return N;
        }

        size_t rows()
        {
            return M;
        }

        inline void setOne()
        {
            setAll(1);
        }

        inline void setNaN()
        {
            setAll(NAN);
        }

        void setIdentity()
        {
            setZero();
            Matrix<Type, M, N> &self = *this;

            const size_t min_i = M > N ? N : M;
            for (size_t i = 0; i < min_i; i++) {
                self(i, i) = 1;
            }
        }

        inline void identity()
        {
            setIdentity();
        }


        inline void swapRows(size_t a, size_t b)
        {
            assert(a < M);
            assert(b < M);

            if (a == b) {
                return;
            }

            Matrix<Type, M, N> &self = *this;

            for (size_t j = 0; j < N; j++) {
                Type tmp = self(a, j);
                self(a, j) = self(b, j);
                self(b, j) = tmp;
            }
        }

        inline void swapCols(size_t a, size_t b)
        {
            assert(a < N);
            assert(b < N);

            if (a == b) {
                return;
            }

            Matrix<Type, M, N> &self = *this;

            for (size_t i = 0; i < M; i++) {
                Type tmp = self(i, a);
                self(i, a) = self(i, b);
                self(i, b) = tmp;
            }
        }


        [[nodiscard]] Matrix<Type, M, N> abs() const
        {
            Matrix<Type, M, N> r;
            for (size_t i=0; i<M; i++) {
                for (size_t j=0; j<N; j++) {
                    r(i,j) = Type(fabs((*this)(i,j)));
                }
            }
            return r;
        }


        [[nodiscard]] Type max() const
        {
            Type max_val = (*this)(0,0);
            for (size_t i=0; i<M; i++) {
                for (size_t j=0; j<N; j++) {
                    Type val = (*this)(i,j);
                    if (val > max_val) {
                        max_val = val;
                    }
                }
            }
            return max_val;
        }


        [[nodiscard]] Type min() const
        {
            Type min_val = (*this)(0,0);
            for (size_t i=0; i<M; i++) {
                for (size_t j=0; j<N; j++) {
                    Type val = (*this)(i,j);
                    if (val < min_val) {
                        min_val = val;
                    }
                }
            }
            return min_val;
        }

        [[nodiscard]] Type mean() const
        {
            Type sum = 0;
            for (size_t i=0; i<M; i++) {
                for (size_t j=0; j<N; j++) {
                    sum += (*this)(i,j);
                }
            }

            return sum / (M*N);
        }

        [[nodiscard]] bool isAllNan() const {
            const Matrix<Type, M, N> &self = *this;
            bool result = true;
            for (size_t i = 0; i < M; i++) {
                for (size_t j = 0; j < N; j++) {
                    result = result && isnan(self(i, j));
                }
            }
            return result;
        }

        [[nodiscard]] bool isAllsafe() const {
            const Matrix<Type, M, N> &self = *this;
            bool result = true;
            for (size_t i = 0; i < M; i++) {
                for (size_t j = 0; j < N; j++) {
                    result = result && Impl::is_safe(self(i, j));;
                }
            }
            return result;
        }

        [[nodiscard]] bool isAllzero() const {
            const Matrix<Type, M, N> &self = *this;
            bool result = true;
            for (size_t i = 0; i < M; i++) {
                for (size_t j = 0; j < N; j++) {
                    result = result && std::abs(self(i, j)) < 1E-12;;
                }
            }
            return result;
        }

        /* matrix operation in linear algebra */
        Vector<Type, (M > N ? N : M)> diag()
        {

            const Matrix<Type, M, N> &self = *this;

            const size_t min_i = M > N ? N : M;
            Vector<Type, min_i> res;

            for (size_t i = 0; i < M; i++) {
                res(i) = self(i, i);
            }

            return res;
        }
//        void Print_Binary(unsigned int x) {
//            if (x > 1) {
//                Print_Binary(x >> 1);
//            }
//            putchar((x & 1) ? '1' : '0');
//        }

        /* For Square Matrix */

        SquareMatrix<Type, N> power(size_t s)
        {
            CLAIM_SQUARE_MATRIX;
            // calc order t
            int t = 0;
            size_t temp = s;
            while (temp != 0)
            {
                t++;
                temp = temp / 2;
            }
            SquareMatrix<Type, N> Z = *this;

            auto beta = [s](size_t k) -> bool {
                return (s >> k) & 0x01;
            };

            size_t q = 0;
            while (beta(q) == 0)
            {
                Z = Z * Z;
                q++;
            }
            auto F = Z;
            for (size_t k = q + 1; k < t; ++k)
            {
                Z = Z * Z;
                if (((s >> k) & 0x01) != 0)
                {
                    F *= Z;
                }
            }

            return F;
        }

        SquareMatrix<Type, N> sqrt()
        {
            CLAIM_SQUARE_MATRIX;
            SquareMatrix<Type, N> X;
            SquareMatrix<Type, N> Y;

            SquareMatrix<Type, N> X_new = *this;
            SquareMatrix<Type, N> Y_new = SquareMatrix<Type, N>::Identity();

            int it = 0;
            while ( (X-X_new).norm_inf() > 1E-14 && it < 20)
            {
                X = X_new;
                Y = Y_new;
                X_new = (X + Y.inv())/2;
                Y_new = (Y + X.inv())/2;
                it++;
            }
            return X;
        }

        // get matrix upper right triangle in a row-major vector format
        Vector<Type, M * (M + 1) / 2> upper_right_triangle() const
        {
            CLAIM_SQUARE_MATRIX;

            Vector<Type, M * (M + 1) / 2> res;
            const SquareMatrix<Type, M> &self = *this;

            unsigned idx = 0;
            for (size_t x = 0; x < M; x++) {
                for (size_t y = x; y < M; y++) {
                    res(idx) = self(x, y);
                    ++idx;
                }
            }

            return res;
        }


        Type trace() const
        {
            CLAIM_SQUARE_MATRIX;

            Type res = 0;
            const SquareMatrix<Type, M> &self = *this;

            for (size_t i = 0; i < M; i++) {
                res += self(i, i);
            }
            return res;
        }


        Type det() {
            CLAIM_SQUARE_MATRIX;

            if constexpr (M == 2)
            {
                return (_data[0][0]*_data[1][1] - _data[0][1]*_data[1][0]);
            }
            if constexpr (M == 3)
            {
                return (
                        _data[0][0]*_data[1][1]*_data[2][2] -
                        _data[0][0]*_data[1][2]*_data[2][1] -
                        _data[0][1]*_data[1][0]*_data[2][2] +
                        _data[0][1]*_data[1][2]*_data[2][0] +
                        _data[0][2]*_data[1][0]*_data[2][1] -
                        _data[0][2]*_data[1][1]*_data[2][0]
                );
            }

            double A[M][N];
            Impl::MatCopy<M, N, Type>(MatRef(_data), MatRef(A));
            return Impl::sdet(MatRef(A), M);
        }

        inline double norm_inf()
        {
            const Matrix<Type, M, N> &mat = *this;
            double norm = 0.0;
            double sum = 0.0;
            for (size_t i = 0; i < M; ++i)
            {
                sum = 0.0;
                for (size_t j = 0; j < N; ++j)
                {
                    sum += std::abs(mat(i,j));
                }
                norm = (norm > sum) ? norm : sum;
            }
            return norm;
        }

        inline double norm_one()
        {
            const Matrix<Type, M, N> &mat = *this;
            double norm = 0.0;
            double sum = 0.0;
            for (size_t j = 0; j < N; ++j)
            {
                sum = 0.0;
                for (size_t i = 0; i < M; ++i)
                {
                    sum += std::abs(mat(i,j));
                }
                norm = (norm > sum) ? norm : sum;
            }
            return norm;
        }

        // matrix L2 norm
        Type norm() {
            const Matrix<Type, M, N>& mat = *this;
            if constexpr (is_scaler<M, N>::value)
            {
                return std::abs(mat(0));
            }
            else if constexpr (is_vector<M,N>::value)
            {
                Type square_sum = 0;
                for (int k = 0; k < M; ++k)
                {
                    square_sum += mat(k) * mat(k);
                }

                return std::sqrt(square_sum);
            }
            // 使用幂迭代方法求最大特征值
            Matrix<Type, N, N> A = mat.T() * mat;
            if (A.abs().max() < FastMath::M_EPS)
                return 0.0;
            Vector<Type, N> v(1);
            Type lambda = Impl::eig_top1<N>(MatRef(A._data), MatRef(v._data));
            return std::sqrt(lambda);
        }

        inline void normalize()
        {
            Matrix<Type, M, N>& Mine = *this;
            Type length = Mine.norm();
            if (length > Impl::M_EPS)
                Mine /= length;
        }

        int rank()
        {
            Matrix<Type, M, N> A = *this;
            return Impl::rank(MatRef(A._data), M, N);
        }

        /* For Vector */
        [[nodiscard]] Vector3d cross(const Vector3d& b) const {
            CLAIM_COL_VECTOR;
            static_assert(M==3, "only support cross product for dim 3");
            const Vector3d& a = *this;
            return Vector3d({a(1)*b(2)-a(2)*b(1), a(2)*b(0)-a(0)*b(2), a(0)*b(1)-a(1)*b(0)});
        }


        [[nodiscard]] Type dot(const Matrix<Type, M, N>& b) const {
            CLAIM_VECTOR;
            const Matrix<Type, M, N>& a = *this;
            Type sum = 0;
            for (size_t k = 0; k < M*N; ++k)
            {
                sum += a(k) * b(k);
            }
            return sum;
        }



        /* For Scalar */
        Type Scalar()
        {
            return (*this)(0);
        }

        SquareMatrix<Type, M> inv()
        {
            CLAIM_SQUARE_MATRIX;
            if constexpr (M == 2)
            {
                Matrix<Type, M, M> & A = *this;
                Matrix<Type, M, M> inv;
                Type det = A(0, 0) * A(1, 1) - A(1, 0) * A(0, 1);

                if(std::fabs(static_cast<float>(det)) < FLT_EPSILON || !Impl::is_finite(det)) {
                    inv.setZero();
                    return inv;
                }

                inv(0, 0) = A(1, 1);
                inv(1, 0) = -A(1, 0);
                inv(0, 1) = -A(0, 1);
                inv(1, 1) = A(0, 0);
                inv /= det;
                return inv;
            }

            if constexpr (M == 3)
            {
                Matrix<Type, M, M> & A = *this;
                Matrix<Type, M, M> inv;
                Type det = A(0, 0) * (A(1, 1) * A(2, 2) - A(2, 1) * A(1, 2)) -
                           A(0, 1) * (A(1, 0) * A(2, 2) - A(1, 2) * A(2, 0)) +
                           A(0, 2) * (A(1, 0) * A(2, 1) - A(1, 1) * A(2, 0));

                if(std::fabs(static_cast<float>(det)) < FLT_EPSILON || !Impl::is_finite(det)) {
                    inv.setZero();
                    return inv;
                }

                inv(0, 0) = A(1, 1) * A(2, 2) - A(2, 1) * A(1, 2);
                inv(0, 1) = A(0, 2) * A(2, 1) - A(0, 1) * A(2, 2);
                inv(0, 2) = A(0, 1) * A(1, 2) - A(0, 2) * A(1, 1);
                inv(1, 0) = A(1, 2) * A(2, 0) - A(1, 0) * A(2, 2);
                inv(1, 1) = A(0, 0) * A(2, 2) - A(0, 2) * A(2, 0);
                inv(1, 2) = A(1, 0) * A(0, 2) - A(0, 0) * A(1, 2);
                inv(2, 0) = A(1, 0) * A(2, 1) - A(2, 0) * A(1, 1);
                inv(2, 1) = A(2, 0) * A(0, 1) - A(0, 0) * A(2, 1);
                inv(2, 2) = A(0, 0) * A(1, 1) - A(1, 0) * A(0, 1);
                inv /= det;
                return inv;
            }

            SquareMatrix<Type, M> A = *this;
            auto flag = Impl::inv<M>(MatRef(A._data));
            if (flag == 0)
                A.setZero();
            return A;
        }

        bool inv(SquareMatrix<Type, M> & inv, size_t rank = M)
        {
            CLAIM_SQUARE_MATRIX;
            SquareMatrix<Type, M>& A = *this;
            SquareMatrix<Type, M> L;
            L.setIdentity();
            SquareMatrix<Type, M> U = A;
            SquareMatrix<Type, M> P;
            P.setIdentity();

            //printf("A:\n"); A.print();

            // for all diagonal elements
            for (size_t n = 0; n < rank; n++) {

                // if diagonal is zero, swap with row below
                if (fabs(U(n, n)) < Type(FLT_EPSILON)) {
                    //printf("trying pivot for row %d\n",n);
                    for (size_t i = n + 1; i < rank; i++) {

                        //printf("\ttrying row %d\n",i);
                        if (fabs(U(i, n)) > Type(FLT_EPSILON)) {
                            //printf("swapped %d\n",i);
                            U.swapRows(i, n);
                            P.swapRows(i, n);
                            L.swapRows(i, n);
                            L.swapCols(i, n);
                            break;
                        }
                    }
                }

                // failsafe, return zero matrix
                if (std::fabs(static_cast<float>(U(n, n))) < FLT_EPSILON) {
                    return false;
                }

                // for all rows below diagonal
                for (size_t i = (n + 1); i < rank; i++) {
                    L(i, n) = U(i, n) / U(n, n);

                    // add i-th row and n-th row
                    // multiplied by: -a(i,n)/a(n,n)
                    for (size_t k = n; k < rank; k++) {
                        U(i, k) -= L(i, n) * U(n, k);
                    }
                }
            }

            //printf("L:\n"); L.print();
            //printf("U:\n"); U.print();

            // solve LY=P*I for Y by forward subst
            //SquareMatrix<Type, M> Y = P;

            // for all columns of Y
            for (size_t c = 0; c < rank; c++) {
                // for all rows of L
                for (size_t i = 0; i < rank; i++) {
                    // for all columns of L
                    for (size_t j = 0; j < i; j++) {
                        // for all existing y
                        // subtract the component they
                        // contribute to the solution
                        P(i, c) -= L(i, j) * P(j, c);
                    }

                    // divide by the factor
                    // on current
                    // term to be solved
                    // Y(i,c) /= L(i,i);
                    // but L(i,i) = 1.0
                }
            }

            //printf("Y:\n"); Y.print();

            // solve Ux=y for x by back subst
            //SquareMatrix<Type, M> X = Y;

            // for all columns of X
            for (size_t c = 0; c < rank; c++) {
                // for all rows of U
                for (size_t k = 0; k < rank; k++) {
                    // have to go in reverse order
                    size_t i = rank - 1 - k;

                    // for all columns of U
                    for (size_t j = i + 1; j < rank; j++) {
                        // for all existing x
                        // subtract the component they
                        // contribute to the solution
                        P(i, c) -= U(i, j) * P(j, c);
                    }

                    // divide by the factor
                    // on current
                    // term to be solved
                    //
                    // we know that U(i, i) != 0 from above
                    P(i, c) /= U(i, i);
                }
            }

            //check sanity of results
            for (size_t i = 0; i < rank; i++) {
                for (size_t j = 0; j < rank; j++) {
                    if (!Impl::is_finite(P(i,j))) {
                        return false;
                    }
                }
            }
            //printf("X:\n"); X.print();
            inv = P;
            return true;
        }


        /**
         * cholesky decomposition
         *
         * Note: A must be positive definite
         */

        SquareMatrix<Type, M> cholesky()
        {
            CLAIM_SQUARE_MATRIX;
            SquareMatrix<Type, M> A = *this;
            int flag = Impl::chol(MatRef(A._data), M);
            if (flag == 0)
                A.setZero();
            return A;
        }


        /**
         * Full rank Cholesky factorization of A
         */
        SquareMatrix<Type, N> fullRankCholesky(size_t* rank = nullptr)
        {
            CLAIM_SQUARE_MATRIX;
            Matrix<Type, N, N> A = *this;
            // Loses one ulp accuracy per row of diag, relative to largest magnitude
            const Type tol = N * FLT_EPSILON * A.diag().max();

            Matrix<Type, N, N> L;

            size_t r = 0;
            for (size_t k = 0; k < N; k++) {

                if (r == 0) {
                    for (size_t i = k; i < N; i++) {
                        L(i, r) = A(i, k);
                    }

                } else {
                    for (size_t i = k; i < N; i++) {
                        // Compute LL = L[k:n, :r] * L[k, :r].T
                        Type LL = Type();
                        for (size_t j = 0; j < r; j++) {
                            LL += L(i, j) * L(k, j);
                        }
                        L(i, r) = A(i, k) - LL;
                    }
                }
                if (L(k, r) > tol) {
                    L(k, r) = std::sqrt(L(k, r));

                    if (k < N - 1) {
                        for (size_t i = k + 1; i < N; i++) {
                            L(i, r) = L(i, r) / L(k, r);
                        }
                    }

                    r = r + 1;
                }
            }

            // Return rank
            if (rank != nullptr)
                *rank = r;

            return L;
        }



        int LU(Matrix<Type, N, N>& L, Matrix<Type, N, N>& U)
        {
            CLAIM_SQUARE_MATRIX;
            SquareMatrix<Type, M> A = *this;
            int flag = Impl::lluu(MatRef(A._data), N, MatRef(L._data), MatRef(U._data));
            return flag;
        }

        int QR(Matrix<Type, M, M>& Q, Matrix<Type, M, N>& R)
        {
            static_assert(is_greater_or_equal<M,N>::value);
            R = *this;
            int flag = Impl::maqr(MatRef(R._data), M, N, MatRef(Q._data));
            return flag;
        }


        /**
         * 奇异值分解
         * @tparam orth 是否需要U矩阵正交化，默认为true，大于最后M-N列会进行单位正交化
         * @param U
         * @param S
         * @param V
         * @param eps
         * @return 成为为1 失败为0
         */

        template<bool orth = true>
        int SVD(Matrix<Type, M, M>& U, Matrix<Type, M, N>& S, Matrix<Type, N, N>& V, double eps=1E-10)
        {
            S = *this;
            int flag = Impl::muav<M, N>(MatRef(S._data), MatRef(U._data), MatRef(V._data), eps);
            V = V.transpose();

            if constexpr (!orth)
            {
                return flag;
            }
            // 单位正交化U剩余部分
            Vector<Type, M> vec;
            Vector<Type, M> temp;
            for (size_t k = N; k < M; ++k)
            {
                vec.setZero();
                vec(k) = 1;
                for (size_t i = 0; i < k; ++i)
                {
                    temp = U.col(i).mat();
                    vec -= temp * vec.dot(temp);
                }
                vec.normalize();
                U.col(k) = vec;
            }

            return flag;

        }


        /**
        * Geninv
        * Fast pseudoinverse based on full rank cholesky factorisation
        *
        * Courrieu, P. (2008). Fast Computation of Moore-Penrose Inverse Matrices, 8(2), 25–29. http://arxiv.org/abs/0804.4809
        */
        bool Geninv(Matrix<Type, N, M>& res)
        {
            const Matrix<Type, M, N> & G = *this;
            size_t rank;
            if (M <= N) {
                SquareMatrix<Type, M> A = G * G.transpose();
                SquareMatrix<Type, M> L = A.fullRankCholesky(&rank);

                A = L.transpose() * L;
                SquareMatrix<Type, M> X;
                if (!A.inv(X, rank)) {
                    res = Matrix<Type, N, M>();
                    return false; // LCOV_EXCL_LINE -- this can only be hit from numerical issues
                }
                // doing an intermediate assignment reduces stack usage
                A = X * X * L.transpose();
                res = G.transpose() * (L * A);

            } else {
                SquareMatrix<Type, N> A = G.transpose() * G;
                SquareMatrix<Type, N> L = A.fullRankCholesky(&rank);

                A = L.transpose() * L;
                SquareMatrix<Type, N> X;
                if(!A.inv(X, rank)) {
                    res = Matrix<Type, N, M>();
                    return false; // LCOV_EXCL_LINE -- this can only be hit from numerical issues
                }
                // doing an intermediate assignment reduces stack usage
                A = X * X * L.transpose();
                res = (L * A) * G.transpose();
            }
            return true;
        }

        // SVD求伪逆
        inline Matrix<Type, N, M> Geninv()
        {
            Matrix<Type, N, M> res;
            this->Geninv(res);
            return res;
        }

        // SVD求伪逆
        inline bool pinv(Matrix<Type, N, M>& res)
        {
            Matrix<Type, M, N> A = *this;
            Matrix<Type, M, M> U;
            Matrix<Type, N, N> VT;
            int flag = Impl::ginv<M, N>(FastMatRef(A), FastMatRef(res), FastMatRef(U), FastMatRef(VT), 1E-10);
            return flag;
        }

        // SVD求伪逆
        inline Matrix<Type, N, M> pinv()
        {
            Matrix<Type, N, M> res;
            this->pinv(res);
            return res;
        }




        SquareMatrix<Type, M> expm()
        {
            CLAIM_SQUARE_MATRIX;
//            const Matrix<Type, M, M> &A = *this;
            // ref: https://people.sc.fsu.edu/~jburkardt/cpp_src/matrix_exponential/matrix_exponential.html
            Matrix<Type, M, M> A2 = *this;
            auto a_norm = A2.norm_inf();
            int ee = ( int ) ( log2( a_norm ) ) + 1;
            int s = (0 > ee + 1) ? 0 : (ee + 1);
            double t = 1.0 / pow(2, s);    // t = 2^(-s)
            A2 = A2*t;
            auto x = A2;
            auto c = 0.5;
            auto e = SquareMatrix<Type, M>::Identity();
            e = e + A2*c;
            auto d = SquareMatrix<Type, M>::Identity();
            d = d - A2*c;
            int p = 1;
            int q = 6;
            for (int k = 2; k <= q; ++k)
            {
                c = c * ( double ) ( q - k + 1 ) / ( double ) ( k * ( 2 * q - k + 1 ) );
                x = A2*x;
                e = x*c + e;
                if (p)
                    d = x*c + d;
                else
                    d = -x*c + d;
                p = !p;
            }
            e = d.inv() * e;
            for (int k = 1; k <= s; ++k)
            {
                e = e*e;
            }

            return e;
        }

        // 使用反向缩放平方法计算logm
        SquareMatrix<Type, M> logm()
        {
            CLAIM_SQUARE_MATRIX;
            SquareMatrix<Type, M> A = *this;
            SquareMatrix<Type, M> I;
            I.setIdentity();
            // 平方收缩
            int k = 0;
            while ((A-I).norm_inf() > 1e-2)
            {
                k++;
                A = A.sqrt();
            }
            // Pade近似
            SquareMatrix<Type, M> A_I = A - I;
            SquareMatrix<Type, M> A_I2 = A_I * A_I;
            SquareMatrix<Type, M> A_I3 = A_I2 * A_I;
            SquareMatrix<Type, M> DA = 60*I + 90*A_I + 36*A_I2 + 3*A_I3;
            SquareMatrix<Type, M> NA = 60*A_I + 60*A_I2 + 11*A_I3;

            A = DA.inv() * NA * std::pow(2.0, k);

            return A;
        }



        // make block diagonal symmetric by taking the average of the two corresponding off diagonal values
        template <size_t Width>
        void makeBlockSymmetric(size_t first)
        {
            CLAIM_SQUARE_MATRIX;
            static_assert(Width <= M, "Width bigger than matrix");
            assert(first + Width <= M);

            SquareMatrix<Type, M> &self = *this;
            if(Width>1) {
                for (size_t row_idx = first+1; row_idx < first+Width; row_idx++) {
                    for (size_t col_idx = first; col_idx < row_idx; col_idx++) {
                        Type tmp = (self(row_idx,col_idx) + self(col_idx,row_idx)) / Type(2);
                        self(row_idx,col_idx) = tmp;
                        self(col_idx,row_idx) = tmp;
                    }
                }
            }
        }

        // make rows and columns symmetric by taking the average of the two corresponding off diagonal values
        template <size_t Width>
        void makeRowColSymmetric(size_t first)
        {
            CLAIM_SQUARE_MATRIX;
            static_assert(Width <= M, "Width bigger than matrix");
            assert(first + Width <= M);

            SquareMatrix<Type, M> &self = *this;
            self.template makeBlockSymmetric<Width>(first);
            for (size_t row_idx = first; row_idx < first+Width; row_idx++) {
                for (size_t col_idx = 0; col_idx < first; col_idx++) {
                    Type tmp = (self(row_idx,col_idx) + self(col_idx,row_idx)) / Type(2);
                    self(row_idx,col_idx) = tmp;
                    self(col_idx,row_idx) = tmp;
                }
                for (size_t col_idx = first+Width; col_idx < M; col_idx++) {
                    Type tmp = (self(row_idx,col_idx) + self(col_idx,row_idx)) / Type(2);
                    self(row_idx,col_idx) = tmp;
                    self(col_idx,row_idx) = tmp;
                }
            }
        }

        // checks if block diagonal is symmetric
        template <size_t Width>
        bool isBlockSymmetric(size_t first, const Type eps = Type(1e-8f))
        {
            CLAIM_SQUARE_MATRIX;
            static_assert(Width <= M, "Width bigger than matrix");
            assert(first + Width <= M);

            SquareMatrix<Type, M> &self = *this;
            if(Width>1) {
                for (size_t row_idx = first+1; row_idx < first+Width; row_idx++) {
                    for (size_t col_idx = first; col_idx < row_idx; col_idx++) {
                        if(!Impl::isEqualF(self(row_idx,col_idx), self(col_idx,row_idx), eps)) {
                            return false;
                        }
                    }
                }
            }
            return true;
        }


        // checks if rows and columns are symmetric
        template <size_t Width>
        bool isRowColSymmetric(size_t first, const Type eps = Type(1e-8f))
        {
            CLAIM_SQUARE_MATRIX;
            static_assert(Width <= M, "Width bigger than matrix");
            assert(first + Width <= M);

            SquareMatrix<Type, M> &self = *this;
            for (size_t row_idx = first; row_idx < first+Width; row_idx++) {
                for (size_t col_idx = 0; col_idx < first; col_idx++) {
                    if(!Impl::isEqualF(self(row_idx,col_idx), self(col_idx,row_idx), eps)) {
                        return false;
                    }
                }
                for (size_t col_idx = first+Width; col_idx < M; col_idx++) {
                    if(!isEqualF(self(row_idx,col_idx), self(col_idx,row_idx), eps)) {
                        return false;
                    }
                }
            }
            return self.template isBlockSymmetric<Width>(first, eps);
        }



        // 特殊矩阵生成
        static Matrix<Type, M, N> Identity()
        {
            Matrix<Type, M, N> A;
            A.setIdentity();
            return A;
        }

        static Matrix<Type, M, N> Zeros()
        {
            Matrix<Type, M, N> A;
            A.setZero();
            return A;
        }

        static Matrix<Type, M, N> Ones()
        {
            Matrix<Type, M, N> A;
            A.setAll(1);
            return A;
        }


        static Matrix<Type, M, N> Nans() {
            Matrix<Type, M, N> m;
            m.setNaN();
            return m;
        }

        static Matrix<Type, M, N> Rand() {
            Matrix<Type, M, N> m;
            for (int i = 0; i < M; ++i)
                for (int j = 0; j < N; ++j)
                    m(i, j) = Impl::rand();
            return m;
        }



    };


    template <typename Type, size_t P, size_t Q, size_t M, size_t N>
    class Slice {
    public:
        friend Matrix<Type, P, Q>;

        Slice(size_t x0, size_t y0, const Matrix<Type, M, N>* data) :
                _x0(x0),
                _y0(y0),
                _data(const_cast<Matrix<Type, M, N>*>(data)) {
            static_assert(P <= M, "Slice rows bigger than backing matrix");
            static_assert(Q <= N, "Slice cols bigger than backing matrix");
            assert(x0 + P <= M);
            assert(y0 + Q <= N);
        }

        const Type &operator()(size_t i, size_t j) const
        {
            assert(i < P);
            assert(j < Q);

            return (*_data)(_x0 + i, _y0 + j);
        }

        Type &operator()(size_t i, size_t j)

        {
            assert(i < P);
            assert(j < Q);

            return (*_data)(_x0 + i, _y0 + j);
        }

        template<size_t MM, size_t NN>
        Slice<Type, P, Q, M, N>& operator=(const Slice<Type, P, Q, MM, NN>& other)
        {
            Slice<Type, P, Q, M, N>& self = *this;
            for (size_t i = 0; i < P; i++) {
                for (size_t j = 0; j < Q; j++) {
                    self(i, j) = other(i, j);
                }
            }
            return self;
        }

        Slice<Type, P, Q, M, N>& operator=(const Matrix<Type, P, Q>& other)
        {
            Slice<Type, P, Q, M, N>& self = *this;
            for (size_t i = 0; i < P; i++) {
                for (size_t j = 0; j < Q; j++) {
                    self(i, j) = other(i, j);
                }
            }
            return self;
        }

        Slice<Type, P, Q, M, N>& operator=(const Type& other)
        {
            Slice<Type, P, Q, M, N>& self = *this;
            for (size_t i = 0; i < P; i++) {
                for (size_t j = 0; j < Q; j++) {
                    self(i, j) = other;
                }
            }
            return self;
        }

        // allow assigning vectors to a slice that are in the axis
        template <size_t DUMMY = 1> // make this a template function since it only exists for some instantiations
        Slice<Type, 1, Q, M, N>& operator=(const Vector<Type, Q>& other)
        {
            Slice<Type, 1, Q, M, N>& self = *this;
            for (size_t j = 0; j < Q; j++) {
                self(0, j) = other(j);
            }
            return self;
        }

        template<size_t MM, size_t NN>
        Slice<Type, P, Q, M, N>& operator+=(const Slice<Type, P, Q, MM, NN>& other)
        {
            Slice<Type, P, Q, M, N>& self = *this;
            for (size_t i = 0; i < P; i++) {
                for (size_t j = 0; j < Q; j++) {
                    self(i, j) += other(i, j);
                }
            }
            return self;
        }

        Slice<Type, P, Q, M, N>& operator+=(const Matrix<Type, P, Q>& other)
        {
            Slice<Type, P, Q, M, N>& self = *this;
            for (size_t i = 0; i < P; i++) {
                for (size_t j = 0; j < Q; j++) {
                    self(i, j) += other(i, j);
                }
            }
            return self;
        }

        Slice<Type, P, Q, M, N>& operator+=(const Type& other)
        {
            Slice<Type, P, Q, M, N>& self = *this;
            for (size_t i = 0; i < P; i++) {
                for (size_t j = 0; j < Q; j++) {
                    self(i, j) += other;
                }
            }
            return self;
        }

        template<size_t MM, size_t NN>
        Slice<Type, P, Q, M, N>& operator-=(const Slice<Type, P, Q, MM, NN>& other)
        {
            Slice<Type, P, Q, M, N>& self = *this;
            for (size_t i = 0; i < P; i++) {
                for (size_t j = 0; j < Q; j++) {
                    self(i, j) -= other(i, j);
                }
            }
            return self;
        }

        Slice<Type, P, Q, M, N>& operator-=(const Matrix<Type, P, Q>& other)
        {
            Slice<Type, P, Q, M, N>& self = *this;
            for (size_t i = 0; i < P; i++) {
                for (size_t j = 0; j < Q; j++) {
                    self(i, j) -= other(i, j);
                }
            }
            return self;
        }

        Slice<Type, P, Q, M, N>& operator-=(const Type& other)
        {
            Slice<Type, P, Q, M, N>& self = *this;
            for (size_t i = 0; i < P; i++) {
                for (size_t j = 0; j < Q; j++) {
                    self(i, j) -= other;
                }
            }
            return self;
        }

        Slice<Type, P, Q, M, N>& operator*=(const Type& other)
        {
            Slice<Type, P, Q, M, N>& self = *this;
            for (size_t i = 0; i < P; i++) {
                for (size_t j = 0; j < Q; j++) {
                    self(i, j) *= other;
                }
            }
            return self;
        }

        Slice<Type, P, Q, M, N>& operator/=(const Type& other)
        {
            return operator*=(Type(1) / other);
        }

        Matrix<Type, P, Q> operator*(const Type& other) const
        {
            const Slice<Type, P, Q, M, N>& self = *this;
            Matrix<Type, P, Q> res;
            for (size_t i = 0; i < P; i++) {
                for (size_t j = 0; j < Q; j++) {
                    res(i, j) = self(i, j) * other;
                }
            }
            return res;
        }

        Matrix<Type, P, Q> operator/(const Type& other) const
        {
            const Slice<Type, P, Q, M, N>& self = *this;
            return self * (Type(1) / other);
        }

        template<size_t R, size_t S>
        const Slice<Type, R, S, M, N> slice(size_t x0, size_t y0) const
        {
            return Slice<Type, R, S, M, N>(x0 + _x0, y0 + _y0, _data);
        }

        template<size_t R, size_t S>
        Slice<Type, R, S, M, N> slice(size_t x0, size_t y0)
        {
            return Slice<Type, R, S, M, N>(x0 + _x0, y0 + _y0, _data);
        }

        void setVal(Type val)
        {
            Slice<Type, P, Q, M, N> &self = *this;
            for (size_t i = 0; i < P; i++) {
                for (size_t j = 0; j < Q; j++) {
                    self(i, j) = val;
                }
            }
        }

        void setZeros()
        {
            this->setVal(0);
        }

        void copyTo(Type dst[P*Q]) const
        {
            const Slice<Type, P, Q, M, N> &self = *this;

            for (size_t i = 0; i < P; i++) {
                for (size_t j = 0; j < Q; j++) {
                    dst[i*N+j] = self(i, j);
                }
            }
        }

        void copyToColumnMajor(Type dst[P*Q]) const
        {
            const Slice<Type, P, Q, M, N> &self = *this;

            for (size_t i = 0; i < P; i++) {
                for (size_t j = 0; j < Q; j++) {
                    dst[i+(j*M)] = self(i, j);
                }
            }
        }

        Vector<Type, P<Q?P:Q> diag()
        {
            const Slice<Type, P, Q, M, N>& self = *this;
            Vector<Type,P<Q?P:Q> res;
            for (size_t j = 0; j < (P<Q?P:Q); j++) {
                res(j) = self(j,j);
            }
            return res;
        }

        Type max() const
        {
            Type max_val = (*this)(0,0);

            for (size_t i = 0; i < P; i++) {
                for (size_t j = 0; j < Q; j++) {
                    Type val = (*this)(i,j);

                    if (val > max_val) {
                        max_val = val;
                    }
                }
            }

            return max_val;
        }

        Type min() const
        {
            Type min_val = (*this)(0,0);

            for (size_t i = 0; i < P; i++) {
                for (size_t j = 0; j < Q; j++) {
                    Type val = (*this)(i,j);

                    if (val < min_val) {
                        min_val = val;
                    }
                }
            }

            return min_val;
        }

        std::string format()
        {
            Matrix<Type, P, Q> mat(*this);
            return mat.format();
        }

        Matrix<Type, P, Q> mat()
        {
            return static_cast<Matrix<Type, P, Q>>(*this);
        }

    private:
        size_t _x0, _y0;
        Matrix<Type,M,N>* _data;
    };



    // extra functions
    template<typename Type, size_t  M, size_t N>
    inline std::ostream& operator<<(std::ostream& os,
                             const Matrix<Type, M, N>& matrix)
    {
        os << std::fixed << std::setprecision(4); // 保留4位小数

        // 计算每个元素的最大宽度
        size_t maxWidth = 0;
        for (size_t i = 0; i < M; ++i) {
            for (size_t j = 0; j < N; ++j) {
                std::ostringstream oss;
                oss << std::fixed << std::setprecision(4) << matrix(i, j);
                maxWidth = std::max(maxWidth, oss.str().length());
            }
        }

        // 打印矩阵
        for (size_t i = 0; i < M; ++i) {
            os << "[ ";
            for (size_t j = 0; j < N; ++j) {
                os << std::setw(maxWidth + 1) << matrix(i, j);
            }
            os << " ]" << std::endl;
        }
        return os;

    }

    template<typename Type, size_t  P, size_t Q, size_t  M, size_t N>
    inline std::ostream& operator<<(std::ostream& os,
                             Slice<Type, P, Q, M, N> slice)
    {
        os << slice.mat();
        return os;
    }


    template<typename Type, size_t M, size_t N1, size_t N2>
    inline Matrix<Type, M, N1+N2> cat_row(Matrix<Type, M, N1> mat1, Matrix<Type, M, N2> mat2)
    {
        Matrix<Type, M, N1+N2> mat;
        for (size_t i = 0; i < M; ++i)
        {
            for (size_t j = 0; j < N1; ++j) {
                mat(i,j) = mat1(i,j);
            }
        }
        for (size_t i = 0; i < M; ++i)
        {
            for (size_t j = 0; j < N2; ++j) {
                mat(i,j+N1) = mat2(i,j);
            }
        }

        return mat;

    }

    template<typename Type, size_t M1, size_t M2, size_t N>
    inline Matrix<Type, M1+M2, N> cat_col(Matrix<Type, M1, N> mat1, Matrix<Type, M2, N> mat2)
    {
        Matrix<Type, M1+M2, N> mat;
        for (size_t i = 0; i < M1; ++i)
        {
            for (size_t j = 0; j < N; ++j) {
                mat(i,j) = mat1(i,j);
            }
        }
        for (size_t i = 0; i < M2; ++i)
        {
            for (size_t j = 0; j < N; ++j) {
                mat(i+M1,j) = mat2(i,j);
            }
        }

        return mat;

    }

    template <typename Type, size_t N>
    inline SquareMatrix<Type, N> diag(Vector<Type, N> vec)
    {
        SquareMatrix<Type, N> A;
        A.setZero();
        for (size_t k = 0; k < N; ++k)
            A(k,k) = vec(k);
        return A;
    }


    // 允许隐式类型转换
    template<typename Type1, typename Type2, size_t  M, size_t N>
    inline Matrix<Type2, M, N> operator*(Type1 scalar, const Matrix<Type2, M, N> &other)
    {
        return other * scalar;
    }

    // 允许隐式类型转换
    template<typename Type1, typename Type2, size_t  M, size_t N>
    inline Matrix<Type2, M, N> operator+(Type1 scalar, const Matrix<Type2, M, N> &other)
    {
        return other + scalar;
    }

    template<typename Type1, typename Type2, size_t  M, size_t N>
    inline Matrix<Type2, M, N> operator-(Type1 scalar, const Matrix<Type2, M, N> &other)
    {
        return -other + scalar;
    }

    template<typename Type1, typename Type2, size_t  M, size_t N>
    inline bool operator<(const Matrix<Type1, M, N> &A, const Matrix<Type2, M, N> &B)
    {
        for (size_t i = 0; i < M; ++i)
            for (size_t j = 0; j < N; ++j)
                if (A(i,j) > B(i,j))
                    return false;
        return true;
    }

    template<typename Type1, typename Type2, size_t  M, size_t N>
    inline bool operator>(const Matrix<Type1, M, N> &A, const Matrix<Type2, M, N> &B)
    {
        return !(A < B);
    }


    namespace typeFunction
    {
        template<typename Type>
        inline Type min(const Type x, const Type y) {
            bool x_is_nan = std::isnan(x);
            bool y_is_nan = std::isnan(y);
            // take the non-nan value if there is one
            if (x_is_nan || y_is_nan) {
                if (x_is_nan && !y_is_nan) {
                    return y;
                }
                // either !x_is_nan && y_is_nan or both are NAN anyways
                return x;
            }
            return (x < y) ? x : y;
        }

        template<typename Type>
        inline Type max(const Type x, const Type y) {
            bool x_is_nan = std::isnan(x);
            bool y_is_nan = std::isnan(y);
            // take the non-nan value if there is one
            if (x_is_nan || y_is_nan) {
                if (x_is_nan && !y_is_nan) {
                    return y;
                }
                // either !x_is_nan && y_is_nan or both are NAN anyways
                return x;
            }
            return (x > y) ? x : y;
        }

        template<typename Type>
        inline Type constrain(const Type x, const Type lower_bound, const Type upper_bound) {
            if (lower_bound > upper_bound) {
                return NAN;
            } else if(std::isnan(x)) {
                return NAN;
            } else {
                return typeFunction::max(lower_bound, typeFunction::min(upper_bound, x));
            }
        }
    }

    template<typename Type, size_t  M, size_t N>
    inline Matrix<Type, M, N> min(const Matrix<Type, M, N> &x, const Type scalar_upper_bound) {
        Matrix<Type,M,N> m;
        for (size_t i = 0; i < M; i++) {
            for (size_t j = 0; j < N; j++) {
                m(i,j) = typeFunction::min(x(i,j),scalar_upper_bound);
            }
        }
        return m;
    }

    template<typename Type, size_t  M, size_t N>
    inline Matrix<Type, M, N> min(const Type scalar_upper_bound, const Matrix<Type, M, N> &x) {
        return min(x, scalar_upper_bound);
    }

    template<typename Type, size_t  M, size_t N>
    inline Matrix<Type, M, N> min(const Matrix<Type, M, N> &x1, const Matrix<Type, M, N> &x2) {
        Matrix<Type,M,N> m;
        for (size_t i = 0; i < M; i++) {
            for (size_t j = 0; j < N; j++) {
                m(i,j) = typeFunction::min(x1(i,j),x2(i,j));
            }
        }
        return m;
    }

    template<typename Type, size_t  M, size_t N>
    inline Matrix<Type, M, N> max(const Matrix<Type, M, N> &x, const Type scalar_lower_bound) {
        Matrix<Type,M,N> m;
        for (size_t i = 0; i < M; i++) {
            for (size_t j = 0; j < N; j++) {
                m(i,j) = typeFunction::max(x(i,j),scalar_lower_bound);
            }
        }
        return m;
    }

    template<typename Type, size_t  M, size_t N>
    inline Matrix<Type, M, N> max(const Type scalar_lower_bound, const Matrix<Type, M, N> &x) {
        return max(x, scalar_lower_bound);
    }

    template<typename Type, size_t  M, size_t N>
    inline Matrix<Type, M, N> max(const Matrix<Type, M, N> &x1, const Matrix<Type, M, N> &x2) {
        Matrix<Type,M,N> m;
        for (size_t i = 0; i < M; i++) {
            for (size_t j = 0; j < N; j++) {
                m(i,j) = typeFunction::max(x1(i,j),x2(i,j));
            }
        }
        return m;
    }

    template<typename Type, size_t  M, size_t N>
    inline Matrix<Type, M, N> constrain(const Matrix<Type, M, N> &x,
                                 const Type scalar_lower_bound,
                                 const Type scalar_upper_bound) {
        Matrix<Type,M,N> m;
        if (scalar_lower_bound > scalar_upper_bound) {
            m.setNaN();
        } else {
            for (size_t i = 0; i < M; i++) {
                for (size_t j = 0; j < N; j++) {
                    m(i,j) = typeFunction::constrain(x(i,j), scalar_lower_bound, scalar_upper_bound);
                }
            }
        }
        return m;
    }

    template<typename Type, size_t  M, size_t N>
    inline Matrix<Type, M, N> constrain(const Matrix<Type, M, N> &x,
                                 const Matrix<Type, M, N> &x_lower_bound,
                                 const Matrix<Type, M, N> &x_upper_bound) {
        Matrix<Type,M,N> m;
        for (size_t i = 0; i < M; i++) {
            for (size_t j = 0; j < N; j++) {
                m(i,j) = typeFunction::constrain(x(i,j), x_lower_bound(i,j), x_upper_bound(i,j));
            }
        }
        return m;
    }





}



#endif //FASTMATH2_MATRIX_HPP
