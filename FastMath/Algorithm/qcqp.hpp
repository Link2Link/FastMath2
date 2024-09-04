/******************************************************************************
  文 件 名   : qcqp.hpp
  作    者   : Barry
  生成日期   : 2024-09-03
  最近修改   :
  功能描述   : 
  函数列表   :
  修改历史   :
  1.日    期   : 2024-09-03
    作    者   : Barry
    修改内容   : 创建文件
******************************************************************************/

#ifndef FASTMATH2_QCQP_HPP
#define FASTMATH2_QCQP_HPP

#include "FastMath.hpp"

namespace FastMath::Algorithm
{

    namespace QCQP_details
    {

        /**
         * solve QCQP in 2 dimensions:
         *  min  0.5*x'*A*x + x'*b
         *      s.t.  sum (xi/di)^2 <= r^2
         *
         * @param res   优化结果
         * @param Ain   2维方阵A
         * @param bin   2维向量b
         * @param d     2维向量d
         * @param r     半径r
         * @param max_it 最大迭代次数
         * @return      未触发约束返回0，触发约束返回1
         */
        inline int QCQP2(double* res, const double* Ain, const double* bin,
                  const double* d, double r, int mat_it = 20)
        {
            double A11, A22, A12, b1, b2;
            double P11, P22, P12, det, detinv, v1, v2, la, val, deriv;

            // scale A,b so that constraint becomes x'*x <= r*r
            b1 = bin[0]*d[0];
            b2 = bin[1]*d[1];
            A11 = Ain[0]*d[0]*d[0];
            A22 = Ain[3]*d[1]*d[1];
            A12 = Ain[1]*d[0]*d[1];

            // Newton iteration
            la = 0;
            for (int iter=0; iter < mat_it; iter++) {
                // det(A+la)
                det = (A11+la)*(A22+la) - A12*A12;

                // check SPD, with 1e-10 threshold
                if (det < 1e-10) {
                    res[0] = 0;
                    res[1] = 0;
                    return 0;
                }

                // P = inv(A+la)
                detinv = 1/det;
                P11 = (A22+la)*detinv;
                P22 = (A11+la)*detinv;
                P12 = -A12*detinv;

                // v = -P*b
                v1 = -P11*b1 - P12*b2;
                v2 = -P12*b1 - P22*b2;

                // val = v'*v - r*r
                val = v1*v1 + v2*v2 - r*r;
                // check for convergence, or initial solution inside constraint set
                if (val < 1e-10) {
                    break;
                }

                // deriv = -2 * v' * P * v
                deriv = -2.0*(P11*v1*v1 + 2.0*P12*v1*v2 + P22*v2*v2);

                // compute update, exit if too small
                double delta = -val/deriv;
                if (delta < 1e-10) {
                    break;
                }

                // update
                la += delta;
            }

            // undo scaling
            res[0] = v1*d[0];
            res[1] = v2*d[1];

            return (la != 0);
        }

/**
         * solve QCQP in 3 dimensions:
         *  min  0.5*x'*A*x + x'*b
         *      s.t.  sum (xi/di)^2 <= r^2
         *
         * @param res   优化结果
         * @param Ain   2维方阵A
         * @param bin   2维向量b
         * @param d     2维向量d
         * @param r     半径r
         * @param max_it 最大迭代次数
         * @return      未触发约束返回0，触发约束返回1
         */
        inline int QCQP3(double* res, const double* Ain, const double* bin,
                         const double* d, double r, int mat_it = 20)
        {
            double A11, A22, A33, A12, A13, A23, b1, b2, b3;
            double P11, P22, P33, P12, P13, P23, det, detinv, v1, v2, v3, la, val, deriv;

            // scale A,b so that constraint becomes x'*x <= r*r
            b1 = bin[0]*d[0];
            b2 = bin[1]*d[1];
            b3 = bin[2]*d[2];
            A11 = Ain[0]*d[0]*d[0];
            A22 = Ain[4]*d[1]*d[1];
            A33 = Ain[8]*d[2]*d[2];
            A12 = Ain[1]*d[0]*d[1];
            A13 = Ain[2]*d[0]*d[2];
            A23 = Ain[5]*d[1]*d[2];

            // Newton iteration
            la = 0;
            for (int iter=0; iter < 20; iter++) {
                // unscaled P
                P11 = (A22+la)*(A33+la) - A23*A23;
                P22 = (A11+la)*(A33+la) - A13*A13;
                P33 = (A11+la)*(A22+la) - A12*A12;
                P12 = A13*A23 - A12*(A33+la);
                P13 = A12*A23 - A13*(A22+la);
                P23 = A12*A13 - A23*(A11+la);

                // det(A+la)
                det = (A11+la)*P11 + A12*P12 + A13*P13;

                // check SPD, with 1e-10 threshold
                if (det < 1e-10) {
                    res[0] = 0;
                    res[1] = 0;
                    res[2] = 0;
                    return 0;
                }

                // detinv
                detinv = 1/det;

                // final P
                P11 *= detinv;
                P22 *= detinv;
                P33 *= detinv;
                P12 *= detinv;
                P13 *= detinv;
                P23 *= detinv;

                // v = -P*b
                v1 = -P11*b1 - P12*b2 - P13*b3;
                v2 = -P12*b1 - P22*b2 - P23*b3;
                v3 = -P13*b1 - P23*b2 - P33*b3;

                // val = v'*v - r*r
                val = v1*v1 + v2*v2 + v3*v3 - r*r;

                // check for convergence, or initial solution inside constraint set
                if (val < 1e-10) {
                    break;
                }

                // deriv = -2 * v' * P * v
                deriv = -2.0*(P11*v1*v1 + P22*v2*v2 + P33*v3*v3)
                        -4.0*(P12*v1*v2 + P13*v1*v3 + P23*v2*v3);

                // compute update, exit if too small
                double delta = -val/deriv;
                if (delta < 1e-10) {
                    break;
                }

                // update
                la += delta;
            }

            // undo scaling
            res[0] = v1*d[0];
            res[1] = v2*d[1];
            res[2] = v3*d[2];

            return (la != 0);
        }

        // Cholesky solve
        template<size_t N>
        Vector<double, N> cholSolve(const SquareMatrix<double, N> mat, const Vector<double, N>& vec) {
            Vector<double, N> res = vec;
            double sum = 0;
            // forward substitution: solve L*res = vec
            for (int i=0; i < N; i++) {
                if (i) {
                    for (int k=0; k < i; k++)
                        res(i) -= mat(i, k) * res(k);
                }

                // diagonal
                res(i) /= mat(i, i);
            }

            // backward substitution: solve L'*res = res
            for (int i=N-1; i >= 0; i--) {
                if (i < N-1) {
                    for (int j=i+1; j < N; j++) {
                        res(i) -= mat(j*N+i) * res(j);
                    }
                }
                // diagonal
                res(i) /= mat(i*(N+1));
            }
            return res;
        }


    }




    /**
     * solve QCQP in N dimensions:
     *  min  0.5*x'*A*x + x'*b
     *      s.t.  sum (xi/di)^2 <= r^2
     * @tparam N    优化维度n
     * @param res   优化结果
     * @param Ain   n维方阵A
     * @param bin   n维向量b
     * @param d     n维向量d
     * @param r     半径r
     * @param max_it 最大迭代次数
     * @return      未触发约束返回0，触发约束返回1
     */

    template<size_t N>
    inline int QCQP(FastMath::Vector<double, N>& res,
                    const FastMath::Matrix<double, N, N>& Ain,
                    const FastMath::Vector<double, N>& bin,
                    const FastMath::Vector<double, N>& d,
                    double r, int max_it = 20)
    {
        FastMath::Matrix<double, N, N> A = Ain;
        FastMath::Matrix<double, N, N> Ala;
        FastMath::Vector<double, N> b = bin;
        double la, val, deriv;
        FastMath::Vector<double, N> tmp;


        // scale A,b so that constraint becomes x'*x <= r*r
        for (int i=0; i < N; i++) {
            b(i) = bin(i) * d(i);

            for (int j=0; j < N; j++) {
                A(j+i*N) = Ain(j+i*N) * d(i) * d(j);
            }
        }

        // Newton iteration
        la = 0;
        for (int iter=0; iter < max_it; iter++) {
            // make A+la
            Ala = A;
            for (int i=0; i < N; i++) {
                Ala(i*(N+1)) += la;
            }

            size_t rank;
            Ala = Ala.fullRankCholesky(&rank);
            if (rank < N)
            {
                res.setZero();
                return 0;
            }

            res =  - QCQP_details::cholSolve(Ala, b);
            val = res.dot(res) - r*r;

            // check for convergence, or initial solution inside constraint set
            if (val < 1e-12) {
                break;
            }

            tmp = QCQP_details::cholSolve(Ala, res);
            deriv = -2.0 * res.dot(tmp);

            // compute update, exit if too small
            double delta = -val/deriv;
            if (delta < 1e-12) {
                break;
            }

            // update
            la += delta;
        }

        res = res.emult(d);
        return (la != 0);
    }


    template<>
    inline int QCQP<2>(FastMath::Vector<double, 2>& res,
                    const FastMath::Matrix<double, 2, 2>& Ain,
                    const FastMath::Vector<double, 2>& bin,
                    const FastMath::Vector<double, 2>& d,
                    double r, int max_it)
    {
        return QCQP_details::QCQP2(res._data[0], Ain._data[0], bin._data[0],
                            d._data[0], r, max_it);
    }

    template<>
    inline int QCQP<3>(FastMath::Vector<double, 3>& res,
                       const FastMath::Matrix<double, 3, 3>& Ain,
                       const FastMath::Vector<double, 3>& bin,
                       const FastMath::Vector<double, 3>& d,
                       double r, int max_it)
    {
        return QCQP_details::QCQP3(res._data[0], Ain._data[0], bin._data[0],
                                   d._data[0], r, max_it);
    }
}

#endif //FASTMATH2_QCQP_HPP
