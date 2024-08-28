/******************************************************************************
  文 件 名   : rand.hpp
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

#ifndef FASTMATH_RAND_HPP
#define FASTMATH_RAND_HPP

#include <cmath>
#include "Impl/Common.hpp"

namespace FastMath::Impl
{
    class RND
    {
    private:
        int  R;       //随机数种子
    public:
        RND(int r=123) : R(r)
        {}

        //产生0到1之间均匀分布的一个随机数
        double rnd1()
        {
            int m;
            double s,u,v,p;
            s=65536.0; u=2053.0; v=13849.0;
            m=(int)(R/s); R=R-m*s;
            R=u*R+v;
            m=(int)(R/s); R=(int)(R-m*s);
            p=R/s;
            return(p);
        }

        //产生给定区间[a，b]内均匀分布的一个随机整数
        int rndab(int a, int b)
        {
            int k,j,m,i,p;
            k=b-a+1; j=2;
            while (j<k) j=j+j;
            m=4*j; k=R; i=1;
            while (i<=1)
            {
                k=k+k+k+k+k;
                k=k%m; j=k/4+a;
                if (j<=b) { p=j; i=i+1;}
            }
            R=k;
            return(p);
        }

        //产生给定均值u与方差g*g的正态分布的一个随机数
        double rndg(double u = 0, double g = 1)
        {
            int i,m;
            double s,w,v,t;
            s=65536.0; w=2053.0; v=13849.0;
            t=0.0;
            for (i=1; i<=12; i++)
            {
                R=R*w+v; m=(int)(R/s);
                R=R-m*s; t=t+R/s;
            }
            t=u+g*(t-6.0);
            return(t);
        }

    };

    static RND _rng(123);   // 全局随机种子

    inline double rand(double a = 0, double b = 1)
    {
        return _rng.rnd1() * (b - a) + a;
    }

    inline double randn(double avg, double std_)
    {
        return _rng.rndg(avg, std_);
    }

    inline int randi(int a, int b)
    {
        return _rng.rndab(a, b);
    }



}


#endif //FASTMATH_RAND_HPP
