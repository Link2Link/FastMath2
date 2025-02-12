
# 1. FastMath

FastMath是一个纯hpp文件实现的数学计算库，其设计上尽量接近Eigen的风格，实现对数学运算的高层次支持；同时又力求做到效率兼顾。FastMath不依赖任何第三方库，内部算法实现选取尽数值稳定性较好同时又简洁易懂易于维护的算法。

## 1.1 组织结构
```
├── Algorithm							通用算法
│   ├── qcqp.hpp
│   └── quadprog.hpp
├── Core								核心封装
│   ├── Core.hpp
│   ├── Dual.hpp
│   └── Matrix.hpp
├── FastMath.hpp						单文件引入整个库
└── Impl								部分内部实现
    ├── Common.hpp					
    ├── basic							Impl中基础件
    │   ├── complex.hpp
    │   └── rand.hpp
    ├── impl.hpp
    ├── linear							Impl中线性计算部分
    │   ├── equation.hpp
    │   └── matrix.hpp
    ├── nonlinear						Impl中非线性计算法部分
    │   ├── optimization.hpp
    │   └── solver.hpp
    └── util							Impl中一些通用工具
        ├── dbg.h
        ├── filter.hpp
        └── function.hpp
```

## 1.2 功能介绍
最为常用的矩阵运算支持在Core/Matrix.hpp中，其中实现了绝大部分常用的矩阵计算功能。包括但不限于：
1.矩阵间计算、矩阵与数之间计算
2.切片操作
3.基础线性代数运算：交换行列、取对角线、计算行列式、方阵求逆、伪逆、计算范数
4.矩阵分解：cholesky、全秩cholesky、LU、QR、SVD
5.矩阵函数：矩阵指数、矩阵对数
6.二次优化：QP、QCQP

# 2.FastMath库使用范例
FastMath的导入使用cmake包含FastMath库目录后，使用包含即可。
```
#include "FastMath.hpp"
```

FastMath的命令空间共有三个	FastMath, FastMath::Impl, FastMath::Algorithm。
建议使用方法为使用FastMath的别名fm来调用其中函数。
```
#include <iostream>
#include "FastMath.hpp"
int main() {
    fm::Matrix<double, 3, 4> A = {1,2,3,4};
    std::cout << A << std::endl;
    return 0;
}
```

## 2.1矩阵声明、赋值

```
#include <iostream>
#include "FastMath.hpp"
int main() {

    // 所有初始化均使用行优先，先写第一行，再写第二行....
    // 使用数组{1,2,3,4,5,6}列表初始化为array，进而传入矩阵A
    // 第一推荐使用此方法进行初始化，此方法会进行编译时检查，若输入数量与矩阵尺寸不符会触发编译时错误
    fm::Matrix<double, 2, 3> A({1,2,3,4,5,6});

    // 统一初始化，剩余不足的元素默认为0，若
    fm::Matrix<double, 3, 4> B = {1,2,3};

    // 统一初始化，超出维度部分自动忽略
    fm::Vector3d x{1,2,3,4};

    // 矩阵初始化为单一值
    fm::SquareMatrix<double, 2> C(3);

    // 默认初始化为全0
    fm::Matrix4d D;

    // 使用特殊矩阵初始化
    fm::Matrix3d X = fm::Matrix3d::Rand();

    // 使用矩阵B的第一列初始化向量b
    fm::Vector3d b = B.col(0);

    // 使用矩阵A的第一行初始化列向量a
    // 此处因为 A.row(0) 返回的不是矩阵对象，没有实现transpose方法，因此先调用 mat()将其转换为普通的矩阵，然后在调用transpose将行向量转为列向量，赋值给a。
    fm::Vector3d a = A.row(0).mat().transpose();
    
    // 所有矩阵均支持 operator<<
    std::cout << C << std::endl;



    return 0;
}
```

## 2.2元素索引访问
```
#include <iostream>
#include "FastMath.hpp"
int main() {

   fm::Matrix<double, 3, 3> R;

   double theta = fm::M_PI_F / 4;  // 建议使用FastMath内精度更高的常量定义，数值稳定性更好一些
   R(0,0) = cos(theta);        // R(0,0) 是第1行第1列
   R(1) = -sin(theta);           // R(1) 是整个矩阵的第2个元素，也就是第1行第2列
    R(1,0) = sin(theta);        // R(1,0) 是第2行第1列
    R(4) = cos(theta);            // R(4) 是整个举着你的第5个元素，也就是第2行第2列
    R(2,2) = 1;                    // R(2,2) 第3行第3列

    fm::Vector3d p = {1,2,3};

    fm::Matrix4d T;
    T.slice<3,3>(0,0) = R;      // 将R矩阵赋值给T矩阵的左上角3x3部分 slice<3,3>代表取3x3子块，(0,0)代表这个子块的第一个元素位于T(0,0)这个位置
    T.slice<3, 1>(0, 3) = p;    // 将R矩阵赋值给T矩阵的右上角3x1部分 slice<3,1>代表取3x1子块，(0,3)代表这个子块的第一个元素位于T(0,3)这个位置
    T(3, 3) = 1;                   // 将最后一个元素给1


    std::cout << T << std::endl;
    
    return 0;
}
```

## 2.3与数组交互
```
double R_array[3*3] = {1,2,3,4,56,6,7,8,9};
    fm::Matrix3d R;
    R.loadFrom(R_array);    //loadFrom 会对数组尺寸进行检查
    std::cout << R << std::endl;

    R.loadFromArrayPtr(R_array);
    std::cout << R << std::endl;   // loadFromArrayPtr 只需传入受元素的指针，不会对尺寸进行检查

    R = R * 2;

    double R_array_out[9];
    R.copyTo(R_array_out);     // 将R中内部数据拷贝到数组R_array_out中
    std::cout << fm::Impl::MatToString(&R_array_out[0], 3, 3) << std::endl; //MatToString是对普通数组的打印支持

    R.copyToArrayPtr(R_array_out);     //将R中内部数据拷贝到R_array_out指向的首地址处
    std::cout << fm::Impl::MatToString(&R_array_out[0], 3, 3) << std::endl; //MatToString是对普通数组的打印支持
```
## 2.4基本数学运算
```
#include <iostream>
#include "FastMath.hpp"
int main() {

    fm::Matrix3d A = {1,2,3,4,5,6};
    fm::Matrix<double, 3, 2> B;
    fm::Matrix<double, 3, 2> C;
    B.setIdentity();


    // 矩阵乘法，会在编译期检查维度匹配，维度错误编译报错
    auto D = A*B;

    // 矩阵相加
    C = C + D;
    C += D;

    fm::Matrix3d X;

    // 数乘矩阵
    X = 3*A;
    X = A*3;
    X *= 3;     // X = X*3
    X *= A;     // X = X*A;
    X = X * A;
    
    return 0;
}

```
## 2.5切片操作

```
#include <iostream>
#include "FastMath.hpp"
int main() {

    fm::Matrix4d T;
    T.setIdentity();

    // 将矩阵左上角写为单位正乘3
    T.slice<3,3>(0,0) = fm::Matrix3d::Identity() * 3;
    std::cout << T << std::endl;

    // 将左上角3x3矩阵赋值给A
    fm::Matrix3d A = T.slice<3,3>(0,0);
    std::cout << A << std::endl;

    // slice的结果使用mat()后会生成一个新矩阵，不会改变原本的矩阵
    T.slice<3,3>(0,0).mat() += 3;
    std::cout << T << std::endl;

    // slice后直接进行+=，会改变原本的矩阵中的数据
    T.slice<3,3>(0,0) -= 3;
    std::cout << T << std::endl;

    // 若希望使用切片之后的结果做矩阵运算，例如求逆，需要使用mat()将其转为矩阵
    fm::Matrix3d U, S, V;
    auto svd_success = T.slice<3, 3>(0,0).mat().SVD(U, S, V);
    if (svd_success)
        std::cout << U*S*V.transpose() << std::endl;
    else
        std::cout << "svd failure" << std::endl;

    return 0;
}
```

# 3.常用线性代数操作
## 3.1基本计算
### 3.1.1初等变换
```
fm::Matrix4d A = {1,2,3,4,5,6};
std::cout << A << std::endl;
A.swapCols(0, 1);       //交换12列
std::cout << A << std::endl;
A.swapRows(0, 3);       // 交换14行
std::cout << A << std::endl;
```

### 3.1.2 基本性质
```
fm::Matrix3d A = {1,2,3,4,5,6};
std::cout << A.rank() << std::endl; // 矩阵的秩
std::cout << A.norm() << std::endl;     // 矩阵的L2范数
std::cout << A.norm_one() << std::endl;     // 矩阵的L1范数
std::cout << A.norm_inf() << std::endl;     // 矩阵的L无穷范
std::cout << A.trace() << std::endl;    // 矩阵的迹
std::cout << A.det() << std::endl; // 矩阵的行列式
```

## 3.2.矩阵分解
### 3.2.1.LU分解
对矩阵A进行LU分解，A必须为方阵，分解为下三角矩阵L， 上三角矩阵U。此算法中未进行选主元操作是，因此数值稳定性不好。此分解计算量O(n^3)。
```
 fm::Matrix3d A = {1,2,3,4,5,6};
 fm::Matrix3d L, U;
 auto success = A.LU(L, U);
 
 std::cout << L << std::endl;
 std::cout << U << std::endl;
 std::cout << L*U << std::endl;
```
### 3.2.2.Cholesky分解
对矩阵A进行Cholesky分解，A必须为对称正定方阵，分解为A=Q*Q^T，Q为下三角矩阵。
全秩cholesky的数值稳定性更好，优先使用fullRankCholesky。分解失败时默认返回Q为零矩阵。
```
    fm::Matrix3d A = fm::Matrix3d::Rand();
    A += fm::Matrix3d::Identity() * 10;
    A = A + A.transpose();
  //auto Q = A.cholesky();
    auto Q = A.fullRankCholesky();
    std::cout << A << std::endl;
    std::cout << Q << std::endl;
    std::cout << Q*Q.transpose() << std::endl;
```
### 3.2.3.QR分解
对矩阵A进行QR分解，A矩阵需要行数大于列数，分解为A=Q*R。Q为正交矩阵，R为上三角矩阵。本算法在不满秩时数值稳定性不好，因此在使用时需要判断矩阵满秩。
```
 fm::Matrix<double, 4, 3> A = {1,2,3,4,5,6,1};
 fm::Matrix<double, 4, 4> Q;
 fm::Matrix<double, 4, 3> R;
 auto success = A.QR(Q, R);
 std::cout << success << std::endl;
 std::cout << A << std::endl;
 std::cout << Q << std::endl;
 std::cout << R << std::endl;
std::cout << Q*R << std::endl;
```

### 3.2.4.SVD分解
对矩阵A进行SVD分解，A可以是任意矩阵，该分解使用Household变换。分解为A=U*S*V^T。分解的结果中U和V为正交矩阵，S的对角线为奇异值。
```
   fm::Matrix<double, 4, 3> A = {1,2,3,4,5,6,1};
   fm::Matrix<double, 4, 4> U;
   fm::Matrix<double, 4, 3> S;
   fm::Matrix<double, 3, 3> V;
   auto success = A.SVD(U, S, V);
   std::cout << success << std::endl;
   std::cout << A << std::endl;
   std::cout << U << std::endl;
   std::cout << S << std::endl;
   std::cout << V << std::endl;
   std::cout << U*S*V.transpose() << std::endl;
  
```

### 3.2.5.矩阵求逆、伪逆
A.inv()可以对方阵求逆，求逆失败返回全零矩阵。
A.pinv()对矩阵求伪逆，使用SVD分解计算伪逆。
A.Geninv()对矩阵求伪逆，使用全秩Cholesky分解计算伪逆。
若无特殊需求，对方阵应使用inv()，而非pinv()或Geninv()。
求伪逆时优先使用pinv()
```
 fm::Matrix3d R = fm::Matrix3d::Rand();
 R += fm::Matrix3d::Identity() * 10;
 std::cout << R * R.inv() << std::endl;
 fm::Matrix<double, 4, 3> A = {1,2,3,4,5,6,1};
 std::cout << A*A.pinv() << std::endl;
 std::cout << A*A.Geninv() << std::endl;
```

# 4.矩阵函数
## 4.1.初等函数
### 4.1.1.矩阵的幂
power( int n )计算n个矩阵的连乘，但算次数比直接连乘少。例如对A*A*A*A，内部会计算一个A2=A*A，计算一次A4 = A2*A2。仅需要2次矩阵乘法，比直接连乘所需要的3次矩阵乘法少1次。
```
 fm::Matrix3d R(2);
 std::cout << R.power(3) << std::endl;
 std::cout << (R*R*R == R.power(3)) << std::endl;
```

### 4.1.2.矩阵的平方根
矩阵的平方根是要A=X*X中的X，但是满足条件的X是不唯一的，此处使用的是主平方根，也就是X的特征值全都有正的实部。算法内部使用牛顿迭代的改进版求解此方程。
```
fm::Matrix3d R(2);
R +=  fm::Matrix3d::Identity();
std::cout << R.std::sqrt() << std::endl;
std::cout << (R.std::sqrt() * R.std::sqrt() == R) << std::endl;
```

### 4.1.3.矩阵的指数
矩阵对数expm是控制系统和机器人学中非常常用的计算，此算法使用Pade近似计算expm
```
  fm::Matrix3d so3;
  so3(1, 2) = -fm::M_PI_F;
  so3(2, 1) = fm::M_PI_F;
  
  std::cout << so3.expm() << std::endl;
```

### 4.1.4.矩阵的对数
logm是expm的反函数，此算法使用反向缩放平方法计算logm的Pade近似。
```
 fm::Matrix3d so3;
 so3(1, 2) = -fm::M_PI_F / 2;
 so3(2, 1) = fm::M_PI_F /2 ;
 std::cout << so3 << std::endl;
 std::cout << so3.expm() << std::endl;
 std::cout << so3.expm().logm() << std::endl;
```