/******************************************************************************
  文 件 名   : Common.hpp
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

#ifndef FASTMATH_COMMON_HPP
#define FASTMATH_COMMON_HPP

#include <cmath>
#include <limits>
#include <array>
#include <cassert>
#include <cstdio>
#include <cstring>
#include <type_traits>
#include <iomanip>
#include <chrono>

//#define DEBUG_FASTMATH        // 开启调试打印

#ifdef DEBUG_FASTMATH
#define DEBUG_COUT std::cout
#define FILE_LINE "File : " << __FILE__ << " line : " << __LINE__
#endif

namespace FastMath::Impl{

    /* 常数定义 */
    constexpr double M_E_F			=  2.718281828459045090795598298427648842334747314453125;
    constexpr double M_LOG2E_F		=  1.442695040888963387004650940070860087871551513671875;
    constexpr double M_LOG10E_F		=  0.434294481903251761156781185491126962006092071533203125;
    constexpr double M_LN2_F		=  0.69314718055994528622676398299518041312694549560546875;
    constexpr double M_LN10_F		=  2.30258509299404590109361379290930926799774169921875;
    constexpr double M_PI_F			=  3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117068;
    constexpr double M_TWOPI_F		=  6.283185307179586476925286766559005768394338798750211641949889184615632812572417997256069650684234136;
    constexpr double M_PI_2_F		=  1.570796326794896619231321691639751442098584699687552910487472296153908203143104499314017412671058534;
    constexpr double M_PI_4_F		=  0.785398163397448309615660845819875721049292349843776455243736148076954101571552249657008706335529267;
    constexpr double M_3PI_4_F		=  2.356194490192344928846982537459627163147877049531329365731208444230862304714656748971026119006587801;
    constexpr double M_SQRTPI_F		=  1.772453850905515881919427556567825376987457275390625;
    constexpr double M_1_PI_F		=  0.318309886183790691216444201927515678107738494873046875;
    constexpr double M_2_PI_F		=  0.63661977236758138243288840385503135621547698974609375;
    constexpr double M_2_SQRTPI_F	=  1.1283791670955125585606992899556644260883331298828125;
    constexpr double M_DEG_TO_RAD_F	=  0.01745329251994329576923690768488612713442871888541725456097191440171009114603449443682241569634509482;
    constexpr double M_RAD_TO_DEG_F	=  57.29577951308232286464772187173366546630859375;
    constexpr double M_SQRT2_F		=  1.414213562373095048801688724209698078569671875376948073176679737990732478462107038850387534327641573;
    constexpr double M_SQRT1_2_F	=  0.7071067811865475244008443621048490392848359376884740365883398689953662392310535194251937671638207864;
    constexpr double M_LN2LO_F		=  1.90821484E-10f;
    constexpr double M_LN2HI_F		=  0.69314718f;
    constexpr double M_SQRT3_F		=  1.732050807568877293527446341505872366942805253810380628055806979451933016908800037081146186757248576;
    constexpr double M_IVLN10_F		=  0.434294481903251761156781185491126962006092071533203125;	// 1 / log(10)
    constexpr double M_LOG2_E_F		=  0.69314718055994528622676398299518041312694549560546875;
    constexpr double M_INVLN2_F		=  1.442695040888963387004650940070860087871551513671875; // 1 / log(2)

    constexpr double M_PI_PRECISE   =  3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117068;
    constexpr double M_DEG_TO_RAD   =  0.01745329251994329576923690768488612713442871888541725456097191440171009114603449443682241569634509482;
    constexpr double M_RAD_TO_DEG   =  57.29577951308232286464772187173366546630859375;

    constexpr double M_EPS   =  1e-10;
    #define FLT_EPSILON  __FLT_EPSILON__


    /**
     * Check value is finite
     */

    template<typename Type>
    inline bool is_finite(Type x) {
        return std::isfinite(x);
    }

    /**
     * Compare if two floating point numbers are equal
     *
     * NAN is considered equal to NAN and -NAN
     * INFINITY is considered equal INFINITY but not -INFINITY
     *
     * @param x right side of equality check
     * @param y left side of equality check
     * @param eps numerical tolerance of the check
     * @return true if the two values are considered equal, false otherwise
     */
    template<typename Type>
    inline bool isEqualF(const Type x, const Type y, const Type eps = Type(M_EPS))
    {
        return (std::fabs(x - y) <= eps)
               || (std::isnan(x) && std::isnan(y))
               || (std::isinf(x) && std::isinf(y) && std::isnan(x - y));
    }

    namespace detail
    {

        template<typename Floating>
        inline Floating wrap_floating(Floating x, Floating low, Floating high)
        {
            // already in range
            if (low <= x && x < high) {
                return x;
            }

            const auto range = high - low;
            const auto inv_range = Floating(1) / range; // should evaluate at compile time, multiplies below at runtime
            const auto num_wraps = std::floor((x - low) * inv_range);
            return x - range * num_wraps;
        }

    }  // namespace detail


    /**
     * Wrap single precision floating point value to stay in range [low, high)
     *
     * @param x input possibly outside of the range
     * @param low lower limit of the allowed range
     * @param high upper limit of the allowed range
     * @return wrapped value inside the range
     */
        inline float wrap(float x, float low, float high)
    {
        return detail::wrap_floating(x, low, high);
    }

    /**
     * Wrap double precision floating point value to stay in range [low, high)
     *
     * @param x input possibly outside of the range
     * @param low lower limit of the allowed range
     * @param high upper limit of the allowed range
     * @return wrapped value inside the range
     */
    inline double wrap(double x, double low, double high)
    {
        return detail::wrap_floating(x, low, high);
    }


    /**
 * Wrap integer value to stay in range [low, high)
 *
 * @param x input possibly outside of the range
 * @param low lower limit of the allowed range
 * @param high upper limit of the allowed range
 * @return wrapped value inside the range
 */
    template<typename Integer>
    inline Integer wrap(Integer x, Integer low, Integer high)
    {
        const auto range = high - low;

        if (x < low) {
            x += range * ((low - x) / range + 1);
        }

        return low + (x - low) % range;
    }


    /**
     * Wrap value in range [-π, π)
     */
    template<typename Type>
    inline Type wrap_pi(Type x)
    {
        return wrap(x, Type(-M_PI_PRECISE), Type(M_PI_PRECISE));
    }

    /**
     * Wrap value in range [0, 2π)
     */
    template<typename Type>
    inline Type wrap_2pi(Type x)
    {
        return wrap(x, Type(0), Type((2 * M_PI_PRECISE)));
    }

    /**
     * Unwrap value that was wrapped with range [low, high)
     *
     * @param[in] last_x Last unwrapped value
     * @param[in] new_x New value in range
     * @param low lower limit of the wrapping range
     * @param high upper limit of the wrapping range
     * @return New unwrapped value
     */
    template<typename Type>
    inline Type unwrap(const Type last_x, const Type new_x, const Type low, const Type high)
    {
        return last_x + wrap(new_x - last_x, low, high);
    }


    /**
     * Unwrap value with range [-π, π)
     *
     * @param[in] last_angle Last unwrapped angle [rad]
     * @param[in] new_angle New angle in [-pi, pi] [rad]
     * @return New unwrapped angle [rad]
     */
    template<typename Type>
    inline Type unwrap_pi(const Type last_angle, const Type new_angle)
    {
        return unwrap(last_angle, new_angle, Type(-M_PI_PRECISE), Type(M_PI_PRECISE));
    }

    /**
     * Type-safe sign/signum function
     *
     * @param[in] val Number to take the sign from
     * @return -1 if val < 0, 0 if val == 0, 1 if val > 0
     */
    template<typename T>
    inline int sign(T val)
    {
        return (T(0) < val) - (val < T(0));
    }

    template <typename T>
    inline auto is_positive_infinity(T val) -> bool
    {
        return val == std::numeric_limits<T>::infinity();
    }

    template <typename T>
    inline auto is_negative_infinity(T val) -> bool
    {
        return val == -std::numeric_limits<T>::infinity();
    }

    template <typename T>
    inline auto is_nan(T val) -> bool
    {
        return std::isnan(val);
    }

    template <typename T>
    inline auto is_safe(T val) -> bool
    {
        bool safe = true;
        safe &= !is_positive_infinity(val);
        safe &= !is_negative_infinity(val);
        safe &= !is_nan(val);

        return safe;
    }





}



#endif //FASTMATH_COMMON_HPP
