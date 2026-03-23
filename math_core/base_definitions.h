#ifndef __base_definitions_h__
#define __base_definitions_h__

#include <math.h>

/**
 * @namespace base_math
 * @brief Contains mathematical constants and utility functions for numerical computations.
 */
namespace base_math {
    /**
     * @brief Mathematical constant PI (�).
     * @tparam T Numeric type.
     */
    template <typename T>
    constexpr T PI = T(3.14159265358979323846);

    /**
     * @brief Mathematical constant 2�.
     * @tparam T Numeric type.
     */
    template <typename T>
    constexpr T TWO_PI = T(6.28318530717958647692);

    /**
     * @brief Mathematical constant �/2.
     * @tparam T Numeric type.
     */
    template <typename T>
    constexpr T HALF_PI = T(1.57079632679489661923);

    /**
     * @brief Conversion factor from degrees to radians (�/180).
     * @tparam T Numeric type.
     */
    template <typename T>
    constexpr T PID180 = T(0.01745329251994329577);

    /**
     * @brief Conversion factor from radians to degrees (180/�).
     * @tparam T Numeric type.
     */
    template <typename T>
    constexpr T PID180INV = T(57.2957795130823208768);

    /**
     * @brief Small epsilon value for floating-point comparisons.
     * @tparam T Numeric type.
     */
    template <typename T>
    constexpr T EPSILON = T(1e-6);

    /**
     * @brief Square of the epsilon value.
     * @tparam T Numeric type.
     */
    template <typename T>
    constexpr T EPSILON_SQ = T(1e-12);

    /**
     * @brief Tolerance value for floating-point comparisons.
     * @tparam T Numeric type.
     */
    template <typename T>
    constexpr T TOLLERANCE = T(EPSILON<T>);

    /**
     * @brief Square of the tolerance value.
     * @tparam T Numeric type.
     */
    template <typename T>
    constexpr T TOLLERANCE_SQ = T(EPSILON_SQ<T>);

    /**
     * @brief Returns the sign of a value.
     * @tparam T Numeric type.
     * @param val Value to check.
     * @return 1 if positive, -1 if negative, 0 if zero.
     */
    template <typename T>
    int sign(T val) {
        return (T(0) < val) - (val < T(0));
    }

    /**
     * @brief Converts degrees to radians.
     * @tparam T Numeric type.
     * @param deg Angle in degrees.
     * @return Angle in radians.
     */
    template <typename T>
    T deg_to_rad(T deg) {
        return T(deg * PID180<T>);
    }
#define dtr deg_to_rad

    /**
     * @brief Converts radians to degrees.
     * @tparam T Numeric type.
     * @param rad Angle in radians.
     * @return Angle in degrees.
     */
    template <typename T>
    T rad_to_deg(T rad) {
        return T(rad * PID180INV<T>);
    }
#define rtd rad_to_deg

    /**
     * @brief Approximates a real number to the nearest integer.
     * @tparam T Numeric type.
     * @param d Real value.
     * @return Nearest integer.
     */
    template <typename T>
    int dti(T d) {
        return (int)floor(d + T(0.5));
    }

    /**
     * @brief Returns the sign of a value (alternative to sign).
     * @tparam T Numeric type.
     * @param val Value to check.
     * @return 1 if positive, -1 if negative, 0 if zero.
     */
    template <typename T> 
    int sgn(T val) {
        return (T(0) < val) - (val < T(0));
    }
}

#endif // __base_definitions_h__
