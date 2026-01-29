#ifndef __algebra_h__
#define __algebra_h__

#include <tuple>

namespace base_math {
    /**
     * @brief Solves the quadratic equation a*x^2 + b*x + c = 0.
     *
     * This function computes the real roots of a quadratic equation with coefficients a, b, and c.
     * If the discriminant is negative, the equation has no real roots and the function returns false.
     * Otherwise, the roots are stored in the provided tuple 'res' as (x1, x2).
     *
     * @tparam T Numeric type for coefficients and roots (e.g., float, double).
     * @param a Coefficient of x^2.
     * @param b Coefficient of x.
     * @param c Constant term.
     * @param res Tuple to store the two real roots (x1, x2) if they exist.
     * @return true if real roots exist, false otherwise.
     */
    template <typename T>
    bool solve_quadratic(T a, T b, T c, std::tuple<T, T>& res) {
        T d = b * b - 4 * a * c;
        if (d < 0)
            return false;
        T x1 = (-b + (T)sqrt(d)) / (T(2) * a);
        T x2 = (-b - (T)sqrt(d)) / (T(2) * a);
        if (x1>x2)
			std::swap(x1, x2);
        res = std::make_tuple(x1, x2);
        return true;
    }

    // Solve 3x3 system A x = b using Cholesky (A symmetric 3x3).
    // - A : basematrix<T,3,3> or raw row-major 3x3 via operator[]
    // - b : basematrix<T,3,1> RHS
    // - x : basematrix<T,3,1> output
    // Returns true on success.
    template <typename T>
    bool cholesky_solve_3x3(const basematrix<T, 3, 3>& A,
        const basematrix<T, 3, 1>& b,
        basematrix<T, 3, 1>& x)
    {
        // lower-triangular L (row-major 3x3, store only lower part)
        T L[9] = { 0.0f };
        const T reg = TOLLERANCE<T>; // small regularization value

        // Cholesky decomposition A = L * L^T
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j <= i; ++j) {
                T sum = 0.0f;
                for (int k = 0; k < j; ++k)
                    sum += L[i * 3 + k] * L[j * 3 + k];

                T a_ij = A[i * 3 + j];
                T val = a_ij - sum;

                if (i == j) {
                    // diagonal
                    if (val <= 0.0f) {
                        // Regularize to avoid sqrt of non-positive number.
                        // This preserves numeric stability for near-singular cases.
                        val += reg;
                    }
                    if (val <= 0.0f) // still non-positive -> fail
                        return false;
                    L[i * 3 + j] = sqrtf(val);
                }
                else {
                    // off-diagonal
                    T ljj = L[j * 3 + j];
                    if (fabsf(ljj) < TOLLERANCE<T>) {
                        // ill-conditioned, fail
                        return false;
                    }
                    L[i * 3 + j] = val / ljj;
                }
            }
        }

        // Forward substitution: L * y = b
        T y[3] = { 0.0f };
        for (int i = 0; i < 3; ++i) {
            T s = b[i];
            for (int k = 0; k < i; ++k)
                s -= L[i * 3 + k] * y[k];
            T lii = L[i * 3 + i];
            if (fabsf(lii) < TOLLERANCE<T>) return false;
            y[i] = s / lii;
        }

        // Backward substitution: L^T * x = y
        for (int i = 2; i >= 0; --i) {
            T s = y[i];
            for (int k = i + 1; k < 3; ++k)
                s -= L[k * 3 + i] * x[k];
            T lii = L[i * 3 + i];
            if (fabsf(lii) < TOLLERANCE<T>) return false;
            x[i] = s / lii;
        }

        return true;
    }




};
#endif // __algebra_h__
