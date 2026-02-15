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
                    L[i * 3 + j] = T(sqrt(val));
                }
                else {
                    // off-diagonal
                    T ljj = L[j * 3 + j];
                    if (T(fabs(ljj)) < TOLLERANCE<T>) {
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
            if (T(fabs(lii)) < TOLLERANCE<T>) return false;
            y[i] = s / lii;
        }

        // Backward substitution: L^T * x = y
        for (int i = 2; i >= 0; --i) {
            T s = y[i];
            for (int k = i + 1; k < 3; ++k)
                s -= L[k * 3 + i] * x[k];
            T lii = L[i * 3 + i];
            if (T(fabs(lii)) < TOLLERANCE<T>) return false;
            x[i] = s / lii;
        }

        return true;
    }



    /**
     * @brief Recovers second fundamental form coefficients (e, f, g) from sampled data via least squares.
     *
     * This routine solves a small overdetermined linear system to estimate the per-triangle
     * second fundamental form entries '(e, f, g)' from preassembled measurements stored in
     * 'second_fundamental'. The system is built as:
     *
     *   E * [e f g]^T ? N
     *
     * where:
     *   - 'N' is a '6x1' vector of measured normal differences (finite differences in normal space),
     *   - 'E' is a '6x3' matrix encoding local '(du, dv)' edge parameter differences.
     *
     * The least-squares solution is obtained via the normal equations:
     *
     *   (E^T * E) * [e f g]^T = E^T * N
     *
     * and solved using the Cholesky-based solver 'cholesky_solve_3x3'. If the normal matrix
     * 'E^T * E' is near-singular or ill-conditioned, a tiny diagonal regularization
     * ('TOLLERANCE<T>') is applied (Tikhonov regularization) and the system is solved again.
     *
     * Input layout (second_fundamental):
     *   - 'second_fundamental[0..5]'   : N (6 components of measured normal differences)
     *   - 'second_fundamental[6..11]'  : three '(du, dv)' pairs, one per triangle edge:
     *        [du1, dv1, du2, dv2, du3, dv3]
     *
     * Internal assembly:
     *   - 'N' is assembled as a 'basematrix<T, 6, 1>' from 'second_fundamental[0..5]'.
     *   - 'E' is assembled as a 'basematrix<T, 6, 3>' from the three '(du, dv)' pairs
     *     following the documented row pattern:
     *
     *       [du1 dv1 0  0 du1 dv1;
     *        du2 dv2 0  0 du2 dv2;
     *        du3 dv3 0  0 du3 dv3]
     *
     *     (compactly encoded in row-major form through the constructor initializer list).
     *
     * Output layout (e_f_g):
     *   - On return, 'e_f_g[0]' = e, 'e_f_g[1]' = f, 'e_f_g[2]' = g.
     *
     * @tparam T
     *   Scalar numeric type used for all computations (e.g., 'float', 'double').
     *
     * @param second_fundamental
     *   Pointer to a contiguous array holding the per-triangle second fundamental form
     *   measurements. It must contain at least 12 entries:
     *     - indices '0..5'  : normal-difference samples forming the RHS vector N;
     *     - indices '6..11' : three '(du, dv)' pairs describing local parameter-space edge
     *                         directions used to build the design matrix E.
     *
     * @param e_f_g
     *   Pointer to an array of 3 elements where the estimated second fundamental form
     *   coefficients will be written:
     *     - 'e_f_g[0]' := e
     *     - 'e_f_g[1]' := f
     *     - 'e_f_g[2]' := g
     *
     * @note
     *   This helper assumes that the input data has already been assembled consistently
     *   in the correct local parameterization of the triangle. It does not perform any
     *   validation of geometry or parameterization consistency; it only performs the
     *   small linear-algebra solve.
     *
     * @warning
     *   If both the unregularized and the regularized Cholesky solves fail (e.g., in
     *   highly degenerate geometric configurations), 'e_f_g' will still be written
     *   with the contents of the last attempted solution, which may be undefined or
     *   numerically unstable. Callers may want to extend this function to propagate
     *   a success/failure flag if they need robust error handling.
     */

    template <typename T>
    void leastSquares(const T* second_fundamental, T* e_f_g)
    {
        // N is the 6x1 RHS vector storing measured normal differences second_fundamental[0..5]
        basematrix<T, 6, 1> N({ second_fundamental[0],second_fundamental[1],second_fundamental[2],
                        second_fundamental[3],second_fundamental[4],second_fundamental[5] });
        // E is the 6x3 matrix storing the local edge differences second_fundamental[6..17]
        // arranged as rows: [du1 dv1 0 0 du1 dv1; du2 dv2 0 0 du2 dv2; du3 dv3 0 0 du3 dv3]
        basematrix<T, 6, 3> E({ second_fundamental[6], second_fundamental[7],0,0,second_fundamental[6],second_fundamental[7] ,
                       second_fundamental[8], second_fundamental[9],0,0,second_fundamental[8],second_fundamental[9] ,
                      second_fundamental[10], second_fundamental[11],0,0,second_fundamental[10],second_fundamental[11] });

        // Transpose of E (3x6)
        basematrix<T, 3, 6> Trans_E = E.transpose();
        // Compute normal equations components
        basematrix<T, 3, 3> tr_x_e = Trans_E * E;           // yields the symmetric 3x3 normal matrix -> E^T * E (3x3)
        basematrix<T, 3, 1> rhs = Trans_E * N;              // accumulates  E^T * N (3x1)
        basematrix<T, 3, 1> sol;                            // will hold [ee ff gg]^T
        // Solving tr_x_e * sol = rhs recovers[ee, ff, gg] ^ T.
        // Cholesky is used because E^T*E is symmetric positive semi-definite.
        if (!cholesky_solve_3x3<T>(tr_x_e, rhs, sol)) {
            // If the first solve fails(e.g., because E^T*E is near-singular), 
            // the diagonal gets a tiny TOLLERANCE<double> bump, 
            // effectively applying Tikhonov regularization, and the solve is retried.
            basematrix<T, 3, 3> reg = tr_x_e;
            for (int i = 0; i < 3; ++i) reg[i * 3 + i] += TOLLERANCE<T>;
            cholesky_solve_3x3<T>(reg, rhs, sol);
        }

        e_f_g[0] = sol[0];
        e_f_g[1] = sol[1];
        e_f_g[2] = sol[2];
    }





};
#endif // __algebra_h__
