#include <cassert>
#include <iostream>
#include <cmath>
#include "matrix.h"

using namespace base_math;

// Helper function to compare floating point values
template<typename T>
bool approxEqual(T a, T b, T tolerance = T(1e-6)) {
    return std::fabs(a - b) < tolerance;
}

void test_default_constructor() {
    fmat3 m;
    for (size_t i = 0; i < 9; ++i) {
        assert(m[i] == 0.0f);
    }
    std::cout << "Default constructor test passed\n";
}

void test_copy_constructor() {
    fmat3 m1({ 1, 2, 3, 4, 5, 6, 7, 8, 9 });
    fmat3 m2(m1);

    for (size_t i = 0; i < 9; ++i) {
        assert(m2[i] == m1[i]);
    }
    std::cout << "Copy constructor test passed\n";
}

void test_array_constructor() {
    float data[] = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
    fmat3 m(data);

    for (size_t i = 0; i < 9; ++i) {
        assert(m[i] == data[i]);
    }
    std::cout << "Array constructor test passed\n";
}

void test_initial_value_constructor() {
    fmat3 m(5.0f);

    for (size_t i = 0; i < 9; ++i) {
        assert(m[i] == 5.0f);
    }
    std::cout << "Initial value constructor test passed\n";
}

void test_vector_constructor() {
    std::vector<float> v = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
    fmat3 m(v);

    for (size_t i = 0; i < 9; ++i) {
        assert(m[i] == v[i]);
    }
    std::cout << "Vector constructor test passed\n";
}

void test_rows_cols() {
    fmat3 m;
    assert(m.rows() == 3);
    assert(m.cols() == 3);

    basematrix<float, 2, 4> m2;
    assert(m2.rows() == 2);
    assert(m2.cols() == 4);

    std::cout << "rows() and cols() test passed\n";
}

void test_element_access() {
    fmat3 m({ 1, 2, 3, 4, 5, 6, 7, 8, 9 });

    // Test operator()
    assert(m(0, 0) == 1.0f);
    assert(m(0, 1) == 2.0f);
    assert(m(1, 0) == 4.0f);
    assert(m(2, 2) == 9.0f);

    // Test operator[] 
    assert(m[0] == 1.0f);
    assert(m[4] == 5.0f);
    assert(m[8] == 9.0f);

    // Test modification
    m(1, 1) = 100.0f;
    assert(m(1, 1) == 100.0f);
    assert(m[4] == 100.0f);

    std::cout << "Element access test passed\n";
}

void test_assignment_operator() {
    fmat3 m1({ 1, 2, 3, 4, 5, 6, 7, 8, 9 });
    fmat3 m2;

    m2 = m1;

    for (size_t i = 0; i < 9; ++i) {
        assert(m2[i] == m1[i]);
    }

    // Ensure deep copy
    m1[0] = 999.0f;
    assert(m2[0] == 1.0f);

    std::cout << "Assignment operator test passed\n";
}

void test_addition_operators() {
    fmat3 m1({ 1, 2, 3, 4, 5, 6, 7, 8, 9 });
    fmat3 m2({ 9, 8, 7, 6, 5, 4, 3, 2, 1 });

    // Test operator+=
    fmat3 m3 = m1;
    m3 += m2;
    for (size_t i = 0; i < 9; ++i) {
        assert(m3[i] == m1[i] + m2[i]);
    }

    // Test operator+
    fmat3 m4 = m1 + m2;
    for (size_t i = 0; i < 9; ++i) {
        assert(m4[i] == m1[i] + m2[i]);
    }

    std::cout << "Addition operators test passed\n";
}

void test_subtraction_operators() {
    fmat3 m1({ 9, 8, 7, 6, 5, 4, 3, 2, 1 });
    fmat3 m2({ 1, 2, 3, 4, 5, 6, 7, 8, 9 });

    // Test operator-=
    fmat3 m3 = m1;
    m3 -= m2;
    for (size_t i = 0; i < 9; ++i) {
        assert(approxEqual(m3[i], m1[i] - m2[i]));
    }

    // Test operator-
    fmat3 m4 = m1 - m2;
    for (size_t i = 0; i < 9; ++i) {
        assert(approxEqual(m4[i], m1[i] - m2[i]));
    }

    std::cout << "Subtraction operators test passed\n";
}

void test_unary_negation() {
    fmat3 m1({ 1, -2, 3, -4, 5, -6, 7, -8, 9 });
    fmat3 m2 = -m1;

    for (size_t i = 0; i < 9; ++i) {
        assert(m2[i] == -m1[i]);
    }

    std::cout << "Unary negation test passed\n";
}

void test_scalar_multiplication() {
    fmat3 m1({ 1, 2, 3, 4, 5, 6, 7, 8, 9 });

    // Test operator*=
    fmat3 m2 = m1;
    m2 *= 2.0f;
    for (size_t i = 0; i < 9; ++i) {
        assert(m2[i] == m1[i] * 2.0f);
    }

    // Test operator* (matrix * scalar)
    fmat3 m3 = m1 * 3.0f;
    for (size_t i = 0; i < 9; ++i) {
        assert(m3[i] == m1[i] * 3.0f);
    }

    // Test operator* (scalar * matrix)
    fmat3 m4 = 4.0f * m1;
    for (size_t i = 0; i < 9; ++i) {
        assert(m4[i] == m1[i] * 4.0f);
    }

    // Test scale method
    fmat3 m5 = m1;
    m5.scale(5.0f);
    for (size_t i = 0; i < 9; ++i) {
        assert(m5[i] == m1[i] * 5.0f);
    }

    std::cout << "Scalar multiplication test passed\n";
}

void test_scalar_division() {
    fmat3 m1({ 2, 4, 6, 8, 10, 12, 14, 16, 18 });

    // Test operator/=
    fmat3 m2 = m1;
    m2 /= 2.0f;
    for (size_t i = 0; i < 9; ++i) {
        assert(approxEqual(m2[i], m1[i] / 2.0f));
    }

    // Test operator/
    fmat3 m3 = m1 / 4.0f;
    for (size_t i = 0; i < 9; ++i) {
        assert(approxEqual(m3[i], m1[i] / 4.0f));
    }

    std::cout << "Scalar division test passed\n";
}

void test_load_identity() {
    fmat3 m({ 1, 2, 3, 4, 5, 6, 7, 8, 9 });
    m.loadIdentity();

    assert(m(0, 0) == 1.0f);
    assert(m(1, 1) == 1.0f);
    assert(m(2, 2) == 1.0f);

    for (size_t r = 0; r < 3; ++r) {
        for (size_t c = 0; c < 3; ++c) {
            if (r == c) {
                assert(m(r, c) == 1.0f);
            }
            else {
                assert(m(r, c) == 0.0f);
            }
        }
    }

    std::cout << "loadIdentity() test passed\n";
}

void test_swap_rows() {
    fmat3 m({ 1, 2, 3, 4, 5, 6, 7, 8, 9 });

    m.swap_rows(0, 2);

    assert(m(0, 0) == 7.0f);
    assert(m(0, 1) == 8.0f);
    assert(m(0, 2) == 9.0f);
    assert(m(2, 0) == 1.0f);
    assert(m(2, 1) == 2.0f);
    assert(m(2, 2) == 3.0f);
    assert(m(1, 0) == 4.0f); // Row 1 unchanged

    std::cout << "swap_rows() test passed\n";
}

void test_swap_cols() {
    fmat3 m({ 1, 2, 3, 4, 5, 6, 7, 8, 9 });

    m.swap_cols(0, 2);

    assert(m(0, 0) == 3.0f);
    assert(m(0, 2) == 1.0f);
    assert(m(1, 0) == 6.0f);
    assert(m(1, 2) == 4.0f);
    assert(m(2, 0) == 9.0f);
    assert(m(2, 2) == 7.0f);
    assert(m(0, 1) == 2.0f); // Column 1 unchanged

    std::cout << "swap_cols() test passed\n";
}

void test_matrix_multiplication() {
    // Test 3x3 * 3x3
    fmat3 m1({ 1, 2, 3, 4, 5, 6, 7, 8, 9 });
    fmat3 m2({ 9, 8, 7, 6, 5, 4, 3, 2, 1 });

    fmat3 result = m1 * m2;

    // Expected: [30, 24, 18, 84, 69, 54, 138, 114, 90]
    assert(result(0, 0) == 30.0f);
    assert(result(0, 1) == 24.0f);
    assert(result(0, 2) == 18.0f);
    assert(result(1, 0) == 84.0f);
    assert(result(1, 1) == 69.0f);
    assert(result(1, 2) == 54.0f);

    // Test identity property
    fmat3 identity;
    identity.loadIdentity();
    fmat3 result2 = m1 * identity;

    for (size_t i = 0; i < 9; ++i) {
        assert(approxEqual(result2[i], m1[i]));
    }

    std::cout << "Matrix multiplication test passed\n";
}

void test_transpose() {
    fmat3 m({ 1, 2, 3, 4, 5, 6, 7, 8, 9 });

    fmat3 t = m.transpose();

    assert(t(0, 0) == 1.0f);
    assert(t(0, 1) == 4.0f);
    assert(t(0, 2) == 7.0f);
    assert(t(1, 0) == 2.0f);
    assert(t(1, 1) == 5.0f);
    assert(t(1, 2) == 8.0f);
    assert(t(2, 0) == 3.0f);
    assert(t(2, 1) == 6.0f);
    assert(t(2, 2) == 9.0f);

    // Test non-square matrix
    basematrix<float, 2, 3> rect({ 1, 2, 3, 4, 5, 6 });
    auto rect_t = rect.transpose();

    assert(rect_t.rows() == 3);
    assert(rect_t.cols() == 2);
    assert(rect_t(0, 0) == 1.0f);
    assert(rect_t(1, 0) == 2.0f);
    assert(rect_t(2, 0) == 3.0f);
    assert(rect_t(0, 1) == 4.0f);
    assert(rect_t(1, 1) == 5.0f);
    assert(rect_t(2, 1) == 6.0f);

    std::cout << "Transpose test passed\n";
}

void test_equality_operator() {
    fmat3 m1({ 1, 2, 3, 4, 5, 6, 7, 8, 9 });
    fmat3 m2({ 1, 2, 3, 4, 5, 6, 7, 8, 9 });
    fmat3 m3({ 1, 2, 3, 4, 5.1f, 6, 7, 8, 9 });

    assert(m1 == m2);
    assert(!(m1 == m3));

    std::cout << "Equality operator test passed\n";
}

void test_rectangular_matrices() {
    // Test 2x3 matrix
    basematrix<float, 2, 3> m1({ 1, 2, 3, 4, 5, 6 });

    assert(m1.rows() == 2);
    assert(m1.cols() == 3);
    assert(m1(0, 0) == 1.0f);
    assert(m1(1, 2) == 6.0f);

    // Test 2x3 * 3x2 multiplication
    basematrix<float, 3, 2> m2({ 1, 2, 3, 4, 5, 6 });
    auto result = m1 * m2;

    assert(result.rows() == 2);
    assert(result.cols() == 2);

    std::cout << "Rectangular matrices test passed\n";
}

void test_4x4_matrices() {
    fmat4 m1({ 1, 2, 3, 4, 2, 3, 4, 5, 3, 4, 5, 6, 4, 5, 6, 7 });
    fmat4 m2({ 7, 6, 5, 4, 6, 5, 4, 3, 5, 4, 3, 2, 4, 3, 2, 1 });

    fmat4 result = m1 * m2;

    // Verify dimensions
    assert(result.rows() == 4);
    assert(result.cols() == 4);

    // Test identity
    fmat4 identity;
    identity.loadIdentity();
    fmat4 result2 = m1 * identity;

    for (size_t i = 0; i < 16; ++i) {
        assert(approxEqual(result2[i], m1[i]));
    }

    std::cout << "4x4 matrices test passed\n";
}

void test_double_precision() {
    dmat3 m1({ 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0 });
    dmat3 m2({ 9.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0 });

    dmat3 result = m1 + m2;

    for (size_t i = 0; i < 9; ++i) {
        assert(approxEqual(result[i], 10.0));
    }

    std::cout << "Double precision test passed\n";
}

int run_all_matrix_tests() {
    std::cout << "\n=== Running Matrix Tests ===\n";

    test_default_constructor();
    test_copy_constructor();
    test_array_constructor();
    test_initial_value_constructor();
    test_vector_constructor();
    test_rows_cols();
    test_element_access();
    test_assignment_operator();
    test_addition_operators();
    test_subtraction_operators();
    test_unary_negation();
    test_scalar_multiplication();
    test_scalar_division();
    test_load_identity();
    test_swap_rows();
    test_swap_cols();
    test_matrix_multiplication();
    test_transpose();
    test_equality_operator();
    test_rectangular_matrices();
    test_4x4_matrices();
    test_double_precision();

    std::cout << "\n=== All Matrix Tests Passed ===\n\n";
    return 0;
}
