/**
 * @file vector_tests.cpp
 * @brief Unit tests for the small math library's vector and matrix types.
 *
 * This file implements a small, self-contained test harness for the types
 * defined in `vector.h` (namespace `base_math`). Tests use simple helper
 * assertion functions (`expect_true`, `expect_eq`, `expect_eq_int`) that
 * print failures to `std::cerr` and increment the global `g_failures`
 * counter. The entry point `run_all_vector_tests()` runs all individual
 * test cases and returns a process-style status code (0 on success, 1 on
 * failure).
 *
 * The tests exercise:
 *  - Constructors and component accessors for `basevector<T, C>`
 *  - Length, squared length and normalization
 *  - Dot and cross products (and operator overloads for vec3)
 *  - Arithmetic operators (+, -, *, /, +=)
 *  - Indexing bounds behavior (const vs non-const)
 *  - Matrix * vector multiplication using `basematrix<T, R, C>`
 *
 * These tests are intentionally small and deterministic so they can be run
 * as part of a unit-test step in a build or CI pipeline.
 */

#include <iostream>
#include <cmath>
#include <stdexcept>

#include "vector.h"

using namespace base_math;

/* Global failure counter used by the lightweight test harness. Each failed
 * expectation increments this counter so the runner can return a non-zero
 * exit code when any test fails.
 */
static int g_failures = 0;

/**
 * @brief Check that a boolean condition is true.
 *
 * If `cond` is false, prints a failure message to `std::cerr` and increments
 * `g_failures`.
 *
 * @param cond The condition expected to be true.
 * @param msg   A short message describing the expectation being checked.
 */
inline void expect_true(bool cond, const char* msg) {
    if (!cond) {
        std::cerr << "FAIL: " << msg << std::endl;
        ++g_failures;
    }
}

/**
 * @brief Check that two floating-point values are approximately equal.
 *
 * Uses an absolute tolerance `tol` to allow for small numerical differences.
 * On failure, prints both the expected and actual values and increments
 * `g_failures`.
 *
 * @param a   The actual value produced by the code under test.
 * @param b   The expected value.
 * @param msg A short message describing the expectation.
 * @param tol Allowed absolute difference (default 1e-5).
 */
inline void expect_eq(float a, float b, const char* msg, float tol = 1e-5f) {
    if (std::fabs(a - b) > tol) {
        std::cerr << "FAIL: " << msg << " (expected " << b << ", got " << a << ")" << std::endl;
        ++g_failures;
    }
}

/**
 * @brief Check that two integers are equal.
 *
 * Prints a failure message and increments `g_failures` when `a != b`.
 *
 * @param a   The actual integer value.
 * @param b   The expected integer value.
 * @param msg Description of the expectation being validated.
 */
inline void expect_eq_int(int a, int b, const char* msg) {
    if (a != b) {
        std::cerr << "FAIL: " << msg << " (expected " << b << ", got " << a << ")" << std::endl;
        ++g_failures;
    }
}

/**
 * @brief Test constructors and accessor methods on vector types.
 *
 * Verifies:
 *  - Multi-component constructors (e.g., vec3(x,y,z))
 *  - Single-value broadcast constructor (e.g., vec4(scalar))
 *  - Construction from a smaller vector plus an appended component
 *  - Pointer conversion exposes underlying components in expected order
 */
void test_constructors_and_accessors() {
    // 3-component constructor
    basevector<float, 3> a(1.f, 2.f, 3.f);
    expect_eq(a.x(), 1.f, "a.x");
    expect_eq(a.y(), 2.f, "a.y");
    expect_eq(a.z(), 3.f, "a.z");
    expect_eq_int((int)a.size(), 3, "size() for vec3");

    // single-value constructor
    basevector<float, 4> b(5.f);
    expect_eq(b.x(), 5.f, "b.x single ctor");
    expect_eq(b.w(), 5.f, "b.w single ctor");
    expect_eq_int((int)b.size(), 4, "size() for vec4");

    // construct from smaller vector + extra component
    basevector<float, 2> v2(7.f, 8.f);
    basevector<float, 3> v3(v2, 9.f);
    expect_eq(v3.x(), 7.f, "v3.x from v2");
    expect_eq(v3.y(), 8.f, "v3.y from v2");
    expect_eq(v3.z(), 9.f, "v3.z from appended value");

    // pointer conversion
    float* pdata = (float*)a;
    expect_eq(pdata[0], 1.f, "pointer conversion points to x");
    expect_eq(pdata[1], 2.f, "pointer conversion points to y");
}

/**
 * @brief Test length and normalization helpers.
 *
 * Confirms length squared and length calculations and that `normalize()`
 * scales a vector to unit length. Uses the classical (3,4,0) example to
 * verify expected numeric results.
 */
void test_length_and_normalize() {
    basevector<float, 3> v(3.f, 4.f, 0.f);
    expect_eq(v.length_sq(), 25.f, "length_sq of (3,4,0)");
    expect_eq(v.length(), 5.f, "length of (3,4,0)");

    v.normalize();
    expect_eq(v.length(), 1.f, "normalized length is 1");
    // normalized (3,4,0) -> (0.6, 0.8, 0)
    expect_eq(v.x(), 0.6f, "normalized x");
    expect_eq(v.y(), 0.8f, "normalized y");
}

/**
 * @brief Test dot and cross product operations.
 *
 * Validates:
 *  - Cross product produces orthogonal vector (i x j = k)
 *  - Operator* behavior for vec3 (expected to be cross product)
 *  - Dot product returns the scalar product
 */
void test_dot_and_cross() {
    basevector<float, 3> i(1.f, 0.f, 0.f);
    basevector<float, 3> j(0.f, 1.f, 0.f);
    basevector<float, 3> k = i.cross(j);
    expect_eq(k.x(), 0.f, "i x j -> x");
    expect_eq(k.y(), 0.f, "i x j -> y");
    expect_eq(k.z(), 1.f, "i x j -> z");

    // operator* should be cross for vec3
    basevector<float, 3> k2 = i * j;
    expect_eq(k2.z(), 1.f, "operator* cross z");

    basevector<float, 3> a(1.f, 2.f, 3.f);
    basevector<float, 3> b(4.f, 5.f, 6.f);
    expect_eq(a.dot(b), 32.f, "dot product 1,2,3 . 4,5,6");
}

/**
 * @brief Test arithmetic operators and combinations.
 *
 * Exercises:
 *  - Vector addition/subtraction and their in-place variants
 *  - Scalar division and multiplication
 *  - Component-wise correctness after operations
 */
void test_operators_and_arithmetic() {
    basevector<float, 3> a(1.f, 1.f, 1.f);
    basevector<float, 3> b(2.f, 3.f, 4.f);

    auto c = a + b;
    expect_eq(c.x(), 3.f, "+ operator x");
    expect_eq(c.y(), 4.f, "+ operator y");
    expect_eq(c.z(), 5.f, "+ operator z");

    a += b;
    expect_eq(a.x(), 3.f, "+= operator x");
    expect_eq(a.y(), 4.f, "+= operator y");
    expect_eq(a.z(), 5.f, "+= operator z");

    auto d = c - b;
    expect_eq(d.x(), 1.f, "- operator x");
    expect_eq(d.y(), 1.f, "- operator y");
    expect_eq(d.z(), 1.f, "- operator z");

    auto e = basevector<float, 3>(2.f, 4.f, 6.f) / 2.f;
    expect_eq(e.x(), 1.f, "/ operator x");
    expect_eq(e.y(), 2.f, "/ operator y");
    expect_eq(e.z(), 3.f, "/ operator z");

    auto f = basevector<float, 3>(2.f, 3.f, 4.f) * 2.5f;
    expect_eq(f.x(), 5.f, "* operator x");
    expect_eq(f.y(), 7.5f, "* operator y");
    expect_eq(f.z(), 10.f, "* operator z");

    auto g = -basevector<float, 3>(2.f, 3.f, 4.f);
    expect_eq(g.x(), -2.f, "negation operator x");
    expect_eq(g.y(), -3.f, "negation operator y");
    expect_eq(g.z(), -4.f, "negation operator z");
}

/**
 * @brief Test indexing behavior and bounds checking.
 *
 * Verifies that:
 *  - The non-const indexing operator allows assignment.
 *  - The const indexing operator throws `std::out_of_range` for indices
 *    >= component count (C). This checks defensive behavior for consumers
 *    of the const API.
 */
void test_index_bounds_and_access() {
    basevector<float, 3> v(1.f, 2.f, 3.f);
    // non-const operator() should work
    v(0) = 7.f;
    expect_eq(v.x(), 7.f, "operator() non-const set");

    // const version throws when out of range
    const basevector<float, 3> cv(1.f, 2.f, 3.f);
    bool threw = false;
    try {
        // attempt to access an out-of-range index
        volatile float val = cv(3);
        (void)val;
    }
    catch (const std::out_of_range&) {
        threw = true;
    }
    catch (...) {
        // other exceptions are failures
    }
    expect_true(threw, "const operator() throws out_of_range when pos >= C");
}

/**
 * @brief Test multiplication of a matrix by a vector.
 *
 * Constructs a 3x3 matrix equal to 2 * identity and multiplies it by a
 * vec3 to ensure each component is scaled by 2. This validates `basematrix`
 * indexing and the matrix-vector operator implementation.
 */
void test_matrix_vector_multiplication() {
    // Test m * v where m is 3x3 and v is vec3
    basematrix<float, 3, 3> m;
    // set matrix to identity * 2 for a predictable result
    for (size_t i = 0; i < 3 * 3; ++i) m[i] = 0.f;
    m[0] = 2.f; // row0 col0
    m[4] = 2.f; // row1 col1 (index 4 because row-major: 1*3+1)
    m[8] = 2.f; // row2 col2

    basevector<float, 3> v(1.f, 2.f, 3.f);
    auto r = m * v; // should scale each component by 2
    expect_eq(r.x(), 2.f, "matrix*vector result x");
    expect_eq(r.y(), 4.f, "matrix*vector result y");
    expect_eq(r.z(), 6.f, "matrix*vector result z");
}

/**
 * @brief Run all vector-related tests.
 *
 * Calls each individual test function and reports overall success or failure
 * using `g_failures`. Returns 0 when all tests pass, otherwise returns 1.
 *
 * @return int 0 on success, 1 when any tests failed.
 */
int run_all_vector_tests() {
    std::cout << "Running vector tests..." << std::endl;

    test_constructors_and_accessors();
    test_length_and_normalize();
    test_dot_and_cross();
    test_operators_and_arithmetic();
    test_index_bounds_and_access();
    test_matrix_vector_multiplication();

    if (g_failures == 0) {
        std::cout << "All vector tests PASSED." << std::endl;
        return 0;
    }
    else {
        std::cerr << g_failures << " test(s) FAILED." << std::endl;
        return 1;
    }
}
