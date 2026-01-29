#ifndef __geom_h__
#define __geom_h__

#include "base_definitions.h"
#include "vector.h"
#include "algebra.h"

namespace base_math {
    /**
     * @brief Calculates the unit normal vector of a plane defined by three points in counter-clockwise order.
     * @tparam T Numeric type.
     * @param vP1 First point.
     * @param vP2 Second point.
     * @param vP3 Third point.
     * @return Unit normal vector.
     */
    template <typename T>
    basevector<T, 3> calc_normal(const basevector<T, 3>& vP1, const basevector<T, 3>& vP2, const basevector<T, 3>& vP3) {
        basevector<T, 3> vNormal;
        basevector<T, 3> vV1(vP2 - vP1);
        basevector<T, 3> vV2(vP3 - vP1);
        vNormal = vV1.cross(vV2);
        vNormal.normalize();
        return vNormal;
    }

    /**
     * @brief Calculates the area of a triangle using Heron's formula.
     * @tparam T Numeric type.
     * @param e1 First edge vector.
     * @param e2 Second edge vector.
     * @param e3 Third edge vector.
     * @return Area of the triangle.
     */
    template <typename T>
    T triangle_area(const basevector<T, 3>& e1, const basevector<T, 3>& e2, const basevector<T, 3>& e3)
    {
        T l1 = e1.length();
        T l2 = e2.length();
        T l3 = e3.length();

        T s = 0.5f * (l1 + l2 + l3); // s ---> semiperimeter
        T v = s * (s - l1) * (s - l2) * (s - l3);
        T area = (T)sqrt(v);
        if (isnan(area))
            return 0;
        return area;
    }

    /**
     * @brief Computes the tangent of an angle, with a check for values near pi/2.
     * @tparam T Numeric type.
     * @param a Angle in radians.
     * @return Tangent value or 1000 if near pi/2.
     */
    template <typename T>
    T safe_tan(T a) {
        if (fabs(fabs(a) - HALF_PI<T>) > TOLLERANCE<T>)
            return (T)tan(a);
        return T(1000);
    }

    /**
     * @brief Calculates the angle between two vectors.
     * @tparam T Numeric type.
     * @param edge1 First vector.
     * @param edge2 Second vector.
     * @return Angle in radians.
     */
    template <typename T>
    T vector_angle(const basevector<T, 3>& edge1, const basevector<T, 3>& edge2)
    {
        T cos_of_angle = edge1.dot(edge2) * (1 / (edge1.length() * edge2.length()));
        T angle1 = T(acos(cos_of_angle));
        return angle1;
    }

    /**
 * @brief Computes the cotangent of the angle at point a formed by points b and c.
 *
 * This function interprets the input points as forming two vectors ba = b - a' and bc = c - a'
 * that share the common vertex at 'a'. It then computes the cotangent of the angle between
 * these two vectors using the relation:
 *
 *     cot(theta) = (ba * bc) / ||ba ? bc||
 *
 * A small cross-product length (below 'TOLLERANCE<T>') is treated as a degenerate case
 * and results in a return value of 'T(0)'.
 *
 * @tparam T Numeric type used for coordinates and computations.
 * @param a Vertex point at which the angle is measured.
 * @param b First point defining the first edge of the angle.
 * @param c Second point defining the second edge of the angle.
 * @return Cotangent of the angle at point 'a', or 'T(0)' for degenerate configurations.
 */
    template <typename T>
    T cotangent(const basevector<T, 3>& a, const basevector<T, 3>& b, const basevector<T, 3>& c) {
        basevector<T, 3> ba = b - a;
        basevector<T, 3> bc = c - a;
        basevector<T, 3> cross = ba.cross(bc);
        T cross_len = cross.length();
        if (cross_len < TOLLERANCE<T>)
            return T(0); // degenerate case
        return ba.dot(bc) / cross_len;
    }

    /**
     * @brief Computes the interior angle at vertex pa for the triangle (pa, pb, pc).
     *
     * The function constructs two edge vectors 'u = pb - pa' and 'v = pc - pa', both
     * originating from the vertex 'pa'. It then computes the angle between them using
     * the standard dot-product relation:
     *
     *     cos(theta) = (u * v) / (||u|| * ||v||)
     *
     * The cosine value is clamped to the range [-1, 1] to improve numerical robustness
     * before applying 'std::acos'. If either edge vector has zero length, the
     * configuration is treated as degenerate and 'T(0)' is returned.
     *
     * @tparam T Numeric type used for coordinates and computations.
     * @param pa Vertex point at which the corner angle is measured.
     * @param pb Second point defining the first incident edge.
     * @param pc Third point defining the second incident edge.
     * @return Corner angle at 'pa' in radians, or 'T(0)' for degenerate configurations.
     */
    template <typename T>
    T corner_angle(const basevector<T, 3>& pa, const basevector<T, 3>& pb, const basevector<T, 3>& pc) {
        basevector<T, 3> u = pb - pa;
        basevector<T, 3> v = pc - pa;

        T lu = u.length();
        T lv = v.length();

        if (lu == T(0) || lv == T(0))
            return T(0); // degenerate case

        T cosang = u.dot(v) / (lu * lv);

        // Clamp for numerical safety
        cosang = T(std::max<T>(-1.0, std::min<T>(1.0, cosang)));

        return T(std::acos(cosang));
    }

    /**
     * @brief Tests whether two pairs of vectors define coplanar directions.
     *
     * Given two vector pairs '(u_f, v_f)' and '(u_p, v_p)', this function checks if the
     * planes spanned by each pair are the same (or nearly the same) by comparing their
     * normals. It computes:
     *
     *   - 'cross1 = u_f * v_f' : normal to the first plane
     *   - 'cross2 = u_p * v_p' : normal to the second plane
     *   - 'cross_res = cross1 * cross2'
     *
     * If 'cross_res' has near-zero length (less than 'TOLLERANCE<T>'), the two normals
     * are considered parallel or anti-parallel, and the vector pairs are treated as
     * coplanar.
     *
     * @tparam T Numeric type used for coordinates and computations.
     * @param u_f First vector of the reference pair.
     * @param v_f Second vector of the reference pair.
     * @param u_p First vector of the pair to be tested.
     * @param v_p Second vector of the pair to be tested.
     * @return 'true' if the vector pairs are coplanar within tolerance, otherwise 'false'.
     */
    template <typename T>
    bool isCoplanar(const basevector<T, 3>& u_f, const basevector<T, 3>& v_f, const basevector<T, 3>& u_p, const basevector<T, 3>& v_p)
    {
        basevector<T, 3> cross1 = u_f.cross(v_f);
        basevector<T, 3> cross2 = u_p.cross(v_p);
        basevector<T, 3> cross_res = cross1.cross(cross2);
        return cross_res.length() < TOLLERANCE<T>;
    }

    /**
     * @brief Constructs a local orthonormal tangent frame (u, v) for a given normal.
     *
     * This function builds two tangent vectors `u` and `v` that are:
     *   - orthogonal to the input `normal` (i.e., lie in the plane perpendicular to it),
     *   - mutually orthogonal to each other (forming a right-handed local frame when
     *     combined with `normal`),
     *   - chosen in a numerically stable way by handling near-axis-aligned normals
     *     as special cases.
     *
     * The resulting frame can be used as a 2D parameterization basis (local UV
     * coordinates) on a surface for tasks such as curvature computation, texture
     * mapping, or local differential-geometry operations.
     *
     * For robustness:
     *   - The function first checks that the input `normal` is non-degenerate
     *     (its length must be at least `TOLLERANCE<float>`); otherwise it returns `false`.
     *   - If the normal is nearly aligned with the major axes (Z or Y), fixed
     *     axis-aligned vectors are used for `u` and `v` to avoid numerical issues.
     *   - In the general case, `u` and `v` are constructed algebraically from the
     *     components of `normal` and then normalized.
     *
     * @tparam T Numeric type used for coordinates and computations.
     * @param normal Input normal vector that defines the plane's orientation.
     * @param u Output tangent vector lying in the plane perpendicular to `normal`.
     * @param v Output tangent vector lying in the plane perpendicular to `normal`
     *          and orthogonal to `u`.
     * @return `true` if a valid tangent frame is constructed, `false` if the input
     *         `normal` is too small (degenerate).
     */
    template <typename T>
    bool create_uv_reference_plane(const basevector<T, 3>& normal, basevector<T, 3>& u, basevector<T, 3>& v)
    {
        // Check for degenerate normal
        if (normal.length() < TOLLERANCE<float>) return false;

        // Handle special cases where normal is aligned with major axes
        if (fabs(normal.z()) > 0.95f) {
            u = basevector<T, 3>(1, 0, 0);
            v = basevector<T, 3>(0, 1, 0);
        }
        else if (fabs(normal.y()) > 0.95f) {
            u = basevector<T, 3>(1, 0, 0);
            v = basevector<T, 3>(0, 0, 1);
        }
        // General case
        else {
            u = basevector<T, 3>(normal.y(), -normal.x(), 0);
            v = basevector<T, 3>(normal.x() * normal.z(), normal.y() * normal.z(), -normal.x() * normal.x() - normal.y() * normal.y());
            u.normalize();
            v.normalize();
        }
        return true;
    }

}

#endif // __geom_h__
