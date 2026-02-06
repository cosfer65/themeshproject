#ifndef __quaternion_h__
#define __quaternion_h__

#include "base_definitions.h"
#include "vector.h"

namespace base_math {
    /**
     * @brief Quaternion class for representing rotations and orientations.
     *
     * Inherits from basevector<T, 4> and provides quaternion operations.
     * The quaternion is stored as (x, y, z, w), where w is the scalar part.
     *
     * @tparam T Numeric type (float, double, etc.)
     */
    template <typename T>
    class quaternion : public basevector<T, 4> {
    public:
        /**
         * @brief Default constructor. Initializes all components to zero.
         */
        quaternion() : basevector<T, 4>(0) {
        }

        /**
         * @brief Constructs a quaternion from scalar and vector components.
         * @param scalar Scalar (w) part.
         * @param i X (i) component.
         * @param j Y (j) component.
         * @param k Z (k) component.
         */
        quaternion(T scalar, T i, T j, T k) : basevector<T, 4>(i, j, k, scalar) {
        }

        /**
         * @brief Returns the magnitude (length) of the quaternion.
         * @return Magnitude as T.
         */
        T magnitude(void) {
            return this->length();
        }

        /**
         * @brief Gets the vector (x, y, z) part of the quaternion.
         * @return basevector<T, 3> containing (x, y, z).
         */
        basevector<T, 3> get_vector(void) const {
            return basevector<T, 3>(this->x(), this->y(), this->z());
        }

        /**
         * @brief Gets the scalar (w) part of the quaternion.
         * @return Scalar value as T.
         */
        T get_scalar(void) const {
            return this->w();
        }

        /**
         * @brief Returns the conjugate of the quaternion.
         * @return Conjugated quaternion.
         */
        quaternion operator~(void) const {
            return quaternion(this->w(), -this->x(), -this->y(), -this->z());
        }

        /**
         * @brief Converts the quaternion to a 4x4 rotation matrix.
         *
         * Computes a homogeneous 4x4 rotation matrix from this quaternion using the
         * standard quaternion-to-matrix conversion formula. The resulting matrix
         * represents the same rotation as the quaternion and can be used for 3D
         * transformations in graphics pipelines.
         *
         * The matrix is stored in column-major order and follows the form:
         *
         * [ 1-2(y^2+z^2)   2(xy-wz)       2(xz+wy)       0 ]
         * [ 2(xy+wz)       1-2(x^2+z^2)   2(yz-wx)       0 ]
         * [ 2(xz-wy)       2(yz+wx)       1-2(x^2+y^2)   0 ]
         * [ 0              0              0              1 ]
         *
         * where x, y, z are the vector components and w is the scalar component
         * of the quaternion.
         *
         * @return basematrix<T, 4, 4> A 4x4 homogeneous rotation matrix with no
         *         translation component (translation elements are set to 0).
         *
         * @note The quaternion does not need to be normalized before calling this
         *       method, but non-normalized quaternions will produce scaled rotation
         *       matrices which may not preserve distances.
         */
        void matrix(T M[16]) const {
            T x = this->x(), y = this->y(), z = this->z(), w = this->w();
            T xx = x * x, yy = y * y, zz = z * z;
            T xy = x * y, xz = x * z, yz = y * z;
            T wx = w * x, wy = w * y, wz = w * z;

            M[0] = 1 - 2 * (yy + zz);
            M[1] = 2 * (xy + wz);
            M[2] = 2 * (xz - wy);
            M[3] = 0;

            M[4] = 2 * (xy - wz);
            M[5] = 1 - 2 * (xx + zz);
            M[6] = 2 * (yz + wx);
            M[7] = 0;

            M[8] = 2 * (xz + wy);
            M[9] = 2 * (yz - wx);
            M[10] = 1 - 2 * (xx + yy);
            M[11] = 0;

            M[12] = M[13] = M[14] = 0;
            M[15] = 1;
        }

        /**
         * @brief Convinience function returning a basematrix<T, 4, 4>
         * made from the matrix calculated in the function above.
         * @return basematrix<T, 4, 4> representing the rotation.
         */
        basematrix<T, 4, 4> matrix() const {
            T M[16];
            matrix(M);
            return basematrix<T, 4, 4>(M);
        }

        /**
         * @brief Multiplies two quaternions.
         * @param q2 Second quaternion.
         * @return Product quaternion.
         */
        template <typename T>
        quaternion<T> operator * (const quaternion<T>& q2) {
            return quaternion(
                this->w() * q2.w() - this->x() * q2.x() - this->y() * q2.y() - this->z() * q2.z(),
                this->w() * q2.x() + this->x() * q2.w() + this->y() * q2.z() - this->z() * q2.y(),
                this->w() * q2.y() - this->x() * q2.z() + this->y() * q2.w() + this->z() * q2.x(),
                this->w() * q2.z() + this->x() * q2.y() - this->y() * q2.x() + this->z() * q2.w());
        }
    };

    /**
     * @brief Constructs a unit quaternion from an axis-angle representation.
     *
     * Creates a quaternion that represents a rotation of `angle` radians
     * around the given 3D `axis`. The axis is interpreted as the rotation
     * axis in 3D space and is expected to be normalized (unit length)
     * for the resulting quaternion to represent a pure rotation without
     * additional scaling.
     *
     * The quaternion components are computed as:
     *   - w = cos(angle / 2)
     *   - (x, y, z) = axis * sin(angle / 2)
     *
     * @tparam T Numeric scalar type (e.g., `float`, `double`).
     * @param axis Rotation axis as a 3D vector. Should be normalized for
     *             correct rotational behavior.
     * @param angle Rotation angle in radians, measured around `axis`
     *              using the right-hand rule.
     * @return Quaternion representing the specified axis-angle rotation.
     *
     * @note If `axis` is not normalized, the resulting quaternion will be
     *       scaled, and any rotations derived from it may introduce
     *       non-uniform scaling artifacts.
     */
    template <typename T>
    inline quaternion<T> fromAxisAngle(const basevector<T, 3>& axis, T angle) {
        T s = T(std::sin(angle * T(0.5)));
        return quaternion<T>(T(std::cos(angle * T(0.5))), axis.x() * s, axis.y() * s, axis.z() * s);
    }


    /**
     * @brief Multiplies a 3D vector by a quaternion.
     * @param v The vector.
     * @param q The quaternion.
     * @return Resulting quaternion.
     */
    template <typename T>
    quaternion<T> operator * (const basevector<T, 3>& v, const quaternion<T>& q) {
        return quaternion(-(q.x() * v.x() + q.y() * v.y() + q.z() * v.z()),
            q.w() * v.x() + q.z() * v.y() - q.y() * v.z(),
            q.w() * v.y() + q.x() * v.z() - q.z() * v.x(),
            q.w() * v.z() + q.y() * v.x() - q.x() * v.y());
    }

    /**
     * @brief Multiplies a quaternion by a 3D vector.
     * @param q The quaternion.
     * @param v The vector.
     * @return Resulting quaternion.
     */
    template <typename T>
    quaternion<T> operator * (const quaternion<T>& q, const basevector<T, 3>& v) {
        return    quaternion(-(q.x() * v.x() + q.y() * v.y() + q.z() * v.z()),
            q.w() * v.x() + q.y() * v.z() - q.z() * v.y(),
            q.w() * v.y() + q.z() * v.x() - q.x() * v.z(),
            q.w() * v.z() + q.x() * v.y() - q.y() * v.x());
    }

    /**
     * @brief Gets the rotation angle represented by a quaternion.
     * @param q The quaternion.
     * @return Angle in radians.
     */
    template <typename T>
    T q_get_angle(const quaternion<T>& q) {
        return (T)(2 * acos(q.w()));
    }

    /**
     * @brief Gets the rotation axis from a quaternion.
     * @param q The quaternion.
     * @return Normalized axis as basevector<T, 3>.
     */
    template <typename T>
    basevector<T, 3> q_get_axis(const quaternion<T>& q) {
        basevector<T, 3> v;
        T m;

        v = q.get_vector();
        m = v.length();

        if (m <= TOLLERANCE<T>)
            return basevector<T, 3>();
        else
            return v / m;
    }

    /**
     * @brief Rotates quaternion q2 by quaternion q1.
     * @param q1 Rotation quaternion.
     * @param q2 Quaternion to rotate.
     * @return Rotated quaternion.
     */
    template <typename T>
    quaternion<T> q_rotate(const quaternion<T>& q1, const quaternion<T>& q2) {
        return q1 * q2 * (~q1);
    }

    /**
     * @brief Rotates a vector by the inverse of a quaternion.
     * @param q The quaternion.
     * @param v The vector.
     * @return Rotated vector.
     */
    template <typename T>
    basevector<T, 3> qv_rotate_inv(const quaternion<T>& q, const basevector<T, 3>& v) {
        return ((~q) * v * q).get_vector();
    }

    /**
     * @brief Rotates a vector by a quaternion.
     * @param q The quaternion.
     * @param v The vector.
     * @return Rotated vector.
     */
    template <typename T>
    basevector<T, 3> qv_rotate(const quaternion<T>& q, const basevector<T, 3>& v) {
        return (q * v * (~q)).get_vector();
    }

    /**
     * @brief Creates a quaternion from Euler angles (in degrees).
     * @param x Pitch angle.
     * @param y Yaw angle.
     * @param z Roll angle.
     * @return Quaternion representing the rotation.
     */
    template <typename T>
    quaternion<T> make_q_from_euler_angles(T x, T y, T z) {
        double roll = deg_to_rad(z);
        double pitch = deg_to_rad(x);
        double yaw = deg_to_rad(y);

        double cyaw, cpitch, croll, syaw, spitch, sroll;
        double cyawcpitch, syawspitch, cyawspitch, syawcpitch;

        cyaw = cos(0.5f * yaw);
        cpitch = cos(0.5f * pitch);
        croll = cos(0.5f * roll);
        syaw = sin(0.5f * yaw);
        spitch = sin(0.5f * pitch);
        sroll = sin(0.5f * roll);

        cyawcpitch = cyaw * cpitch;
        syawspitch = syaw * spitch;
        cyawspitch = cyaw * spitch;
        syawcpitch = syaw * cpitch;

        quaternion<T> q;
        q.w() = (T)(cyawcpitch * croll + syawspitch * sroll);
        q.x() = (T)(cyawcpitch * sroll - syawspitch * croll);
        q.y() = (T)(cyawspitch * croll + syawcpitch * sroll);
        q.z() = (T)(syawcpitch * croll - cyawspitch * sroll);
        q.normalize();

        return q;
    }

    /**
     * @brief Converts a quaternion to Euler angles (in radians).
     * @param q The quaternion.
     * @return basevector<T, 3> containing (roll, pitch, yaw).
     */
    template <typename T>
    basevector<T, 3> make_euler_angles_from_q(const quaternion<T>& q) {
        double r11, r21, r31, r32, r33, r12, r13;
        double q00, q11, q22, q33;
        double tmp;
        basevector<T, 3> u;

        q00 = q.w() * q.w();
        q11 = q.x() * q.x();
        q22 = q.y() * q.y();
        q33 = q.z() * q.z();

        r11 = q00 + q11 - q22 - q33;
        r21 = 2 * (q.x() * q.y() + q.w() * q.z());
        r31 = 2 * (q.x() * q.z() - q.w() * q.y());
        r32 = 2 * (q.y() * q.z() + q.w() * q.x());
        r33 = q00 - q11 - q22 + q33;

        tmp = fabs(r31);
        if (tmp > 0.999999) {
            r12 = 2 * (q.x() * q.y() - q.w() * q.z());
            r13 = 2 * (q.x() * q.z() + q.w() * q.y());

            u.x() = 0.f; //roll
            u.y() = (T)(-(PI<T> / 2) * r31 / tmp); // pitch
            u.z() = (T)atan2(-r12, -r31 * r13); // yaw
            return u;
        }

        u.x() = (T)atan2(r32, r33); // roll
        u.y() = (T)asin(-r31); // pitch
        u.z() = (T)atan2(r21, r11); // yaw
        return u;
    }

    /**
     * @brief Rotates a vector around a given axis by an angle.
     * @param v Vector to rotate (modified in place).
     * @param rotation_axis Axis to rotate around.
     * @param angle Angle in radians.
     */
    template <typename T>
    void rotate_vector(basevector<T, 3>& v, const basevector<T, 3>& rotation_axis, T angle) {
        // build the rotation quaternion
        quaternion<T> q;
        q.w() = (T)cos(angle / 2);
        q.x() = (T)sin(angle / 2) * rotation_axis.x();
        q.y() = (T)sin(angle / 2) * rotation_axis.y();
        q.z() = (T)sin(angle / 2) * rotation_axis.z();
        // rotate the vector
        v = qv_rotate(q, v);
    }

    /**
     * @brief Calculates a transformation matrix from position and orientation.
     * @param transformMatrix Output 4x4 matrix.
     * @param position Translation vector.
     * @param orientation Rotation quaternion.
     */
    template <typename T>
    void calculateTransformationMatrix(basematrix<T, 4, 4>& transformMatrix, const basevector<T, 3>& position, const quaternion<T>& orientation)
    {
        transformMatrix[0] = 1 - 2 * orientation.y() * orientation.y() -
            2 * orientation.z() * orientation.z();
        transformMatrix[1] = 2 * orientation.x() * orientation.y() -
            2 * orientation.w() * orientation.z();
        transformMatrix[2] = 2 * orientation.x() * orientation.z() +
            2 * orientation.w() * orientation.y();
        transformMatrix[3] = position.x();

        transformMatrix[4] = 2 * orientation.x() * orientation.y() +
            2 * orientation.w() * orientation.z();
        transformMatrix[5] = 1 - 2 * orientation.x() * orientation.x() -
            2 * orientation.z() * orientation.z();
        transformMatrix[6] = 2 * orientation.y() * orientation.z() -
            2 * orientation.w() * orientation.x();
        transformMatrix[7] = position.y();

        transformMatrix[8] = 2 * orientation.x() * orientation.z() -
            2 * orientation.w() * orientation.y();
        transformMatrix[9] = 2 * orientation.y() * orientation.z() +
            2 * orientation.w() * orientation.x();
        transformMatrix[10] = 1 - 2 * orientation.x() * orientation.x() -
            2 * orientation.y() * orientation.y();
        transformMatrix[11] = position.z();
    }

    /**
     * @brief Extracts Euler angles from a 3x3 rotation matrix.
     *
     * Interprets the input array `R` as a 3x3 rotation matrix stored in
     * row-major order:
     *
     *   R = [ r00 r01 r02
     *         r10 r11 r12
     *         r20 r21 r22 ]
     *
     * The function assumes a rotation order consistent with the rest of this
     * module (yaw–pitch–roll convention) and returns the three Euler angles
     * (stored in a `basevector<T, 3>`) in radians:
     *
     *   - `e.x()` : yaw angle
     *   - `e.y()` : pitch angle
     *   - `e.z()` : roll angle
     *
     * It computes the pitch from element r02 and then derives yaw and roll
     * from the remaining elements. When the cosine of the pitch angle is
     * close to zero (i.e., pitch is near ±90 degrees), the function detects
     * gimbal lock and switches to an alternative computation path to avoid
     * numerical instability. In this singular case, roll is set to zero and
     * yaw is derived from r10 and r11.
     *
     * @tparam T Numeric scalar type (e.g., `float`, `double`).
     * @param R Pointer to an array of at least 9 elements representing a 3x3
     *          rotation matrix in row-major layout.
     * @return `basevector<T, 3>` containing the Euler angles in radians:
     *         `(yaw, pitch, roll)`.
     *
     * @note The input matrix is assumed to be a valid rotation matrix
     *       (orthonormal with determinant +1). Invalid or noisy matrices may
     *       produce undefined or unstable angle outputs.
     */
    template <typename T>
    inline basevector<T,3> eulerAnglesFromRotationMatrix(const T* R)
    {
        basevector<T, 3> e;

        // Extract matrix elements for readability
        const T r00 = R[0], r01 = R[1], r02 = R[2];
        const T r10 = R[3], r11 = R[4], r12 = R[5];
        const T r20 = R[6], r21 = R[7], r22 = R[8];

        // Pitch (around Y)
        e.y() = T(std::asin(r02));

        // Check for gimbal lock
        const T cos_pitch = T(std::cos(e.y()));

        if (std::fabs(cos_pitch) > 1e-6)
        {
            // Standard case
            e.x() = T(std::atan2(-r01, r00));  // yaw
            e.z() = T(std::atan2(-r12, r22));  // roll
        }
        else
        {
            // Gimbal lock: pitch is ±90 degrees
            e.x() = T(std::atan2(r10, r11));
            e.z() = T(0.0);
        }

        return e;
    }


    /**
     * @brief Typedef for float quaternion.
     */
    typedef quaternion<float> fquaternion;
    /**
     * @brief Typedef for double quaternion.
     */
    typedef quaternion<double> dquaternion;
}
#endif // __quaternion_h__
