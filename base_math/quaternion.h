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
         * @return basematrix<T, 4, 4> representing the rotation.
         */
        basematrix<T, 4, 4> matrix() const {
            basematrix<T, 4, 4> R;
            R[0] = 1 - 2 * this->y() * this->y() - 2 * this->z() * this->z();
            R[1] = 2 * this->x() * this->y() - 2 * this->w() * this->z();
            R[2] = 2 * this->x() * this->z() + 2 * this->w() * this->y();
            R[3] = 0;

            R[4] = 2 * this->x() * this->y() + 2 * this->w() * this->z();
            R[5] = 1 - 2 * this->x() * this->x() - 2 * this->z() * this->z();
            R[6] = 2 * this->y() * this->z() - 2 * this->w() * this->x();
            R[7] = 0;

            R[8] = 2 * this->x() * this->z() - 2 * this->w() * this->y();
            R[9] = 2 * this->y() * this->z() + 2 * this->w() * this->x();
            R[10] = 1 - 2 * this->x() * this->x() - 2 * this->y() * this->y();
            R[11] = 0;

            R[12] = R[13] = R[14] = 0;
            R[15] = 1;
            return R;
        }
    };

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
     * @brief Multiplies two quaternions.
     * @param q1 First quaternion.
     * @param q2 Second quaternion.
     * @return Product quaternion.
     */
    template <typename T>
    quaternion<T> operator * (const quaternion<T>& q1, const quaternion<T>& q2) {
        return quaternion(q1.w() * q2.w() - q1.x() * q2.x() - q1.y() * q2.y() - q1.z() * q2.z(),
            q1.w() * q2.x() + q1.x() * q2.w() + q1.y() * q2.z() - q1.z() * q2.y(),
            q1.w() * q2.y() + q1.y() * q2.w() + q1.z() * q2.x() - q1.x() * q2.z(),
            q1.w() * q2.z() + q1.z() * q2.w() + q1.x() * q2.y() - q1.y() * q2.x());
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
     * @brief Typedef for float quaternion.
     */
    typedef quaternion<float> fquaternion;
    /**
     * @brief Typedef for double quaternion.
     */
    typedef quaternion<double> dquaternion;
}
#endif // __quaternion_h__
