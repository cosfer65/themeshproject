#include "gl_graphics.h"
#include "gl_camera.h"

namespace base_opengl {

    /**
     * @brief Rotates the camera around its local X axis.
     *
     * This function rotates both the forward and up vectors of the camera
     * around the axis perpendicular to the current forward and up vectors
     * (the local X axis), by the specified angle in radians.
     *
     * @param fAngle The angle in radians to rotate around the local X axis.
     */
    void gl_camera::rotate_locX(float fAngle)
    {
        fmat4 mRotation;
        fvec3 vCross;

        vCross = vForward.cross(vUp);
        mRotation = rotation_matrix(fAngle, vCross.x(), vCross.y(), vCross.z());

        fvec3 vNewVect;
        // Inline 3x3 matrix multiply for rotation only
        vNewVect.x() = mRotation[0] * vForward.x() + mRotation[4] * vForward.y() + mRotation[8] * vForward.z();
        vNewVect.y() = mRotation[1] * vForward.x() + mRotation[5] * vForward.y() + mRotation[9] * vForward.z();
        vNewVect.z() = mRotation[2] * vForward.x() + mRotation[6] * vForward.y() + mRotation[10] * vForward.z();
        // memcpy(vForward, vNewVect, sizeof(float) * 3);
        vForward = vNewVect;

        // Update pointing up vector
        vNewVect.x() = mRotation[0] * vUp.x() + mRotation[4] * vUp.y() + mRotation[8] * vUp.z();
        vNewVect.y() = mRotation[1] * vUp.x() + mRotation[5] * vUp.y() + mRotation[9] * vUp.z();
        vNewVect.z() = mRotation[2] * vUp.x() + mRotation[6] * vUp.y() + mRotation[10] * vUp.z();
        // memcpy(vUp, vNewVect, sizeof(float) * 3);
        vUp = vNewVect;
    }

    /**
     * @brief Rotates the camera around its local Y axis.
     *
     * This function rotates the forward vector of the camera around the up vector
     * (the local Y axis), by the specified angle in radians.
     *
     * @param angle The angle in radians to rotate around the local Y axis.
     */
    void gl_camera::rotate_locY(float angle)
    {
        fmat4 mRotation;

        mRotation = rotation_matrix(angle, vUp.x(), vUp.y(), vUp.z());

        vForward = rotate_vector(mRotation, vForward);
    }

    /**
     * @brief Rotates the camera around its local Z axis.
     *
     * This function rotates the up vector of the camera around the forward vector
     * (the local Z axis), by the specified angle in radians.
     *
     * @param fAngle The angle in radians to rotate around the local Z axis.
     */
    void gl_camera::rotate_locZ(float fAngle)
    {
        fmat4 mRotation;

        // Only the up vector needs to be rotated
        mRotation = rotation_matrix(-fAngle, vForward.x(), vForward.y(), vForward.z());

        fvec3 vNewVect;
        vNewVect.x() = mRotation[0] * vUp.x() + mRotation[4] * vUp.y() + mRotation[8] * vUp.z();
        vNewVect.y() = mRotation[1] * vUp.x() + mRotation[5] * vUp.y() + mRotation[9] * vUp.z();
        vNewVect.z() = mRotation[2] * vUp.x() + mRotation[6] * vUp.y() + mRotation[10] * vUp.z();
        vUp = vNewVect;
    }

    /**
     * @brief Sets up the camera with a position, look-at point, and up vector.
     *
     * Initializes the camera's position, forward, and up vectors. The forward vector
     * is calculated as the normalized direction from the position to the look-at point.
     * The up vector is re-orthogonalized to ensure it is perpendicular to the forward vector.
     *
     * @param pos The camera position.
     * @param lookat The point the camera is looking at.
     * @param up The up direction vector.
     */
    void gl_camera::setup(const fvec3& pos, const fvec3& lookat, const fvec3& up)
    {
        vLocation = pos;
        vForward = lookat - pos;
        vUp = up;
        vForward.normalize();

        fvec3 right = vForward.cross(vUp);
        vUp = right.cross(vForward);
        vUp.normalize();
    }

    /**
     * @brief Re-orthogonalizes the camera's up vector.
     *
     * Normalizes the forward vector and recalculates the up vector to ensure
     * it is perpendicular to the forward vector.
     */
    void gl_camera::setup_d() {
        vForward.normalize();
        fvec3 right = vForward.cross(vUp);
        vUp = right.cross(vForward);
        vUp.normalize();
    }

    /**
     * @brief Sets up the camera with explicit position, forward, and up vectors.
     *
     * Initializes the camera's position, forward, and up vectors. The forward vector
     * is normalized, and the up vector is re-orthogonalized to ensure it is perpendicular
     * to the forward vector.
     *
     * @param pos The camera position.
     * @param forw The forward direction vector.
     * @param up The up direction vector.
     */
    void gl_camera::setup_d(const fvec3& pos, const fvec3& forw, const fvec3& up) {
        vLocation = pos;
        vForward = forw;
        vUp = up;
        vForward.normalize();

        fvec3 right = vForward.cross(vUp);
        vUp = right.cross(vForward);
        vUp.normalize();
    }
}