#include "gl_math.h"

namespace base_opengl {
    //////////////////////////////////////////////////////////////////////////
    // view matrices, column major as OpenGL expects them

    /**
     * @brief Constructs a view matrix for a camera in 3D space.
     *
     * This function creates a column-major view matrix as expected by OpenGL, positioning
     * the camera at `eyePosition3D`, looking towards `center3D`, with the given `upVector3D`.
     *
     * @param eyePosition3D The position of the camera (eye) in world coordinates.
     * @param center3D The point in world space the camera is looking at.
     * @param upVector3D The up direction vector.
     * @return fmat4 The resulting view matrix.
     */
    fmat4 lookAt(const fvec3& eyePosition3D, const fvec3& center3D, const fvec3& upVector3D) {
        fvec3 forward(center3D - eyePosition3D);
        forward.normalize();
        // --------------------
        // right = forward x up
        fvec3 right = forward.cross(upVector3D);
        right.normalize();
        // --------------------
        // Recompute up as: up = right x forward
        fvec3 up = right.cross(forward);
        up.normalize();
        // --------------------
        fmat4 result;
        result[0] = right[0];
        result[4] = right[1];
        result[8] = right[2];
        result[12] = -right.dot(eyePosition3D);
        // --------------------
        result[1] = up[0];
        result[5] = up[1];
        result[9] = up[2];
        result[13] = -up.dot(eyePosition3D);
        // --------------------
        result[2] = -forward[0];
        result[6] = -forward[1];
        result[10] = -forward[2];
        result[14] = forward.dot(eyePosition3D);
        // --------------------
        result[3] = result[7] = result[11] = 0.0;
        result[15] = 1.0;
        // --------------------
        return result;
    }

    /**
     * @brief Constructs a perspective frustum projection matrix.
     *
     * Fills the provided 4x4 matrix (column-major) with a frustum projection, as used in OpenGL.
     *
     * @param matrix Pointer to a 16-element float array to receive the matrix.
     * @param left Left vertical clipping plane.
     * @param right Right vertical clipping plane.
     * @param bottom Bottom horizontal clipping plane.
     * @param top Top horizontal clipping plane.
     * @param znear Near depth clipping plane (must be positive).
     * @param zfar Far depth clipping plane (must be positive).
     */
    void frustum_matrix(float* matrix, float left, float right, float bottom, float top, float znear, float zfar) {
        float n2, rml, tmb, fmn;
        n2 = 2.0f * znear;
        rml = right - left;
        tmb = top - bottom;
        fmn = zfar - znear;

        memset(matrix, 0, 16 * sizeof(float));
        matrix[0] = n2 / rml;
        matrix[5] = n2 / tmb;
        matrix[8] = (right + left) / rml;
        matrix[9] = (top + bottom) / tmb;
        matrix[10] = -zfar / fmn;
        matrix[11] = -1.0;
        matrix[14] = -(znear * zfar) / fmn;
    }

    /**
     * @brief Constructs a perspective projection matrix.
     *
     * Fills the provided 4x4 matrix (column-major) with a perspective projection matrix
     * based on the specified field of view, aspect ratio, and near/far planes.
     *
     * @param matrix Pointer to a 16-element float array to receive the matrix.
     * @param fovyInRadians Field of view angle in the y direction, in radians.
     * @param aspectRatio Aspect ratio (width/height).
     * @param znear Near depth clipping plane (must be positive).
     * @param zfar Far depth clipping plane (must be positive).
     */
    void perspective_matrix(float* matrix, float fovyInRadians, float aspectRatio, float znear, float zfar) {
        float ymax, xmax;
        // half of the angle!!
        ymax = znear * (float)tan(fovyInRadians / 2.f);
        xmax = ymax * aspectRatio;
        frustum_matrix(matrix, -xmax, xmax, -ymax, ymax, znear, zfar);
    }

    /**
     * @brief Constructs and returns a perspective projection matrix.
     *
     * Returns a column-major 4x4 perspective projection matrix based on the specified
     * field of view, aspect ratio, and near/far planes.
     *
     * @param fovyInRadians Field of view angle in the y direction, in radians.
     * @param aspectRatio Aspect ratio (width/height).
     * @param znear Near depth clipping plane (must be positive).
     * @param zfar Far depth clipping plane (must be positive).
     * @return fmat4 The resulting perspective projection matrix.
     */
    fmat4 perspective_matrix(float fovyInRadians, float aspectRatio, float znear, float zfar) {
        fmat4 ret;
        perspective_matrix((float*)ret, fovyInRadians, aspectRatio, znear, zfar);
        return ret;
    }

    /**
     * @brief Constructs an orthographic projection matrix.
     *
     * Fills the provided 4x4 matrix (column-major) with an orthographic projection matrix
     * as used in OpenGL.
     *
     * @param matrix Pointer to a 16-element float array to receive the matrix.
     * @param l Left vertical clipping plane.
     * @param r Right vertical clipping plane.
     * @param b Bottom horizontal clipping plane.
     * @param t Top horizontal clipping plane.
     * @param n Near depth clipping plane.
     * @param f Far depth clipping plane.
     */
    void ortho(float* matrix, float l, float r, float b, float t, float n, float f)
    {
        // set OpenGL orthographic projection matrix
        matrix[0] = 2 / (r - l);
        matrix[1] = 0;
        matrix[2] = 0;
        matrix[3] = 0;

        matrix[4] = 0;
        matrix[5] = 2 / (t - b);
        matrix[6] = 0;
        matrix[7] = 0;

        matrix[8] = 0;
        matrix[9] = 0;
        matrix[10] = -2 / (f - n);
        matrix[11] = 0;

        matrix[12] = -(r + l) / (r - l);
        matrix[13] = -(t + b) / (t - b);
        matrix[14] = -(f + n) / (f - n);
        matrix[15] = 1;
    }

    /**
     * @brief Constructs and returns an orthographic projection matrix.
     *
     * Returns a column-major 4x4 orthographic projection matrix as used in OpenGL.
     *
     * @param l Left vertical clipping plane.
     * @param r Right vertical clipping plane.
     * @param b Bottom horizontal clipping plane.
     * @param t Top horizontal clipping plane.
     * @param n Near depth clipping plane.
     * @param f Far depth clipping plane.
     * @return fmat4 The resulting orthographic projection matrix.
     */
    fmat4 ortho(float l, float r, float b, float t, float n, float f)
    {
        fmat4 o;
        ortho(o, l, r, b, t, n, f);
        return o;
    }
}  // namespace base_opengl
