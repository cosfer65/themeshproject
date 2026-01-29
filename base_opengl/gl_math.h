#ifndef __gl_math_h__
#define __gl_math_h__

#include "vector.h"

using namespace base_math;

namespace base_opengl {
    //////////////////////////////////////////////////////////////////////////
    /**
     * @brief Constructs a view matrix for a camera in 3D space.
     * 
     * The resulting matrix is column-major, as expected by OpenGL.
     * @param eyePosition3D The position of the camera (eye) in world coordinates.
     * @param center3D The point in world space the camera is looking at.
     * @param upVector3D The up direction for the camera.
     * @return fmat4 The resulting view matrix.
     */
    fmat4 lookAt(const fvec3& eyePosition3D, const fvec3& center3D, const fvec3& upVector3D);

    /**
     * @brief Creates a frustum (perspective) projection matrix.
     * 
     * The matrix is stored in column-major order.
     * @param matrix Pointer to a 16-element array to receive the matrix.
     * @param left Left vertical clipping plane.
     * @param right Right vertical clipping plane.
     * @param bottom Bottom horizontal clipping plane.
     * @param top Top horizontal clipping plane.
     * @param znear Distance to the near depth clipping plane (must be positive).
     * @param zfar Distance to the far depth clipping plane (must be positive).
     */
    void frustum_matrix(float* matrix, float left, float right, float bottom, float top, float znear, float zfar);

    /**
     * @brief Creates a perspective projection matrix using field of view.
     * 
     * The matrix is stored in column-major order.
     * @param matrix Pointer to a 16-element array to receive the matrix.
     * @param fovyInRadians Field of view angle in the y direction, in radians.
     * @param aspectRatio Aspect ratio that determines the field of view in the x direction.
     * @param znear Distance to the near depth clipping plane (must be positive).
     * @param zfar Distance to the far depth clipping plane (must be positive).
     */
    void perspective_matrix(float* matrix, float fovyInRadians, float aspectRatio, float znear, float zfar);

    /**
     * @brief Constructs and returns a perspective projection matrix.
     * 
     * @param fovyInRadians Field of view angle in the y direction, in radians.
     * @param aspectRatio Aspect ratio that determines the field of view in the x direction.
     * @param znear Distance to the near depth clipping plane (must be positive).
     * @param zfar Distance to the far depth clipping plane (must be positive).
     * @return fmat4 The resulting perspective projection matrix.
     */
    fmat4 perspective_matrix(float fovyInRadians, float aspectRatio, float znear, float zfar);

    /**
     * @brief Creates an orthographic projection matrix.
     * 
     * The matrix is stored in column-major order.
     * @param matrix Pointer to a 16-element array to receive the matrix.
     * @param l Left vertical clipping plane.
     * @param r Right vertical clipping plane.
     * @param b Bottom horizontal clipping plane.
     * @param t Top horizontal clipping plane.
     * @param n Distance to the near depth clipping plane.
     * @param f Distance to the far depth clipping plane.
     */
    void ortho(float* matrix, float l, float r, float b, float t, float n, float f);

    /**
     * @brief Constructs and returns an orthographic projection matrix.
     * 
     * @param l Left vertical clipping plane.
     * @param r Right vertical clipping plane.
     * @param b Bottom horizontal clipping plane.
     * @param t Top horizontal clipping plane.
     * @param n Distance to the near depth clipping plane.
     * @param f Distance to the far depth clipping plane.
     * @return fmat4 The resulting orthographic projection matrix.
     */
    fmat4 ortho(float l, float r, float b, float t, float n, float f);
}

#endif // __gl_math_h__
