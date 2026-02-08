/// \file gl_camera.h
/// \brief Declares the `base_opengl::gl_camera` class, a simple OpenGL camera wrapper.

#ifndef __camera_h__
#define __camera_h__

#include "vector.h"
#include "gl_math.h"

namespace base_opengl {

    /// \brief Simple perspective camera for OpenGL rendering.
    ///
    /// Encapsulates camera position, orientation and projection parameters,
    /// and provides convenience methods to build view/projection matrices,
    /// configure the OpenGL viewport and perform basic camera panning.
    class gl_camera {
        /// \brief Camera position in world space.
        base_math::fvec3 location;

        /// \brief Point in world space the camera is looking at.
        base_math::fvec3 target;

        /// \brief Up direction of the camera in world space.
        base_math::fvec3 up = base_math::fvec3(0, 1, 0);

        /// \brief Vertical field of view in radians.
        float fov = PI<float> / 6;

        /// \brief Distance to the near clipping plane.
        float nearPlane = 0.1f;

        /// \brief Distance to the far clipping plane.
        float farPlane = 1000.0f;

        /// \brief Aspect ratio of the viewport (width / height).
        float aspect = 1.f;

        /// \brief Window width and height in pixels.
        int width, height;  ///< Window width and height.

        /// \brief Left coordinate of the viewport within the window.
        int left = 0;

        /// \brief Bottom coordinate of the viewport within the window.
        int bottom = 0; ///< Viewport position within the window.
    public:

        /// \brief Constructs a default camera with uninitialized position/target.
        gl_camera() = default;

        /// \brief Constructs a camera with the given position, target and up vector.
        ///
        /// \param _location Camera position in world space.
        /// \param _target   Look-at target point in world space.
        /// \param _up       Up direction vector (defaults to +Y).
        gl_camera(const base_math::fvec3& _location,
                  const base_math::fvec3& _target,
                  const base_math::fvec3& _up = base_math::fvec3(0, 1, 0)) :
            location(_location), target(_target), up(_up) {
        }

        /// \brief Sets the viewport size and updates the aspect ratio.
        ///
        /// \param _width  Viewport width in pixels.
        /// \param _height Viewport height in pixels.
        void set_aspect(int _width, int _height) {
            width = _width;
            height = _height;
            aspect = float(width) / float(height);
        }

        /// \brief Initializes camera position, target and up vectors.
        ///
        /// \param _location Camera position in world space.
        /// \param _target   Look-at target point in world space.
        /// \param _up       Up direction vector (defaults to +Y).
        void setup(const base_math::fvec3& _location,
                   const base_math::fvec3& _target,
                   const base_math::fvec3& _up = base_math::fvec3(0, 1, 0))
        {
            location = _location;
            target = _target;
            up = _up;
        }

        /// \brief Sets the vertical field of view.
        ///
        /// \param fov_in_radians Field of view in radians.
        void set_fov(float fov_in_radians) {
            fov = fov_in_radians;
        }

        /// \brief Sets viewport rectangle and updates aspect ratio.
        ///
        /// \param _left   Left coordinate of the viewport within the window.
        /// \param _bottom Bottom coordinate of the viewport within the window.
        /// \param _width  Viewport width in pixels.
        /// \param _height Viewport height in pixels.
        void set_viewport(int _left, int _bottom, int _width, int _height) {
            left = _left;
            bottom = _bottom;
            width = _width;
            height = _height;
            aspect = float(width) / float(height);
        }

        /// \brief Sets the near and far clipping planes.
        ///
        /// \param nearP Distance to the near plane.
        /// \param farP  Distance to the far plane.
        void set_depth_range(float nearP, float farP) {
            nearPlane = nearP;
            farPlane = farP;
        }

        /// \brief Returns the normalized forward direction vector.
        ///
        /// Defined as the normalized vector from `location` to `target`.
        base_math::fvec3 forward() const { return (target - location).normalize(); }

        /// \brief Returns the normalized right direction vector.
        ///
        /// Computed as the cross product of `forward()` and `up`.
        base_math::fvec3 right()   const { return forward().cross(up).normalize(); }

        /// \brief Returns the camera up vector derived from forward and right.
        ///
        /// Computed as the cross product of `right()` and `forward()`.
        base_math::fvec3 upVec()   const { return right().cross(forward()); }

        /// \brief Returns the distance between camera location and target.
        float distance() const { return (target - location).length(); }

        /// \brief Builds the perspective projection matrix for the camera.
        ///
        /// Uses the current `fov`, `aspect`, `nearPlane` and `farPlane`.
        base_math::fmat4 view_perspective() {
            return perspective_matrix(fov, aspect, nearPlane, farPlane);
        }

        /// \brief Builds the view matrix using a look-at transform.
        ///
        /// Uses the current `location`, `target` and `up`.
        base_math::fmat4 camera_perspective() {
            return lookAt(location, target, up);
        }

        /// \brief Returns the combined view-projection matrix.
        ///
        /// Computed as `camera_perspective() * view_perspective()`.
        base_math::fmat4 perspective() {
            return camera_perspective() * view_perspective();
        }

        /// \brief Applies the currently configured viewport to OpenGL.
        ///
        /// Calls `glViewport(left, bottom, width, height)`.
        void set_viewport() {
            glViewport(left, bottom, width, height);
        }

        /// \brief Pans the camera parallel to the view plane.
        ///
        /// Translates both `location` and `target` along the right and up
        /// directions, scaled by the current distance to the target.
        ///
        /// \param dx Horizontal pan amount (screen space units).
        /// \param dy Vertical pan amount (screen space units).
        void pan(float dx, float dy) {
            float s = distance() * 0.002f;

            base_math::fvec3 delta = (-dx * right() + dy * upVec()) * s;
            delta *= 0.1f; // adjust the pan speed
            location += delta;
            target += delta;
        }

        void zoom(float delta) {
            float s = distance() * 0.001f;
            base_math::fvec3 delta_vec = forward() * delta * s;
            location += delta_vec;
        }

        /// \brief Returns the current camera position in world space.
        const base_math::fvec3& get_location() const { return location; }

        /// \brief Returns the current look-at target in world space.
        const base_math::fvec3& get_target() const { return target; }
    };
}

#endif // __camera_h__
