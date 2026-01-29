#ifndef __camera_h__
#define __camera_h__

#include "vector.h"
#include "gl_math.h"

using namespace base_math;

namespace base_opengl {

    /**
     * @class gl_camera
     * @brief Represents a 3D camera for OpenGL rendering, supporting movement and rotation.
     *
     * The gl_camera class encapsulates camera position, orientation, and provides
     * methods for moving, rotating, and generating a view matrix.
     */
    class gl_camera {
        /**
         * @brief Normalizes the camera's forward and up vectors.
         * 
         * Ensures the direction vectors remain unit length after transformations.
         */
        void update_camera() {
            vForward.normalize();
            vUp.normalize();
        }
    public:
        /// Camera position in world coordinates.
        fvec3 vLocation;
        /// Forward direction vector.
        fvec3 vForward;
        /// Up direction vector.
        fvec3 vUp;

        /**
         * @brief Default constructor. Initializes camera at (0,0,3) looking at (0,0,0) with up (0,1,0).
         */
        gl_camera() {
            setup(fvec3(0, 0, 3), fvec3(0, 0, 0), fvec3(0, 1, 0));
        }

        /**
         * @brief Constructs a camera with specified position, forward, and up vectors.
         * @param _l Camera position.
         * @param _f Forward direction.
         * @param _u Up direction.
         */
        gl_camera(const fvec3& _l, const fvec3 _f, const fvec3& _u) {
            setup(_l, _f, _u);
        }

        /**
         * @brief Copy constructor.
         * @param cp Camera to copy from.
         */
        gl_camera(gl_camera& cp) {
            setup(cp.vLocation, cp.vForward, cp.vUp);
        }

        /**
         * @brief Destructor.
         */
        ~gl_camera() {
        }

        /**
         * @brief Returns the view matrix for the camera.
         * @return View matrix (lookAt).
         */
        fmat4 perspective() {
            return lookAt(vLocation, vLocation + vForward, vUp);
        }

        /**
         * @brief Moves the camera by the specified offsets.
         * @param fwdx Offset along X.
         * @param fwdy Offset along Y.
         * @param fwdz Offset along Z.
         */
        void moveby(float fwdx, float fwdy, float fwdz) {
            vLocation.x() += fwdx;
            vLocation.y() += fwdy;
            vLocation.z() += fwdz;
            setup(vLocation, vForward, vUp);
        }

        /**
         * @brief Moves the camera to the specified position.
         * @param x New X position.
         * @param y New Y position.
         * @param z New Z position.
         */
        void moveto(float x, float y, float z) {
            vLocation.x() = x;
            vLocation.y() = y;
            vLocation.z() = z;
            setup(vLocation, vForward, vUp);
        }

        /**
         * @brief Rotates the camera around its local axes.
         * @param x Angle to rotate around right axis.
         * @param y Angle to rotate around up axis.
         * @param z Angle to rotate around forward axis.
         */
        void rotate(float x, float y, float z) {
            // x - right
            // y - up
            // z - forward
            rotate_locX(x);
            rotate_locY(y);
            rotate_locZ(z);
        }

        //void rotate_vec(const fvec3& v, float x, float y, float z);

        /**
         * @brief Rotates the camera around its local X (right) axis.
         * @param angle Angle in radians.
         */
        void rotate_locX(float angle);

        /**
         * @brief Rotates the camera around its local Y (up) axis.
         * @param angle Angle in radians.
         */
        void rotate_locY(float angle);

        /**
         * @brief Rotates the camera around its local Z (forward) axis.
         * @param angle Angle in radians.
         */
        void rotate_locZ(float angle);

        /**
         * @brief Sets up the camera using position, look-at, and up vectors.
         * @param pos Camera position.
         * @param lookat Target point to look at.
         * @param up Up direction.
         */
        void setup(const fvec3& pos, const fvec3& lookat, const fvec3& up);

        /**
         * @brief Sets up the camera using position, forward, and up vectors.
         * @param pos Camera position.
         * @param forw Forward direction.
         * @param up Up direction.
         */
        void setup_d(const fvec3& pos, const fvec3& forw, const fvec3& up);

        /**
         * @brief Sets up the camera using current member variables.
         */
        void setup_d();
    };

    /**
     * @class gl_viewport
     * @brief Manages OpenGL viewport and projection parameters.
     *
     * The gl_viewport class encapsulates viewport size, aspect ratio, field of view,
     * and provides methods for configuring the OpenGL viewport and projection matrix.
     */
    class gl_viewport {
        int width, height;  ///< Window width and height.
        float near_plane, far_plane; ///< Near and far clipping planes.
        float zoom_angle; ///< Field of view angle in degrees.
        float aspect; ///< Aspect ratio (width/height).
        int left, bottom; ///< Viewport position within the window.

    public:
        /**
         * @brief Constructs a viewport with specified parameters.
         * @param _pa Field of view angle in degrees.
         * @param _w Window width.
         * @param _h Window height.
         * @param _n Near clipping plane.
         * @param _f Far clipping plane.
         */
        gl_viewport(float _pa = 45, int _w = 800, int _h = 600, float _n = 0.1f, float _f = 1000.f) :
            zoom_angle(_pa), near_plane(_n), far_plane(_f), width(_w), height(_h), aspect((float)_w / (float)_h),
            left(0), bottom(0) {
        }

        /**
         * @brief Returns the perspective projection matrix.
         * @return Perspective projection matrix.
         */
        fmat4 perspective() {
            return perspective_matrix(zoom_angle, aspect, near_plane, far_plane);
        }

        /**
         * @brief Sets the perspective projection parameters.
         * @param angle Field of view angle in degrees.
         * @param np Near clipping plane.
         * @param fp Far clipping plane.
         */
        void set_perspective(float angle, float np, float fp) {
            zoom_angle = angle;
            near_plane = np;
            far_plane = fp;
        }

        // fmat4 ortho() {
        //     return ortho_matrix(zoom_angle, aspect, near_plane, far_plane);
        // }

        /**
         * @brief Sets the window size and updates the aspect ratio.
         * @param _w Window width.
         * @param _h Window height.
         */
        void set_window_aspect(int _w, int _h) {
            width = _w;
            height = _h;
            aspect = (float)width / (float)height;
        }

        /**
         * @brief Sets the viewport position within the window.
         * @param _left Left position.
         * @param _bottom Bottom position.
         */
        void set_position(int _left, int _bottom) {
            left = _left;
            bottom = _bottom;
        }

        /**
         * @brief Sets the field of view angle.
         * @param _z Field of view angle in degrees.
         */
        void set_fov(float _z) {
            zoom_angle = _z;
        }

        /**
         * @brief Sets the near and far clipping planes.
         * @param _n Near clipping plane.
         * @param _f Far clipping plane.
         */
        void set_view_field(float _n, float _f) {
            near_plane = _n;
            far_plane = _f;
        }

        /**
         * @brief Applies the viewport settings to OpenGL.
         */
        void set_viewport() {
            glViewport(left, bottom, width, height);
        }
    };
}
#endif // __camera_h__
