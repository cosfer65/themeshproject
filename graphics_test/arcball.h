#ifndef __arcball_h__
#define __arcball_h__
#include "quaternion.h"

using namespace base_math;

/// \class arcball
/// \brief Implements a virtual trackball for intuitive 3D view rotation using mouse drags.
///
/// The `arcball` maps 2D mouse positions to points on a virtual unit sphere and
/// accumulates the induced rotations as a quaternion. Typical usage:
/// - construct with the current viewport size,
/// - call `beginDrag` on mouse-press,
/// - call `drag` on mouse-move while pressed,
/// - call `endDrag` on mouse-release,
/// - query the current rotation via `rotation`.
class arcball {
public:
    /// \brief Construct an arcball controller for a given viewport size.
    ///
    /// \param width  Current viewport width in pixels.
    /// \param height Current viewport height in pixels.
    ///
    /// Initializes the internal quaternion to identity and disables dragging.
    // Constructor takes the viewport size (width, height), used to map mouse coordinates into a normalized [-1, 1] range.
    arcball(float width, float height)
        : W(width), H(height), dragging(false), m_rotation(1, 0, 0, 0) {
    }

    /// \brief Update the viewport size used to normalize mouse coordinates.
    ///
    /// Call this whenever the window is resized so that mouse positions are
    /// correctly mapped to the virtual sphere.
    /// \param width  New viewport width in pixels.
    /// \param height New viewport height in pixels.
    void resize(float width, float height) {
        W = width; H = height;
    }

    /// \brief Start a drag operation from the given mouse position.
    ///
    /// Records the starting point on the virtual sphere and enables dragging.
    /// \param mx Mouse X position in pixels (window coordinates).
    /// \param my Mouse Y position in pixels (window coordinates).
    void beginDrag(float mx, float my) {
        dragging = true;
        v0 = project(mx, my);
    }

    /// \brief Update the rotation based on a drag to the given mouse position.
    ///
    /// If dragging is active, the previous and current projected mouse
    /// positions define an incremental rotation around their cross-product
    /// axis. The resulting quaternion is left-multiplied into the accumulated
    /// rotation.
    ///
    /// \param mx Current mouse X position in pixels.
    /// \param my Current mouse Y position in pixels.
    void drag(float mx, float my) {
        if (!dragging) return;

        // Convert mouse position to normalized device coordinates(NDC) in [-1, 1]
        fvec3 v1 = project(mx, my);

        // Computing the Rotation Quaternion
        // Angle of rotation: angle = acos(dot(v0, v1))
        float dot = std::max<float>(-1.0f, std::min<float>(1.0f, v0.dot(v1)));
        float angle = std::acos(dot);

        // Axis of rotation: axis = v0.cross(v1).normalize() prev to current
        fvec3 axis = v0.cross(v1).normalize();
        if (axis.x() == 0 && axis.y() == 0 && axis.z() == 0)
            return;

        // Create rotation quaternion from axis-angle representation
        quaternion<float> q(fromAxisAngle<float>(axis, angle));

        // Update the current rotation
        m_rotation = q * m_rotation;
        m_rotation.normalize();

        // Update v0 for the next drag event
        v0 = v1;
    }

    /// \brief Finish the current drag operation.
    ///
    /// Subsequent calls to `drag` will have no effect until `beginDrag` is
    /// called again.
    void endDrag() {
        dragging = false;
    }

    /// \brief Get the current rotation as a column-major 4x4 matrix.
    ///
    /// \param M Pointer to an array of 16 floats to receive the rotation matrix.
    ///          The layout is determined by `quaternion::matrix(float[16])`.
    void rotation(float M[16]) const {
        m_rotation.matrix(M);
    }

    /// \brief Get the current rotation as an `fmat4`.
    ///
    /// \return A 4x4 float matrix representing the accumulated rotation.
    fmat4 rotation() const {
        return m_rotation.matrix();
    }

    /// \brief Get the current rotation as a unit quaternion.
    ///
    /// \return A const reference to the accumulated rotation quaternion.
    ///         This quaternion represents the same orientation as the matrix
    ///         returned by `rotation()`.
    const quaternion<float>& get_quaternion() const {
        return m_rotation;
    }

private:
    /// \brief Flag indicating whether a drag operation is active.
    bool dragging;

    /// \brief Accumulated orientation represented as a unit quaternion.
    ///
    /// This quaternion is updated incrementally on each `drag` call.
    quaternion<float> m_rotation;

    /// \brief Viewport width in pixels used to normalize mouse coordinates.
    float W;

    /// \brief Viewport height in pixels used to normalize mouse coordinates.
    float H;

    /// \brief Last projected mouse position on the virtual sphere.
    fvec3 v0;

    // Mapping Mouse Coordinates to the Sphere
    // 1. Normalize screen coordinates
    // if the window is WxH:
    // 	x = 2*mousex/W   and y = 1 - 2*mousey/H
    // 	
    // 2. Project onto a sphere
    // let d = x^2+y^2
    // 	if d <= 1:
    // 		z = sqrt(1-d)
    // 		v = fvec3(x,y,z);
    // 	else (if d > 1)
    // 		v = fvec3(x,y,0).normalize()
    /// \brief Project a mouse position to a point on the virtual unit sphere.
    ///
    /// The mouse coordinates are first mapped to normalized device coordinates
    /// in [-1, 1] using the current viewport size. Points inside the unit
    /// circle are lifted onto the hemisphere; points outside are projected
    /// back onto the circle and normalized.
    ///
    /// \param mx Mouse X position in pixels.
    /// \param my Mouse Y position in pixels.
    /// \return A 3D vector on or on the boundary of the unit sphere.
    fvec3 project(float mx, float my) const {
        float x = 2.0f * mx / W - 1.0f;
        float y = 1.0f - 2.0f * my / H;
        float d = x * x + y * y;

        float z = d <= 1.0f ? std::sqrt(1.0f - d) : 0.0f;

        if (d > 1.0f) {
            float len = std::sqrt(d);
            x /= len;
            y /= len;
        }
        return fvec3(x, y, z);
    }
};

#endif // __arcball_h__

