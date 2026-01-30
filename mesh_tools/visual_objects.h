#ifndef __visual_objects_h__
#define __visual_objects_h__

#include "vector.h"

/// @brief Forward declarations for core OpenGL helper types used by visual objects.
namespace base_opengl {
    /// @brief Represents a simple GPU-resident mesh (vertex/index buffers, layout, etc.).
    class gl_mesh;
    /// @brief Represents a drawable OpenGL primitive associated with a mesh and material/shader state.
    class gl_prim;
    /// @brief Represents an OpenGL shader program used for rendering.
    class gl_shader;
}

/// @brief Manages a collection of small auxiliary visual objects rendered on top of a main scene.
///
/// The `visual_objects` class owns or references a simple mesh and an associated
/// OpenGL primitive used to visualize helper elements such as debug vectors or
/// other miscellaneous objects. Visibility can be toggled at runtime.
class visual_objects {
    base_opengl::gl_mesh* m_mesh; ///< Pointer to the associated simple mesh.
    base_opengl::gl_prim* m_prim; ///< Pointer to the associated OpenGL primitive.
    bool m_visible;    ///< Visibility flag for the miscellaneous objects.
public:
    /// @brief Constructs an empty visual object container with no mesh/primitive attached.
    ///
    /// The internal mesh/primitive pointers are initialized, and the visibility flag
    /// is set to a default state (implementation-defined, typically visible).
    visual_objects();

    /// @brief Destroys the visual object container and releases any owned GPU resources.
    ///
    /// Cleans up the associated mesh/primitive if they are owned by this object.
    ~visual_objects();

    /// @brief Adds a line segment defined by two 3D points to the underlying mesh.
    ///
    /// This is typically used to add debug vectors or simple guide lines.
    ///
    /// @param vstart Start position of the vector in world or object space.
    /// @param vend   End position of the vector in world or object space.
    /// @return An implementation-defined identifier or index of the added vector,
    ///         or a negative value on failure.
    int add_vector(const base_math::fvec3& vstart, const base_math::fvec3& vend);

    /// @brief Creates and initializes the OpenGL primitive that will draw the mesh.
    ///
    /// Sets up the `gl_prim` instance to reference the internal mesh and prepares it
    /// for rendering with a shader.
    ///
    /// @return Pointer to the created `gl_prim` instance. Ownership semantics are
    ///         implementation-defined but typically remain with this object.
    base_opengl::gl_prim* create_prim();

    /// @brief Returns the underlying OpenGL primitive used for rendering.
    ///
    /// This allows external code to further configure the primitive (e.g., draw mode).
    ///
    /// @return The `gl_prim` instance associated with this visual object, or `nullptr`
    ///         if it has not been created yet.
    base_opengl::gl_prim* get_prim() {
        return m_prim;
    }

    /// @brief Renders all visual objects using the given shader.
    ///
    /// If the visibility flag is disabled or the primitive is not initialized,
    /// this call may be a no-op.
    ///
    /// @param sh Shader program to use when drawing the underlying primitive.
    void render(base_opengl::gl_shader* sh);

    /// @brief Toggles the visibility state of the visual objects.
    ///
    /// When visibility is disabled, `render()` should skip drawing.
    void toggle_visibility() {
        m_visible = !m_visible;
    }
};

/// @brief Forward declaration of implementation details for `UCS_view`.
///
/// The actual data and logic are hidden behind a PIMPL (pointer-to-implementation)
/// to keep the header lightweight and reduce compile-time dependencies.
struct UCS_view_private;

/// @brief Manages visualization and interaction with a UCS (User Coordinate System) view.
///
/// `UCS_view` is responsible for initializing, rendering, and updating the orientation
/// of a local coordinate system, typically rendered as a small axis gizmo.
class UCS_view {
    UCS_view_private* m_private_data; ///< Opaque pointer to implementation-specific data.
public:
    /// @brief Constructs an uninitialized UCS view.
    ///
    /// The underlying implementation data is allocated but not fully initialized
    /// until `initialize()` is called.
    UCS_view();

    /// @brief Destroys the UCS view and releases implementation data and GPU resources.
    ~UCS_view();

    /// @brief Initializes OpenGL resources and internal state required for rendering the UCS.
    ///
    /// Must be called before `render()` or any rotation methods are used.
    void initialize();

    /// @brief Renders the UCS indicator to the current OpenGL context.
    ///
    /// Typically called once per frame after the main scene has been drawn.
    void render(const base_math::fmat4& user_rotation);

    /// @brief Updates the UCS view to match a new window size.
    ///
    /// This is usually called from a window-resize callback to adjust viewport,
    /// projection, or placement of the UCS overlay.
    ///
    /// @param width  New width of the rendering surface, in pixels.
    /// @param height New height of the rendering surface, in pixels.
    void resize_window(int width, int height);

    /// @brief Applies an incremental rotation to the UCS around the X, Y, and Z axes.
    ///
    /// @param x Rotation delta around the X-axis, in degrees (or radians, implementation-defined).
    /// @param y Rotation delta around the Y-axis, in degrees (or radians, implementation-defined).
    /// @param z Rotation delta around the Z-axis, in degrees (or radians, implementation-defined).
    void rotate_ucs_by(float x, float y, float z);

    /// @brief Sets the absolute orientation of the UCS around the X, Y, and Z axes.
    ///
    /// @param x Target rotation around the X-axis, in degrees (or radians, implementation-defined).
    /// @param y Target rotation around the Y-axis, in degrees (or radians, implementation-defined).
    /// @param z Target rotation around the Z-axis, in degrees (or radians, implementation-defined).
    void rotate_ucs_to(float x, float y, float z);
};


#endif // __visual_objects_h__
