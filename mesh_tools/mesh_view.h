/// \file mesh_view.h
/// \brief Declares the `mesh_view` class, responsible for visualizing 3D mesh data.
///
/// This header defines the public interface for the mesh viewing component,
/// including camera configuration, rendering entry points, and user interaction
/// hooks (mouse and command handling). The actual rendering and mesh-processing
/// logic are implemented in the corresponding source file.

#ifndef __mesh_view_h__
#define __mesh_view_h__

#include "matrix.h"
// #include "gl_shaders.h"

// namespace base_math {
//     class fmat4; ///< Forward declaration of a 4x4 float matrix class.
// }
namespace base_opengl {
    class gl_shader; ///< Forward declaration of an OpenGL shader wrapper class.
}

/// \brief Forward declaration of the private implementation for `mesh_view`.
///
/// The `mesh_view_private` struct contains the internal data structures and
/// implementation details needed by `mesh_view`. It is intentionally hidden
/// from users of this header to keep the public interface minimal and to reduce
/// compile-time dependencies (PIMPL idiom).
struct mesh_view_private;

struct view_state {
    bool show_wireframe = false;  ///< Whether to render the mesh in wireframe mode.
    bool show_curvature = false; ///< Whether to visualize curvature directions.
    bool show_normals = false;   ///< Whether to visualize face normals.
    bool show_gauss_map = false;   ///< Whether to visualize the Gaussian curvature map.

    void set_model_changed() {
        show_curvature = false; // disable curvature visualization until it is regenerated for the new model
        show_normals = false;   // disable normals visualization until it is regenerated for the new model
        show_wireframe = false;  // disable wireframe visualization until the user toggles it
        show_gauss_map = false;  // disable Gauss map visualization until the user toggles it
    }
};

struct model_state {
    bool m_curvatures_calculated = false; ///< Flag indicating whether curvature has been calculated for the current model.
    bool m_gauss_curvature_view_generated = false; ///< Flag indicating whether the Gaussian curvature visualization has been generated for the current model.
    void set_model_changed() {
            m_curvatures_calculated = false; // reset curvature calculation state when the model changes
            m_gauss_curvature_view_generated = false; // reset Gauss curvature view state when the model changes
    }
};

/// \brief Viewer and controller for displaying and interacting with a 3D mesh.
///
/// `mesh_view` encapsulates the logic for:
/// - Initializing the rendering context and associated resources.
/// - Rendering the current mesh with optional visualization modes
///   (e.g., curvature, normals).
/// - Responding to window resize events.
/// - Handling mouse input for camera control and interaction.
/// - Loading and switching between mesh models.
///
/// Typical usage:
/// - Construct an instance.
/// - Call `initialize()` once before any rendering.
/// - On each frame, call `render()`.
/// - Route window size and input events to `resize()`, `mouse_move()`,
///   `mouse_wheel()`, and `onCommand()`.
class mesh_view {
    /// \brief Pointer to the private implementation.
    ///
    /// All heavy state (mesh data, OpenGL objects, shaders, etc.) is stored
    /// in this private object. This keeps the header stable and reduces
    /// rebuilds when implementation details change.
    mesh_view_private* m_private;                    ///< Private implementation details

    view_state m_view_state;                         ///< Current visualization state (e.g., whether to show curvature or normals)
    model_state m_model_state;                       ///< State related to the currently loaded model (e.g., whether curvature has been calculated)

    bool dragging = false;                           ///< Indicates whether the user is currently dragging with the mouse (for panning)
    int last_mouse_x = 0;                            ///< Last recorded mouse X position (used for calculating deltas during dragging)
    int last_mouse_y = 0;                            ///< Last recorded mouse Y position (used for calculating deltas during dragging)

    /// \brief Indicates whether mouse events should currently be ignored.
    ///
    /// When `true`, mouse input handlers (such as `mouse_move()` and
    /// `mouse_wheel()`) should treat incoming events as blocked. This is
    /// typically enabled during operations like loading a new model to
    /// avoid inconsistent state while resources are being updated.
    bool block_mouse_event = false;                   ///< Block mouse event flag (used when loading a new model)

    /// \brief Constructs or refreshes the core GPU-side representation of the mesh.
    ///
    /// Uploads the mesh geometry and any required per-vertex/per-face
    /// attributes (positions, normals, indices, etc.) into graphics API
    /// resources such as vertex/index buffers or vertex array objects.
    /// This is typically called after loading a new model or when the
    /// mesh data has changed to ensure subsequent rendering uses an
    /// up-to-date representation.
    void build_model_representation();

    /// \brief Builds or updates the visualization of mesh curvature.
    ///
    /// Uses precomputed curvature data (see `do_curvature_calculation()`)
    /// to generate buffers or other renderable resources needed to display
    /// curvature as a shading, color map, or overlay on the mesh surface.
    void generate_model_curvature_view();

    /// \brief Builds or updates the Gaussian curvature visualization for the current mesh.
    ///
    /// Uses curvature data associated with the loaded model to construct a
    /// renderable representation of the Gaussian curvature field. Typical
    /// outputs include color-mapped surfaces or auxiliary geometry that
    /// can be toggled via the view state (see `view_state::show_gauss_map`).
    /// This method assumes curvature information is already available;
    /// if not, it may trigger or depend on prior curvature computation.
    void generate_model_gauss_curvature_view();

    /// \brief Computes curvature information for the currently loaded mesh.
    ///
    /// Performs the numerical or geometric analysis necessary to estimate
    /// curvature per vertex or per face. The results are stored in the
    /// private implementation and can later be consumed by
    /// `generate_model_curvature_view()`.
    void do_curvature_calculation();

    /// \brief Generates visualization data for mesh normals.
    ///
    /// Creates renderable representations of vertex or face normals, such as
    /// line segments or color encodings, which can be used to debug or inspect
    /// the orientation of the mesh surface.
    void generate_model_normals();

public:
    /// \brief Constructs a new `mesh_view` with default configuration.
    ///
    /// The constructor sets up initial internal state but does not allocate
    /// GPU resources or prepare the rendering context. Call `initialize()`
    /// before the first call to `render()`.
    mesh_view();

    /// \brief Destroys the `mesh_view` instance and releases resources.
    ///
    /// Ensures that any resources owned by the private implementation,
    /// such as GPU buffers and loaded mesh data, are properly released.
    virtual ~mesh_view();

    /// \brief Initializes rendering-related resources for the mesh view.
    ///
    /// Call this once after the rendering context is available and before
    /// performing any rendering. Typical responsibilities include creating
    /// shaders, vertex buffers, and any other GPU-side state required by
    /// the view.
    void initialize();

    /// \brief Renders the 3D scene using externally provided camera and rotation matrices.
    ///
    /// \param cam_matrix Reference to a 4x4 camera/view matrix that defines the camera's
    ///                   position and orientation in the scene.
    /// \param rot_mat    Reference to a 4x4 rotation matrix that controls the orientation
    ///                   of the mesh being rendered.
    /// \param shdr       Pointer to the OpenGL shader program to use for rendering the scene.
    ///
    /// This is an advanced rendering method that allows external control over the
    /// camera transformation, object rotation, and shader selection. Unlike `render()`,
    /// which uses internal camera state and default shaders, this method accepts
    /// pre-configured matrices and shader, making it suitable for multi-pass rendering,
    /// custom camera setups, or integration with external rendering pipelines.
    void render_scene(base_math::fmat4& cam_matrix, base_math::fmat4& rot_mat, base_opengl::gl_shader* shdr);

    /// \brief Renders the currently loaded mesh to the active render target.
    ///
    /// Should be called once per frame after `initialize()` has been
    /// successfully invoked. Uses the current camera settings, including
    /// field of view (`fov`), and any active visualization modes such as
    /// curvature or normal display.
    void render();

    /// \brief Loads a new mesh model from the given file path.
    ///
    /// \param fnm Null-terminated path to the mesh file to load.
    /// \return `true` if the model was loaded successfully; otherwise `false`.
    ///
    /// On success, internal buffers and visualization data (e.g., curvature
    /// view, normals) may be regenerated. Mouse events may be temporarily
    /// blocked via `block_mouse_event` while the load is in progress.
    bool load_new_model(const char* fnm);

    /// \brief Notifies the view that the output surface has been resized.
    ///
    /// \param width  New width of the viewport or window, in pixels.
    /// \param height New height of the viewport or window, in pixels.
    ///
    /// This should be called whenever the hosting window is resized so that
    /// the projection matrix and viewport can be updated accordingly.
    void resize(int width, int height);

    /// \brief Handles high-level commands directed at the mesh view.
    ///
    /// \param cmd Integer identifier representing the command to execute.
    /// \return Implementation-defined integer result (e.g., success code).
    ///
    /// This function provides a generic hook for menu actions, toolbar
    /// commands, or other UI triggers that affect the mesh viewer, such as
    /// toggling visualization modes or resetting the camera.
    int onCommand(int cmd);

    /// \brief Handles mouse move events over the mesh view.
    ///
    /// \param x      Current mouse X position in window/client coordinates.
    /// \param y      Current mouse Y position in window/client coordinates.
    /// \param extra  Additional mouse state, typically a bitmask containing
    ///               button and/or modifier key flags provided by the windowing system.
    ///
    /// This is typically used to update camera orientation or perform
    /// interactive operations (e.g., rotation, panning, selection) while one
    /// or more mouse buttons are pressed. When mouse input is blocked
    /// (see `block_mouse_event`), this handler should ignore incoming events.
    void onMouseMove(int x, int y, unsigned __int64 extra);

    /// \brief Handles left mouse button press events.
    ///
    /// \param x      Mouse X position at the time of the click.
    /// \param y      Mouse Y position at the time of the click.
    /// \param extra  Additional mouse state, typically a bitmask of modifier
    ///               keys or other button states.
    ///
    /// Common uses include starting camera interactions (orbit / pan),
    /// initiating drag operations, or beginning mesh element selection.
    void onLMouseDown(int x, int y, unsigned __int64 extra);

    /// \brief Handles left mouse button release events.
    ///
    /// \param x      Mouse X position at the time the button is released.
    /// \param y      Mouse Y position at the time the button is released.
    /// \param extra  Additional mouse state, typically a bitmask of modifier
    ///               keys or other button states.
    ///
    /// Used to finalize drag/rotate operations, complete selection, and
    /// generally stop any interaction that was started in `onLMouseDown()`.
    void onLMouseUp(int x, int y, unsigned __int64 extra);

    /// \brief Processes mouse wheel events to adjust zoom or field of view.
    ///
    /// \param delta      Wheel delta value indicating scroll direction and
    ///                   magnitude.
    /// \param extra_btn  Bitmask or flag set representing extra mouse buttons
    ///                   and/or modifier keys active during the event.
    ///
    /// Typically used to update `fov`, resulting in zooming in or out of
    /// the displayed mesh. When `block_mouse_event` is `true`, the event
    /// may be ignored until interaction is re-enabled.
    void onMouseWheel(int delta, unsigned __int64 extra_btn);


    /// \brief Handles right mouse button press events.
    ///
    /// \param x      Mouse X position in window/client coordinates at the time
    ///               the right button is pressed.
    /// \param y      Mouse Y position in window/client coordinates at the time
    ///               the right button is pressed.
    /// \param extra  Additional mouse state, typically a bitmask of modifier
    ///               keys and/or other button states reported by the windowing
    ///               system.
    ///
    /// This is typically used to start secondary interactions such as opening
    /// a context menu, initiating an alternative camera control mode, or
    /// beginning a secondary drag operation.
    void onRMouseDown(int x, int y, unsigned __int64 extra);

    /// \brief Handles right mouse button release events.
    ///
    /// \param x      Mouse X position in window/client coordinates at the time
    ///               the right button is released.
    /// \param y      Mouse Y position in window/client coordinates at the time
    ///               the right button is released.
    /// \param extra  Additional mouse state, typically a bitmask of modifier
    ///               keys and/or other button states reported by the windowing
    ///               system.
    ///
    /// Used to finalize or cancel interactions started in `onRMouseDown`,
    /// such as stopping a secondary drag operation or confirming a
    /// context-sensitive action.
    void onRMouseUp(int x, int y, unsigned __int64 extra);

};

#endif // __mesh_view_h__
