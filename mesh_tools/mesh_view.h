/// \file mesh_view.h
/// \brief Declares the `mesh_view` class, responsible for visualizing 3D mesh data.
///
/// This header defines the public interface for the mesh viewing component,
/// including camera configuration, rendering entry points, and user interaction
/// hooks (mouse and command handling). The actual rendering and mesh-processing
/// logic are implemented in the corresponding source file.

#ifndef __mesh_view_h__
#define __mesh_view_h__

/// \brief Forward declaration of the private implementation for `mesh_view`.
///
/// The `mesh_view_private` struct contains the internal data structures and
/// implementation details needed by `mesh_view`. It is intentionally hidden
/// from users of this header to keep the public interface minimal and to reduce
/// compile-time dependencies (PIMPL idiom).
struct mesh_view_private;

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
    mesh_view_private* m_private;                     ///< Private implementation details

    /// \brief Field of view, in degrees, for the perspective projection.
    ///
    /// This value controls how "zoomed in" the view is. It is adjusted by
    /// mouse wheel events in `mouse_wheel()` and used when configuring the
    /// camera/projection matrix during rendering.
    float fov = 15.f;                                 ///< Field of view in degrees (changed with mouse wheel and affects zoom)

    /// \brief Indicates whether mouse events should currently be ignored.
    ///
    /// When `true`, mouse input handlers (such as `mouse_move()` and
    /// `mouse_wheel()`) should treat incoming events as blocked. This is
    /// typically enabled during operations like loading a new model to
    /// avoid inconsistent state while resources are being updated.
    bool block_mouse_event = false;                   ///< Block mouse event flag (used when loading a new model)

    /// \brief Builds or updates the visualization of mesh curvature.
    ///
    /// Uses precomputed curvature data (see `do_curvature_calculation()`)
    /// to generate buffers or other renderable resources needed to display
    /// curvature as a shading, color map, or overlay on the mesh surface.
    void generate_model_curvature_view();

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

    /// \brief Updates the current rotation and panning offsets of the viewed objects.
    ///
    /// Applies any accumulated user interaction (e.g., mouse drag, orbit,
    /// or pan operations) to the internal camera/object transform state.
    /// This method is typically invoked from input handlers or the render
    /// loop to convert raw input deltas into stable rotation and translation
    /// values used when drawing the mesh.
    void update_objects_rotation_and_pan();

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
};

#endif // __mesh_view_h__
