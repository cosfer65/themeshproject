/// @file mesh_view.cpp
/// @brief Implements the `mesh_view` class responsible for displaying and interacting with 3D mesh models.
///
/// This module sets up an OpenGL-based viewport, camera, lighting, and shader pipeline to render
/// mesh models. It also provides visualization for face normals and principal curvature directions,
/// and exposes user interaction through mouse input and command handling.

#include "mesh_view.h"
#include "resource.h"
#include "arcball.h"

#include "model.h"

#include "gl_graphics.h"
#include "gl_camera.h"
#include "gl_shaders.h"
#include "gl_light.h"
#include "gl_prim.h"

#include "common_dialogs.h"
#include "visual_objects.h"

using namespace base_math;
using namespace base_opengl;

/// @brief Internal data container for `mesh_view`.
///
/// This struct encapsulates all OpenGL-related objects and visualization helpers used by `mesh_view`.
/// It is allocated on the heap and owned by `mesh_view` to hide implementation details from clients.
struct mesh_view_private {
    std::unique_ptr<gl_viewport> m_view;                ///< Viewport configuration (FOV, aspect, and window region).
    std::unique_ptr<gl_camera> m_cam;                   ///< Camera used to view the scene and compute view/projection.
    std::unique_ptr<gl_shader> m_shader;                ///< Shader program used for mesh and helper rendering.
    std::unique_ptr<gl_light> m_light;                  ///< Primary scene light affecting shading.

    std::unique_ptr<UCS_view> m_ucs_view;               ///< UCS (user coordinate system) view used to render the axes widget.

    // {{ the loaded model ant the acompanying visualization objects
    std::unique_ptr<model> m_model;                     ///< Currently loaded mesh model (can contain multiple parts).
    std::vector<std::unique_ptr<gl_prim>> m_draw_parts; ///< Drawable mesh parts converted to OpenGL primitives.
    // }}

    std::vector<std::unique_ptr<visual_objects>> m_face_normals;        ///< Visualization of per-face normals as line segments.
    std::vector<std::unique_ptr<visual_objects>> m_model_curvatures_min;///< Visualization of minimum principal curvature directions per vertex.
    std::vector<std::unique_ptr<visual_objects>> m_model_curvatures_max;///< Visualization of maximum principal curvature directions per vertex.

    std::unique_ptr<arcball> m_arcball;                                 ///< Arcball for mouse interaction
};

/// @brief Constructs a new `mesh_view` and initializes basic helpers.
///
/// Allocates the private implementation struct and creates a default viewport and UCS view.
/// OpenGL-specific state (camera, shaders, light) is created later in `initialize()`.
mesh_view::mesh_view() {
    m_private = new mesh_view_private();
    m_private->m_view.reset(new base_opengl::gl_viewport());
    m_private->m_ucs_view.reset(new UCS_view());
    m_private->m_arcball.reset(new arcball(800, 600));
}

/// @brief Destroys the `mesh_view` and releases all associated resources.
///
/// Deletes the private implementation struct which in turn releases all owned OpenGL and model objects.
mesh_view::~mesh_view() {
    delete m_private;
}

/// @brief Initializes OpenGL-related state for rendering.
///
/// This method must be called once before any rendering:
/// - Initializes the UCS view.
/// - Configures the viewport FOV using the current `fov`.
/// - Creates and positions the camera.
/// - Loads and links the mesh rendering shader program.
/// - Sets up a spotlight and its lighting parameters.
/// - Enables common OpenGL capabilities such as face culling, depth testing, and multisampling.
void mesh_view::initialize() {
    m_private->m_ucs_view->initialize();
    m_private->m_view->set_fov(dtr(fov));
    m_private->m_cam.reset(new gl_camera(fvec3(0, 0, 50), fvec3(0, 0, 0), fvec3(0, 1, 0)));

    m_private->m_shader.reset(new gl_shader);
    m_private->m_shader->add_file(GL_VERTEX_SHADER, "resources/shaders/mesh_tools_VertexShader.glsl");
    m_private->m_shader->add_file(GL_FRAGMENT_SHADER, "resources/shaders/mesh_tools_FragmentShader.glsl");
    m_private->m_shader->load();

    m_private->m_light.reset(new gl_light(gl_light::SPOTLIGHT));
    m_private->m_light->set_position(fvec3(-20, 20, 20));
    m_private->m_light->set_ambient(fvec3(0.75f));
    m_private->m_light->set_diffuse(fvec3(0.5f));
    m_private->m_light->set_specular(fvec3(0.1f));

    // OpenGL initialization
    glEnable(GL_CULL_FACE);
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LEQUAL);
    glEnable(GL_MULTISAMPLE);
}

/// @brief Renders the current scene.
///
/// Performs the full rendering pass:
/// - Configures the viewport to cover the entire window.
/// - Clears the color and depth buffers.
/// - Builds the combined camera matrix from the camera and viewport.
/// - Activates the shader and applies lighting and camera uniforms.
/// - Draws all mesh parts and active visualization overlays (normals and curvature).
/// - Renders the UCS axes widget on top.
void mesh_view::render() {
    // set the viewport to the whole window
    m_private->m_view->set_viewport();

    // clear screen
    glClearColor(.25f, .25f, .25f, 1.f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    fmat4 rot_mat(m_private->m_arcball->rotation());

    // combine the view and camera matrices into one
    // these matrices are COLUMN MAJOR!
    fmat4 cam_matrix = rot_mat * m_private->m_cam->perspective() * m_private->m_view->perspective();

    // enable the shader
    m_private->m_shader->use();
    m_private->m_light->apply(m_private->m_shader.get());
    // set the combined view matrix
    m_private->m_shader->set_mat4("camera", cam_matrix);

    // render the objects using color
    m_private->m_shader->set_vec4("object_color", fvec4(.9f, .9f, .9f, 1.f));

    for (const auto& part : m_private->m_draw_parts) {
        part->render(m_private->m_shader.get());
    }
    for (const auto& part : m_private->m_face_normals) {
        part->render(m_private->m_shader.get());
    }
    for (const auto& part : m_private->m_model_curvatures_min) {
        part->render(m_private->m_shader.get());
    }
    for (const auto& part : m_private->m_model_curvatures_max) {
        part->render(m_private->m_shader.get());
    }

    m_private->m_shader->end();

    // render coordinate system arrows
    if (m_private->m_ucs_view.get() != nullptr) {
        m_private->m_ucs_view->render(rot_mat);
    }
}

/// @brief Builds visualizations of vertex principal curvature directions for the current model.
///
/// For each mesh part of the loaded model:
/// - Uses the precomputed principal directions and the average edge length to derive a display scale.
/// - Creates two `visual_objects` instances: one for minimum curvature (red) and one for maximum curvature (green).
/// - For every vertex, adds line segments representing the principal directions.
/// - Creates corresponding `gl_prim` objects, aligns them with the base mesh rotation, sets colors,
///   and initializes them as initially hidden/visible (depending on `toggle_visibility()` default).
void mesh_view::generate_model_curvature_view() {
    m_private->m_model_curvatures_min.clear();
    m_private->m_model_curvatures_max.clear();
    if (m_private->m_model) {
        const auto& parts = m_private->m_model->get_parts();
        m_private->m_model_curvatures_min.reserve(parts.size());
        m_private->m_model_curvatures_max.reserve(parts.size());

        for (base_math::half_edge_mesh<float>* part : parts) {
            float vl = part->average_edge_length * 0.5f;
            auto mo_min = std::make_unique<visual_objects>();
            auto mo_max = std::make_unique<visual_objects>();
            // for each vertex, add curvature vectors
            for (const auto& v : part->get_vertices()) {
                fvec3 p = v.second->coords;
                fvec3 k1_end = p + v.second->principal_dir_min * vl;
                fvec3 k2_end = p + v.second->principal_dir_max * vl;
                // add min curvature vector in red
                mo_min->add_vector(p, k1_end);
                // add max curvature vector in green
                mo_max->add_vector(p, k2_end);
            }
            mo_min->create_prim();
            mo_min->get_prim()->rotate_to(m_private->m_draw_parts[0]->get_rotation());
            mo_min->get_prim()->set_color(fvec3(1.f, 0.f, 0.f));
            mo_min->toggle_visibility();
            mo_max->create_prim();
            mo_max->get_prim()->rotate_to(m_private->m_draw_parts[0]->get_rotation());
            mo_max->get_prim()->set_color(fvec3(0.f, 1.f, 0.f));
            mo_max->toggle_visibility();
            m_private->m_model_curvatures_min.push_back(std::move(mo_min));
            m_private->m_model_curvatures_max.push_back(std::move(mo_max));
        }
    }
}

/// @brief Computes per-vertex curvature information for a half-edge mesh.
///
/// External function used by `mesh_view::do_curvature_calculation()` to fill curvature data
/// (principal directions and magnitudes) for each vertex of the given mesh.
void compute_vertex_curvatures(half_edge_mesh<float>* mesh);

/// @brief Triggers curvature computation for the current model and updates its visualization.
///
/// Iterates over all mesh parts in the loaded model, computes per-vertex curvature using
/// `compute_vertex_curvatures()`, and then rebuilds the curvature visualization geometry
/// by calling `generate_model_curvature_view()`.
void mesh_view::do_curvature_calculation() {
    if (!m_private->m_model)
        return;
    const auto& parts = m_private->m_model->get_parts();
    for (base_math::half_edge_mesh<float>* part : parts) {
        compute_vertex_curvatures(part);
    }
    generate_model_curvature_view();
}

/// @brief Generates per-face normal visualization for the current model.
///
/// For each mesh part:
/// - Computes face normals and geometric properties.
/// - Creates a `visual_objects` instance storing a line segment per face originating from the face center
///   and pointing along the face normal, scaled by half the average edge length.
/// - Creates a blue `gl_prim` from the collected vectors and stores it for overlay rendering.
void mesh_view::generate_model_normals() {
    m_private->m_face_normals.clear();
    if (m_private->m_model) {
        const auto& parts = m_private->m_model->get_parts();
        m_private->m_face_normals.reserve(parts.size());

        for (base_math::half_edge_mesh<float>* part : parts) {
            float vl = part->average_edge_length * 0.5f;
            part->compute_face_normals();
            part->compute_face_properties();
            auto mo = std::make_unique<visual_objects>();
            for (const auto& face : part->faces) {
                fvec3 face_normal = face.second->normal;
                fvec3 face_center = face.second->center;
                fvec3 normal_end = face_center + face_normal * vl;
                mo->add_vector(face_center, normal_end);
            }
            mo->create_prim();
            mo->get_prim()->set_color(fvec3(0.f, 0.f, 1.f));
            m_private->m_face_normals.push_back(std::move(mo));
        }
    }
}

/// @brief Updates the viewport and UCS view when the window size changes.
///
/// @param width  New window width in pixels.
/// @param height New window height in pixels.
///
/// Adjusts the viewport aspect ratio and forwards the resize event to the UCS view (if present).
void mesh_view::resize(int width, int height) {
    m_private->m_arcball->resize((float)width, (float)height);
    m_private->m_view->set_window_aspect(width, height);
    if (m_private->m_ucs_view.get() != nullptr)
        m_private->m_ucs_view->resize_window(width, height);
}

/// @brief Loads a new mesh model from disk and prepares it for rendering.
///
/// @param fnm Path to the mesh file to load (e.g., OBJ).
/// @return `true` if the model was successfully loaded; otherwise `false`.
///
/// This method:
/// - Validates the file name pointer.
/// - Resets interactive state (blocks first mouse event, resets FOV, reorients UCS).
/// - Loads the model using `load_model()`.
/// - Converts each mesh part into a `gl_prim` for rendering.
/// - Clears any existing curvature overlays and regenerates face normal visualization.
bool mesh_view::load_new_model(const char* fnm) {
    if (fnm == nullptr) {
        return false;
    }

    block_mouse_event = true;
    fov = 15.f;
    m_private->m_ucs_view->rotate_ucs_to(0, 0, 0);
    m_private->m_model.reset(load_model(fnm));
    m_private->m_draw_parts.clear();

    if (m_private->m_model) {
        const auto& parts = m_private->m_model->get_parts();
        m_private->m_draw_parts.reserve(parts.size());

        for (const base_math::half_edge_mesh<float>* part : parts) {
            mesh_data md;
            collect_mesh_data(part, md);
            auto prim = std::make_unique<gl_prim>();
            prim->create_from_mesh(&md, GL_FILL);
            m_private->m_draw_parts.push_back(std::move(prim));
        }
    }
    m_private->m_model_curvatures_min.clear();
    m_private->m_model_curvatures_max.clear();
    generate_model_normals();
    return m_private->m_model != nullptr;
}

/// @brief Handles high-level mesh-related commands.
///
/// @param cmd Command identifier (e.g., menu/toolbar command).
/// @return Always returns 1 to indicate the command was processed.
///
/// Supported commands:
/// - `ID_FILE_LOAD`: Opens a file dialog to load a mesh file and initializes the view with it.
/// - `ID_MESH_HIDE_SHOW`: Toggles visibility of the face normal visualization overlay.
/// - `ID_MESH_CURVATURE`: Computes and displays curvature visualization for the current model.
int mesh_view::onCommand(int cmd) {
    switch (cmd) {
    case ID_FILE_LOAD:
    {
        const char* cfname = OpenFileDialog("All Files\0*.*\0Obj Files\0*.obj\0");
        if (cfname) {
            load_new_model(cfname);
        }
    }
    break;

    case ID_MESH_HIDE_SHOW:
        for (const auto& part : m_private->m_face_normals) {
            part->toggle_visibility();
        }
        break;

    case ID_MESH_CURVATURE:
        do_curvature_calculation();
        break;
    }
    return 1;
}


/// @brief Handles mouse wheel input to zoom in and out by changing the field of view.
///
/// @param delta     Wheel delta (positive for wheel up, negative for wheel down).
/// @param extra_btn Bitmask of mouse buttons (currently unused).
///
/// Adjusts the `fov` value within a clamped range [0.5°, 65°] and updates the viewport
/// projection accordingly.
void mesh_view::onMouseWheel(int delta, unsigned __int64 extra_btn) {
    if (delta > 0)
        fov -= 0.5f;
    if (delta < 0)
        fov += 0.5f;
    if (fov < 0.5f)
        fov = 0.5f;
    if (fov > 65.f)
        fov = 65.f;

    m_private->m_view->set_fov(dtr(fov));
}

/// @brief Handles mouse move events while interacting with the mesh view.
///
/// Updates the arcball controller with the current cursor position to continuously
/// update the view rotation during an active drag operation.
///
/// @param x     Current mouse X position in window/client coordinates.
/// @param y     Current mouse Y position in window/client coordinates.
/// @param extra Additional mouse state flags (currently unused).
void mesh_view::onMouseMove(int x, int y, unsigned __int64 extra) {
    m_private->m_arcball->drag(float(x), float(y));
}

/// @brief Handles left mouse button press to start an arcball rotation.
///
/// Captures the initial mouse position and notifies the arcball controller that
/// a drag operation has started, enabling interactive rotation of the mesh view.
///
/// @param x     Mouse X position at the time of button press.
/// @param y     Mouse Y position at the time of button press.
/// @param extra Additional mouse state flags (currently unused).
void mesh_view::onLMouseDown(int x, int y, unsigned __int64 extra) {
    m_private->m_arcball->beginDrag(float(x), float(y));
}

/// @brief Handles left mouse button release to end an arcball rotation.
///
/// Signals the arcball controller to finalize the current drag operation, freezing
/// the current view rotation until a new drag is started.
///
/// @param x     Mouse X position at the time of button release (unused).
/// @param y     Mouse Y position at the time of button release (unused).
/// @param extra Additional mouse state flags (currently unused).
void mesh_view::onLMouseUp(int x, int y, unsigned __int64 extra) {
    m_private->m_arcball->endDrag();
}