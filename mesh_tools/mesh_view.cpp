/// @file mesh_view.cpp
/// @brief Implements the `mesh_view` class responsible for displaying and interacting with 3D mesh models.
///
/// This module sets up an OpenGL-based viewport, camera, lighting, and shader pipeline to render
/// mesh models. It also provides visualization for face normals and principal curvature directions,
/// and exposes user interaction through mouse input and command handling.
#define _CRT_SECURE_NO_WARNINGS
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

#undef max // avoid conflicts with std::max in some Windows headers

using namespace base_math;
using namespace base_opengl;

/// external function to compute vertex curvatures, implemented in mesh_tools/curvature.cpp
void compute_vertex_curvatures(half_edge_mesh<float>* mesh);


/// @brief Internal data container for `mesh_view`.
///
/// This struct encapsulates all OpenGL-related objects and visualization helpers used by `mesh_view`.
/// It is allocated on the heap and owned by `mesh_view` to hide implementation details from clients.
struct mesh_view_private {
    std::unique_ptr<gl_camera> m_cam;                      ///< gl_camera used to view the scene and compute view/projection.
    std::unique_ptr<gl_shader> m_shader;                ///< Shader program used for mesh and helper rendering.
    std::unique_ptr<gl_light> m_light;                  ///< Primary scene light affecting shading.

    std::unique_ptr<UCS_view> m_ucs_view;               ///< UCS (user coordinate system) view used to render the axes widget.

    // {{ the loaded model and the acompanying visualization objects
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
    m_private->m_ucs_view.reset(new UCS_view());
    m_private->m_arcball.reset(new arcball(800, 600));
    m_private->m_cam.reset(new gl_camera(fvec3(0, 0, 50), fvec3(0, 0, 0), fvec3(0, 1, 0)));

    m_model_state.m_curvatures_calculated = false;
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
    m_private->m_cam->set_fov(dtr(30.f));

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
    //glEnable(GL_CULL_FACE);
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LEQUAL);
    glEnable(GL_MULTISAMPLE);
}


/// @brief Renders all mesh geometry and related overlays for a given camera and rotation.
/// 
/// This helper encapsulates the actual scene drawing logic used by `mesh_view::render()`.
/// It assumes that the OpenGL viewport and frame buffers have already been configured.
///
/// Rendering steps:
/// - Binds and activates the provided shader program.
/// - Applies the active light source to the shader.
/// - Uploads the combined projection * view matrix as the `camera` uniform.
/// - Sets the default `object_color` uniform used for shaded mesh rendering.
/// - Draws all mesh parts twice:
///   - First pass in filled polygon mode with depth testing enabled.
///   - Second pass in line (wireframe) mode, using polygon offset so the wireframe
///     is rendered slightly in front of the filled geometry and forcing black color.
/// - Draws all auxiliary visualizations (face normals and principal curvature
///   directions) using the same rotation matrix so they stay aligned with the mesh.
///
/// @param cam_matrix Combined projection * view matrix of the current camera.
/// @param rot_mat    Rotation matrix derived from the arcball, applied to all
///                   mesh geometry and visualization primitives.
/// @param shdr       Active shader used for rendering the mesh parts; must be a
///                   valid, already linked `gl_shader` instance.
void mesh_view::render_scene(fmat4& cam_matrix, fmat4& rot_mat, gl_shader* shdr)
{
    shdr->use();
    m_private->m_light->apply(shdr);
    // set the combined view matrix
    shdr->set_mat4("camera", cam_matrix);

    // render the objects using color
    shdr->set_vec4("object_color", fvec4(.9f, .9f, .9f, 1.f));

    // render the mesh parts with the current rotation applied
    {
        // render filled polygons first
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        glEnable(GL_DEPTH_TEST);
        for (const auto& part : m_private->m_draw_parts) {
            part->set_draw_mode(GL_FILL);
            part->force_black = false;
            part->view_matrix = rot_mat;   // apply the current arcball rotation to the mesh parts
            part->set_use_vertex_color(m_view_state.show_gauss_map ? 1 : 0); // ensure vertex color is disabled by default
            part->render(shdr);
        }
        if (m_view_state.show_wireframe) {
            // then render wireframe on top
            glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
            glEnable(GL_POLYGON_OFFSET_LINE);
            glPolygonOffset(-1.0f, -1.0f); // pull lines toward camera
            for (const auto& part : m_private->m_draw_parts) {
                part->set_draw_mode(GL_LINE);
                part->force_black = true; // render wireframe in black
                part->render(shdr);
            }
        }
    }

    if (m_view_state.show_normals) {
        for (const auto& part : m_private->m_face_normals) {
            part->get_prim()->view_matrix = rot_mat;  // apply the same rotation to the normals visualization
            part->render(m_private->m_shader.get());
        }
    }
    if (m_view_state.show_curvature) {
        for (const auto& part : m_private->m_model_curvatures_min) {
            part->get_prim()->view_matrix = rot_mat;  // apply the same rotation to the curvature visualization
            part->render(m_private->m_shader.get());
        }
        for (const auto& part : m_private->m_model_curvatures_max) {
            part->get_prim()->view_matrix = rot_mat;  // apply the same rotation to the curvature visualization
            part->render(m_private->m_shader.get());
        }
    }

    m_private->m_shader->end();
}

/// @brief Renders the complete mesh view each frame.
///
/// Sets up the OpenGL viewport and clears the color/depth buffers, then renders
/// the main scene (mesh, overlays, lighting) followed by the auxiliary user
/// coordinate system (UCS) widget.
///
// Rendering steps:
/// - Applies the camera's viewport via `gl_camera::set_viewport()`.
/// - Clears the frame using a dark gray background and clears the depth buffer.
/// - Computes the projection * view matrix from the active camera.
/// - Retrieves the current rotation from the arcball controller as a rotation matrix.
/// - Invokes `render_scene()` to draw the mesh, wireframe, normals, and curvature
///   visualizations using the current camera and rotation.
/// - If a UCS view is available, applies the same rotation so its axes match the
///   mesh orientation and renders the UCS overlay in screen space.
///
/// This method is intended to be called once per frame from the main render loop.
void mesh_view::render() {
    m_private->m_cam->set_viewport();
    glClearColor(.25f, .25f, .25f, 1.f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    fmat4 cam_matrix = m_private->m_cam->perspective();
    fmat4 rot_mat(m_private->m_arcball->rotation());

    render_scene(cam_matrix, rot_mat, m_private->m_shader.get());

    // render coordinate system arrows
    if (m_private->m_ucs_view.get() != nullptr) {
        m_private->m_ucs_view->set_user_rotation(rot_mat); // apply the same rotation to the UCS view
        m_private->m_ucs_view->render();
    }
}

/// @brief Builds or refreshes visualization objects for per-vertex principal curvatures of the current model.
///
/// Clears any previously generated curvature visualization and, for each mesh part of the
/// currently loaded model, creates two `visual_objects` collections:
/// - one for minimum principal curvature directions (rendered in red),
/// - one for maximum principal curvature directions (rendered in green).
///
/// For every vertex in a part:
/// - The base position `p` is taken from the vertex coordinates.
/// - Two end points are computed by adding scaled principal direction vectors
///   (`principal_directions[0]` and `principal_directions[1]`) to `p`.
/// - Line segments representing these directions are added to the corresponding
///   `visual_objects` instances.
///
/// The vectors are scaled by half of the part's `average_edge_length` so that
/// they are visually proportional to the mesh size. After all vectors are added,
/// a `gl_prim` is created for each `visual_objects`, aligned to the rotation of
/// the first mesh draw part (`m_draw_parts[0]`), and colored appropriately:
/// - minimum curvature vectors in red (`(1, 0, 0)`),
/// - maximum curvature vectors in green (`(0, 1, 0)`).
///
/// Both visualization sets are initially toggled (typically making them visible),
/// and stored in `m_model_curvatures_min` / `m_model_curvatures_max` so that
/// they can be rendered alongside the mesh in `render_scene()`.
void mesh_view::generate_model_curvature_view() {
    m_view_state.show_curvature = false; // ensure curvature visualization is disbled while we rebuild it

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
                fvec3 k1_end = p + v.second->curvature_data.principal_directions[0] * vl;
                fvec3 k2_end = p + v.second->curvature_data.principal_directions[1] * vl;
                // add min curvature vector in red
                mo_min->add_vector(p, k1_end);
                // add max curvature vector in green
                mo_max->add_vector(p, k2_end);
            }
            mo_min->create_prim();
            mo_min->get_prim()->rotate_to(m_private->m_draw_parts[0]->get_rotation());
            mo_min->get_prim()->set_color(fvec3(1.f, 0.f, 0.f));
            mo_max->create_prim();
            mo_max->get_prim()->rotate_to(m_private->m_draw_parts[0]->get_rotation());
            mo_max->get_prim()->set_color(fvec3(0.f, 1.f, 0.f));
            m_private->m_model_curvatures_min.push_back(std::move(mo_min));
            m_private->m_model_curvatures_max.push_back(std::move(mo_max));
        }
    }
}

/// @brief Generates or refreshes the Gauss curvature color map for the current model.
///
/// This method prepares per-vertex Gauss curvature for visualization by:
/// - Checking that a model is loaded and that the Gauss curvature view has not
///   already been generated for the current curvature state.
/// - Normalizing the absolute Gauss curvature values of all vertices into the
///   [0, 1] range via `model::normalize_absolute_curvature()`. The normalized
///   values are typically stored per vertex and later interpreted as colors
///   (e.g. for a heat map) in the rendering pipeline.
/// - Rebuilding the mesh GPU representation by calling `build_model_representation()`
///   so that any curvature-dependent vertex colors are uploaded to the graphics
///   primitives used for rendering.
///
/// After successful generation, the flag `m_gauss_curvature_view_generated` is
/// set to `true` to avoid redundant recomputation until curvature data changes.
void mesh_view::generate_model_gauss_curvature_view() {
    if (m_private->m_model && !m_model_state.m_gauss_curvature_view_generated) {
        // normalize absolut Gauss curvature values to [0, 1] range for visualization
        m_private->m_model->normalize_absolute_curvature();
        // rebuild the mesh representation to update vertex colors based on the new curvature values
        build_model_representation(); 
        m_model_state.m_gauss_curvature_view_generated = true;
    }
}

/// @brief Triggers curvature computation for the current model and updates its visualization.
///
/// Invokes curvature analysis on each mesh part of the currently loaded model by calling
/// `compute_vertex_curvatures` for every `half_edge_mesh` instance. This function is expected
/// to populate the `curvature_data` for each vertex (principal directions, curvature magnitudes, etc.).
///
/// Once curvature data is available, `generate_model_curvature_view()` is called to rebuild
/// the OpenGL visualization objects used to display the principal curvature directions on screen.
///
/// This function is typically called from command handling (e.g. menu/toolbar actions) to
/// enable or refresh curvature visualization after loading or modifying a mesh.
void mesh_view::do_curvature_calculation() {
    if (!m_private->m_model)
        return;
    const auto& parts = m_private->m_model->get_parts();
    for (base_math::half_edge_mesh<float>* part : parts) {
        compute_vertex_curvatures(part);
    }
    generate_model_curvature_view();
    m_model_state.m_curvatures_calculated = true;
    m_model_state.m_gauss_curvature_view_generated = false;
}

/// @brief Generates per-face normal visualization for the current model.
///
/// Clears any existing face-normal visualizations and, for each mesh part of the
/// current model, computes and visualizes face normals as line segments.
///
/// For each part:
/// - `compute_face_normals()` updates per-face normal vectors.
/// - `compute_face_properties()` updates derived quantities such as face centers.
/// - A `visual_objects` instance is created.
/// - For every face in the part, a line segment is added that starts at the face
///   center and ends at `center + normal * vl`, where `vl` is half the part's
///   `average_edge_length`, ensuring the normal length is scaled to the mesh size.
/// - A `gl_prim` is created from the accumulated vectors and colored blue
///   (`(0, 0, 1)`), then stored in `m_face_normals`.
///
/// The generated primitives are used by `render_scene()` to overlay face normals
/// on top of the shaded mesh, aiding debugging and inspection of mesh quality.
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

/// @brief Updates camera, arcball, and UCS view when the host window is resized.
///
/// Adjusts the interactive controls and projection parameters to match a new
/// viewport size:
/// - Resizes the arcball controller so mouse interactions remain consistent
///   with the new window dimensions.
/// - Updates the camera aspect ratio used for perspective projection.
/// - If a UCS (user coordinate system) view exists, notifies it about the new
///   window size so its overlay rendering stays aligned with the main viewport.
///
/// @param width  New client-area width in pixels.
/// @param height New client-area height in pixels.
void mesh_view::resize(int width, int height) {
    m_private->m_arcball->resize((float)width, (float)height);
    m_private->m_cam->set_aspect(width, height);
    if (m_private->m_ucs_view.get() != nullptr)
        m_private->m_ucs_view->resize_window(width, height);
}

/// @brief Loads a new mesh model from disk and initializes its rendering resources.
///
/// Replaces the currently loaded model, resets view-related state, and rebuilds
/// all associated OpenGL primitives and visualizations.
///
/// Behavior:
/// - Returns `false` immediately if `fnm` is `nullptr`.
/// - Sets `block_mouse_event` flag to temporarily block mouse interaction while
///   the new model is being loaded.
/// - Resets the field of view to a default (15 degrees) and applies it to the camera.
/// - Resets the arcball rotation to the identity orientation.
/// - Calls `load_model(fnm)` to create a new `model` instance and assigns it to
///   `m_model`.
/// - Clears existing `m_draw_parts` and, for each mesh part in the model:
///   - Collects raw vertex/index data via `collect_mesh_data`.
///   - Creates a `gl_prim`, initializes it with the mesh data in `GL_FILL` mode,
///     and stores it in `m_draw_parts`.
/// - Clears any previously generated curvature visualizations.
/// - Regenerates per-face normals visualization via `generate_model_normals()`.
///
/// @param fnm Null-terminated path to the mesh file to load (e.g. OBJ).
/// @return `true` if a model was successfully loaded and assigned, `false` otherwise.
bool mesh_view::load_new_model(const char* fnm) {
    if (fnm == nullptr) {
        return false;
    }

    block_mouse_event = true;
    
    // reset view parameters
    m_private->m_cam->setup(fvec3(0, 0, 50), fvec3(0, 0, 0), fvec3(0, 1, 0));

    // reset arcball rotation
    m_private->m_arcball->reset();

    m_model_state.m_curvatures_calculated = false;

    m_private->m_model.reset(load_model(fnm));

    m_view_state.set_model_changed();  // disable cached view states that depend on the model (e.g. curvature visualization) until they are regenerated

    build_model_representation();

    m_private->m_model_curvatures_min.clear();
    m_private->m_model_curvatures_max.clear();
    generate_model_normals();
    return m_private->m_model != nullptr;
}

/// @brief Rebuilds the GPU-backed representation of the currently loaded model.
///
/// Clears and recreates the collection of drawable OpenGL primitives (`m_draw_parts`)
/// from the logical mesh model stored in `m_model`.
///
/// Behavior:
/// - Empties `m_draw_parts` so that any previous mesh geometry is discarded.
/// - If a model is present:
///   - Retrieves all `half_edge_mesh` parts via `model::get_parts()`.
///   - Reserves space in `m_draw_parts` to avoid repeated reallocations.
///   - For each mesh part:
///     - Populates a local `mesh_data` instance using `collect_mesh_data()`,
///       which extracts vertex positions, indices, normals, and any auxiliary
///       attributes required for rendering.
///     - Allocates a new `gl_prim`, initializes it from the collected mesh data
///       by calling `gl_prim::create_from_mesh(&md, GL_FILL)`, and stores the
///       resulting primitive in `m_draw_parts`.
///
/// The resulting primitives are used by the rendering pipeline (e.g. in
/// `render_scene()`) to draw the mesh in filled-polygon mode, and may later be
/// reused for additional passes (such as wireframe overlay or curvature coloring).
void mesh_view::build_model_representation() {
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

}

/// @brief Handles high-level mesh view commands dispatched from the application UI.
///
/// Interprets and executes menu/toolbar command identifiers that affect the
/// mesh view, such as loading a model, toggling normal visualization, or
/// computing curvature.
///
/// Supported commands:
/// - `ID_FILE_LOAD`:
///   - Opens a file dialog allowing the user to select a mesh file (e.g. OBJ).
///   - If a file is chosen, calls `load_new_model()` with the selected path.
/// - `ID_MESH_HIDE_SHOW`:
///   - Toggles the visibility of all face-normal visualization objects stored
///     in `m_face_normals`.
/// - `ID_MESH_CURVATURE`:
///   - Invokes `do_curvature_calculation()` to compute and display curvature
///     directions for the current model.
///
/// @param cmd Application-specific command identifier (e.g. from a menu or accelerator).
/// @return Always returns 1 to indicate the command was handled.
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

    case ID_VIEW_FACENORMALS:
        m_view_state.show_normals = !m_view_state.show_normals;
        break;

    case ID_VIEW_PRINCIPALCURVATURES:
        if (m_model_state.m_curvatures_calculated == false) {
            MessageBox(NULL, "Curvature must be calculated before generating the curvature view.", "Error", MB_ICONERROR);
            return 1;
        }
        m_view_state.show_curvature = !m_view_state.show_curvature;
        break;

    case ID_VIEW_MESH:
        m_view_state.show_wireframe = !m_view_state.show_wireframe;
        break;

    case ID_VIEW_GAUSSCURVATURE:
        if (m_model_state.m_curvatures_calculated == false) {
            MessageBox(NULL, "Curvature must be calculated before generating the Gauss curvature view.", "Error", MB_ICONERROR);
            return 1;
        }
        generate_model_gauss_curvature_view();
        m_view_state.show_gauss_map = !m_view_state.show_gauss_map;

        break;        

    case ID_MESH_CURVATURE:
        do_curvature_calculation();
        break;
    }
    return 1;
}

/// @brief Handles mouse wheel input to zoom the camera in and out.
/// 
/// Translates mouse wheel delta steps into a floating-point value and forwards
/// it to the underlying `gl_camera` via `gl_camera::zoom`. Positive deltas
/// typically move the camera closer to the scene (zoom in), while negative deltas
/// move it farther away (zoom out), depending on the `gl_camera` implementation.
///
/// @param delta      Mouse wheel delta value provided by the windowing system.
/// @param extra_btn  Additional mouse state flags (currently unused).
void mesh_view::onMouseWheel(int delta, unsigned __int64 extra_btn) {
    m_private->m_cam->zoom(float(delta));
}

/// @brief Handles mouse move events to rotate and pan the view.
///
/// Always forwards the current mouse position to the arcball controller via
/// `arcball::drag` to update the mesh rotation while a left-button drag is active.
///
/// Additionally, when `dragging` is `true` (set on right-button drag):
/// - Computes the mouse delta relative to the last recorded position.
/// - Calls `gl_camera::pan` to translate the camera parallel to the view plane,
///   enabling panning of the scene.
/// - Updates `last_mouse_x` / `last_mouse_y` with the current position so that
///   subsequent moves compute incremental deltas.
///
/// @param x     Current mouse X position in window coordinates.
/// @param y     Current mouse Y position in window coordinates.
/// @param extra Additional mouse state (currently unused).
void mesh_view::onMouseMove(int x, int y, unsigned __int64 extra) {
    m_private->m_arcball->drag(float(x), float(y));
    if (dragging) {
        int deltax = x - last_mouse_x;
        int deltay = y - last_mouse_y;

        m_private->m_cam->pan(float(deltax), float(deltay)); // adjust the camera position based on mouse movement

        last_mouse_x = x;
        last_mouse_y = y;
    }
}

/// @brief Handles left mouse button press to begin arcball rotation.
///
/// Starts an arcball drag operation from the given mouse position, enabling
/// interactive rotation of the mesh around its center during subsequent mouse moves.
///
/// @param x     Mouse X position at the time of button press.
/// @param y     Mouse Y position at the time of button press.
/// @param extra Additional mouse state (currently unused).
void mesh_view::onLMouseDown(int x, int y, unsigned __int64 extra) {
    m_private->m_arcball->beginDrag(float(x), float(y));
}

/// @brief Handles left mouse button release to end arcball rotation.
///
/// Signals the arcball controller that the drag operation has finished, freezing
/// the current rotation until the next drag is started.
///
/// @param x     Mouse X position at the time of button release.
/// @param y     Mouse Y position at the time of button release.
/// @param extra Additional mouse state (currently unused).
void mesh_view::onLMouseUp(int x, int y, unsigned __int64 extra) {
    m_private->m_arcball->endDrag();
}

/// @brief Handles right mouse button press to start camera panning mode.
///
/// Enables panning by:
/// - Setting the `dragging` flag to `true`.
/// - Recording the current mouse coordinates as the reference for subsequent
///   move events handled by `onMouseMove()`.
///
/// While `dragging` is `true`, mouse movements will be converted into camera
/// panning offsets.
///
/// @param x     Mouse X position at the time of button press.
/// @param y     Mouse Y position at the time of button press.
/// @param extra Additional mouse state (currently unused).
void mesh_view::onRMouseDown(int x, int y, unsigned __int64 extra){
    dragging = true;
    last_mouse_x = x;
    last_mouse_y = y;
}

/// @brief Handles right mouse button release to stop camera panning mode.
///
/// Disables panning by clearing the `dragging` flag. Subsequent mouse move
/// events will no longer translate the camera until panning is re-enabled
/// via `onRMouseDown()`.
///
/// @param x     Mouse X position at the time of button release.
/// @param y     Mouse Y position at the time of button release.
/// @param extra Additional mouse state (currently unused).
void mesh_view::onRMouseUp(int x, int y, unsigned __int64 extra) {
    dragging = false;
}
