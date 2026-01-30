/// @file visual_objects.cpp
/// @brief Implements simple visual helper primitives and a User Coordinate System (UCS) view.
//
/// This file provides two main utilities:
/// - `visual_objects`: a lightweight container used to accumulate line segments into an
///   internal `gl_mesh` and build a corresponding `gl_prim` for rendering (typically used
///   for debug or helper visuals such as vectors).
/// - `UCS_view`: a small, self-contained viewport that renders a user coordinate system
///   widget using its own camera, light, shader, and primitive setup.

#include "visual_objects.h"
#include <memory>

#include "gl_graphics.h"
#include "gl_camera.h"
#include "gl_light.h"
#include "gl_mesh.h"
#include "gl_prim.h"

using namespace base_math;
using namespace base_opengl;

/// @class visual_objects
/// @brief Aggregates simple visual line primitives and exposes them as a drawable `gl_prim`.
///
/// `visual_objects` owns an internal `gl_mesh` that accumulates vertices and indices. Each
/// call to `add_vector` adds one line segment (two vertices and two indices). Once all
/// desired vectors have been added, `create_prim` turns the mesh into a line-based
/// `gl_prim` that can be rendered using a given `gl_shader`.
///
/// Typical usage:
/// @code
/// visual_objects vis;
/// vis.add_vector(fvec3(0, 0, 0), fvec3(1, 0, 0));
/// vis.add_vector(fvec3(0, 0, 0), fvec3(0, 1, 0));
/// gl_prim* prim = vis.create_prim();
/// vis.m_visible = true;
/// vis.render(shader);
/// @endcode
visual_objects::visual_objects() {
    // Allocate the backing mesh used to store line vertices and indices.
    m_mesh = new gl_mesh;
    // No primitive exists until `create_prim` is called.
    m_prim = nullptr;
    // Visuals are initially not rendered; caller must enable visibility explicitly.
    m_visible = false;
}

visual_objects::~visual_objects() {
    // Release owned rendering resources.
    if (m_mesh) delete m_mesh;
    if (m_prim) delete m_prim;
}

/// @brief Adds a single 3D vector as a line segment to the internal mesh.
/// @param vstart World-space start position of the vector.
/// @param vend   World-space end position of the vector.
/// @return Always returns 1 on success.
///
/// This function appends two vertices (`vstart`, `vend`) to the internal `gl_mesh` and
/// adds an index pair referencing those vertices as a line segment. The mesh is only
/// converted to a renderable primitive when `create_prim` is invoked.
int visual_objects::add_vector(const fvec3& vstart, const fvec3& vend) {
    unsigned int idx_start = (unsigned int)m_mesh->n_vertices();
    m_mesh->addVertex(vstart);
    m_mesh->addVertex(vend);
    m_mesh->addIndices(idx_start, idx_start + 1);
    return 1;
}

/// @brief Creates (or recreates) the drawable primitive from the accumulated mesh data.
/// @return Pointer to the internally owned `gl_prim` configured for line rendering.
///
/// The current contents of `m_mesh` are collected into a `mesh_data` structure and used
/// to build a `gl_prim` with `GL_LINES` draw mode. Subsequent calls will overwrite the
/// previous primitive instance. Call this after finishing all `add_vector` calls and
/// before invoking `render`.
gl_prim* visual_objects::create_prim() {
    mesh_data mdata;
    collect_mesh_data(m_mesh, mdata);
    m_prim = new gl_prim;
    m_prim->set_draw_type(GL_LINES);
    m_prim->create_from_mesh(&mdata, GL_LINE);
    return m_prim;
}

/// @brief Renders the internal primitive if it exists and is marked visible.
/// @param sh Active shader used to render the line primitive.
///
/// The caller is responsible for configuring the shader (e.g., view/projection matrices,
/// colors) before invoking this function. If `m_visible` is `false` or `m_prim` is null,
/// this function does nothing.
void visual_objects::render(gl_shader* sh) {
    if (m_prim && m_visible) {
        m_prim->render(sh);
    }
}

//////////////////////////////////////////////////////////////

/// @brief Private implementation details for `UCS_view`.
///
/// This structure encapsulates the resources required by `UCS_view`, including:
/// - `m_view`: viewport configuration (projection, size, position).
/// - `m_cam`: camera controlling the UCS view orientation and position.
/// - `m_shader`: shader program used to render the UCS primitive.
/// - `m_light`: simple spotlight used to illuminate the UCS.
/// - `m_ucs`: the actual UCS `gl_prim` (e.g., axis triad).
struct UCS_view_private {
    std::unique_ptr<gl_viewport> m_view;             ///< Viewport configuration
    std::unique_ptr<gl_camera> m_cam;                ///< Camera for scene viewing
    std::unique_ptr<gl_shader> m_shader;             ///< Shader program for rendering

    std::unique_ptr<gl_light> m_light;               ///< Dedicated light for the UCS widget
    std::unique_ptr<gl_prim> m_ucs;                  ///< User coordinate system visualization
};

/// @class UCS_view
/// @brief Small, self-contained view that renders a user coordinate system widget.
///
/// `UCS_view` manages its own viewport, camera, shader, light, and UCS primitive. It is
/// typically rendered into a fixed corner of the main window as an orientation helper
/// (similar to axis widgets in CAD/3D tools). The user can rotate the UCS via
/// `rotate_ucs_by` or `rotate_ucs_to`.
UCS_view::UCS_view() {
    // Allocate the private implementation and initialize the viewport container.
    m_private_data = new UCS_view_private;
    m_private_data->m_view.reset(new gl_viewport());
}

UCS_view::~UCS_view() {
    // Destroy the PIMPL container; smart pointers inside clean up their resources.
    delete m_private_data;
}

/// @brief Initializes the UCS view resources (projection, camera, light, shader, and UCS).
///
/// This must be called once before any calls to `render()`. It:
/// - Configures a perspective projection on the internal viewport.
/// - Creates a camera positioned along +Z looking at the origin.
/// - Sets up a white spotlight illuminating the UCS from above/side.
/// - Loads the mesh tools vertex/fragment shaders.
/// - Creates and positions the UCS primitive at the world origin with a default scale.
void UCS_view::initialize() {
    m_private_data->m_view->set_perspective(dtr<float>(30.f), 0.1f, 10000.f);
    m_private_data->m_cam.reset(new gl_camera(fvec3(0, 0, 30), fvec3(0, 0, 0), fvec3(0, 1, 0)));

    m_private_data->m_light.reset(new gl_light(gl_light::SPOTLIGHT));
    m_private_data->m_light->set_position(fvec3(-3000, 3000, 3000));
    m_private_data->m_light->set_ambient(fvec3(1, 1, 1));
    m_private_data->m_light->set_diffuse(fvec3(1, 1, 1));
    m_private_data->m_light->set_specular(fvec3(1, 1, 1));

    m_private_data->m_shader.reset(new gl_shader);
    m_private_data->m_shader->add_file(GL_VERTEX_SHADER, "resources/shaders/mesh_tools_VertexShader.glsl");
    m_private_data->m_shader->add_file(GL_FRAGMENT_SHADER, "resources/shaders/mesh_tools_FragmentShader.glsl");
    m_private_data->m_shader->load();

    m_private_data->m_ucs.reset(create_UCS());
    m_private_data->m_ucs->move_to(0, 0, 0);
    m_private_data->m_ucs->rotate_to(0, 0, 0);
    m_private_data->m_ucs->set_scale(fvec3(6.f));
}

/// @brief Renders the UCS widget into its configured viewport.
///
/// This sets the viewport, binds the UCS shader, applies the UCS-specific light, uploads
/// the combined camera and projection matrices, and finally draws the UCS primitive. The
/// shader uniform `object_or_vertex_color` is set to 0 to indicate that object-level
/// color should be used instead of per-vertex color.
void UCS_view::render(const fmat4& user_rotation) {
    m_private_data->m_view->set_viewport();
    m_private_data->m_shader->use();
    m_private_data->m_light->apply(m_private_data->m_shader.get());
    fmat4 cam_matrix = user_rotation*m_private_data->m_cam->perspective() * m_private_data->m_view->perspective();
    m_private_data->m_shader->set_mat4("camera", cam_matrix);
    m_private_data->m_shader->set_vec3("cameraPos", m_private_data->m_cam->vLocation);
    m_private_data->m_shader->set_int("object_or_vertex_color", 0);  // use object color
    m_private_data->m_ucs->render(m_private_data->m_shader.get());
    m_private_data->m_shader->end();
}

/// @brief Updates the UCS viewport placement and aspect ratio.
/// @param width  New window width (currently ignored).
/// @param height New window height (currently ignored).
///
/// Currently this implementation pins the UCS viewport to the origin (0,0) and forces a
/// fixed aspect ratio of 1:1 (100x100). The `width` and `height` parameters are not yet
/// used, but are kept for a future implementation that may position and scale the UCS
/// relative to the main window size.
void UCS_view::resize_window(int width, int height) {
    m_private_data->m_view->set_position(0, 0);
    m_private_data->m_view->set_window_aspect(100, 100);
}

/// @brief Rotates the UCS widget incrementally around its local axes.
/// @param x Delta rotation in degrees around the X axis.
/// @param y Delta rotation in degrees around the Y axis.
/// @param z Delta rotation in degrees around the Z axis.
///
/// This function applies a relative rotation to the existing UCS orientation, allowing
/// interactive rotation controls (e.g., mouse drag) to adjust the widget.
void UCS_view::rotate_ucs_by(float x, float y, float z) {
    m_private_data->m_ucs->rotate_by(fvec3(x, y, z));
}

/// @brief Sets the absolute orientation of the UCS widget.
/// @param x Absolute rotation in degrees around the X axis.
/// @param y Absolute rotation in degrees around the Y axis.
/// @param z Absolute rotation in degrees around the Z axis.
///
/// Unlike `rotate_ucs_by`, this function overwrites the current orientation of the UCS.
/// It is useful for snapping the widget to a specific orientation (e.g., reset or
/// predefined views).
void UCS_view::rotate_ucs_to(float x, float y, float z) {
    m_private_data->m_ucs->rotate_to(x, y, z);// fvec3(x, y, z));
}