#include <iostream>
#include "themeshview.h"
#include "gl_context.h"
#include "glew.h"
#include "application.h"
#include "cmd_target.h"

#include "camera.h"
#include "shaders.h"
#include "light.h"
#include "prim.h"
#include "arcball.h"
#include "ucs_view.h"
#include "font.h"

#include "themeshmodel.h"
#include "curvature.h"


static inline fvec3 to_fvec3(const dvec3& v) {
    return fvec3(static_cast<float>(v.x()), static_cast<float>(v.y()), static_cast<float>(v.z()));
}


class visual_object : public gl_prim {
    mesh_data m_mesh_data;
public:
    visual_object() {
        draw_type = GL_LINES;
        draw_mode = GL_LINE;
        draw_elements = true;
    }
    void add_vector(const btm::fvec3& vstart, const btm::fvec3& vend) {
        unsigned int idx_start = (unsigned int)m_mesh_data.num_vertices/3;
        m_mesh_data.add_vertex(vstart);
        m_mesh_data.add_vertex(vend);
        m_mesh_data.add_indices(idx_start, idx_start + 1);
    }
    void create_prim() {
        create_from_mesh(&m_mesh_data, GL_LINE);
    }
    void clear_mesh_data() {
        m_mesh_data.vertices.clear();
        m_mesh_data.normals.clear();
        m_mesh_data.indices.clear();
        m_mesh_data.num_vertices = 0;
        m_mesh_data.num_normals = 0;
        m_mesh_data.num_indices = 0;
    }
};

struct mesh_view_resources {
    std::unique_ptr<gl_camera> m_cam;                   ///< gl_camera used to view the scene and compute view/projection.
    std::unique_ptr<gl_shader> m_shader;                ///< Shader program used for mesh and helper rendering.
    std::unique_ptr<gl_light> m_light;                  ///< Primary scene light affecting shading.

    std::unique_ptr<UCS_view> m_ucs_view;               ///< UCS (user coordinate system) view used to render the axes widget.

    std::unique_ptr<arcball> m_arcball;                                 ///< Arcball for mouse interaction
    gl_font* font2D;

    // {{ the loaded model and the acompanying visualization objects
    std::vector<std::unique_ptr<gl_prim>> m_draw_parts; ///< Drawable mesh parts converted to OpenGL primitives.
    // }}

    std::unique_ptr<visual_object> m_face_normals;
    std::unique_ptr<visual_object> m_vertex_normals;
    std::unique_ptr<visual_object> m_model_curvatures_k1;                  ///< Visualization of minimum principal curvature directions per vertex.
    std::unique_ptr<visual_object> m_model_curvatures_k2;                  ///< Visualization of maximum principal curvature directions per vertex.

    // bool vertex_normals_created = false;                                ///< Flag to track whether vertex normals visualization has been generated.
    // std::vector<std::unique_ptr<visual_objects>> m_vertex_normals;      ///< Visualization of per-vertex normals as line segments.
    // bool curvature_directions_created = false;                          ///< Flag to track whether curvature directions visualization has been generated.
    // std::vector<std::unique_ptr<visual_object>> m_model_curvatures_k2; ///< Visualization of maximum principal curvature directions per vertex.

    void invalidate_visualizations() {
        if (m_face_normals) {
            m_face_normals->clear_vao();
            m_face_normals->clear_mesh_data();
        }
        if (m_vertex_normals) {
            m_vertex_normals->clear_vao();
            m_vertex_normals->clear_mesh_data();
        }
        if (m_model_curvatures_k1) {
            m_model_curvatures_k1->clear_vao();
            m_model_curvatures_k1->clear_mesh_data();
        }
        if (m_model_curvatures_k2) {
            m_model_curvatures_k2->clear_vao();
            m_model_curvatures_k2->clear_mesh_data();
        }
        // face_normals_created = false;
        // vertex_normals_created = false;
        // curvature_directions_created = false;
    }

    void init() {
        m_cam = std::make_unique<gl_camera>(fvec3(0, 0, 50), fvec3(0, 0, 0), fvec3(0, 1, 0));

        m_shader = std::make_unique<gl_shader>();
        m_shader->add_file(GL_VERTEX_SHADER, "resources/shaders/meshprojectVertexShader.glsl");
        m_shader->add_file(GL_FRAGMENT_SHADER, "resources/shaders/meshprojectFragmentShader.glsl");
        m_shader->load();

        m_light = std::make_unique<gl_light>();
        m_light->set_position(fvec3(-20, 20, 50));

        m_ucs_view = std::make_unique<UCS_view>();
        m_ucs_view->initialize();
        m_arcball = std::make_unique<arcball>(800.f, 600.f);

        font2D = get_font_manager().create_font("Consolas", "Consolas", 12);
    }
};
mesh_view_resources* g_view_resources = nullptr;

theMeshView::theMeshView(cModel* model/*=nullptr*/) : m_model(model) {
    init();
}

void theMeshView::set_callbacks() {
    REGISTER_MOUSE_CALLBACK(WM_MOUSEMOVE, [](int x, int y, unsigned __int64 extra) {
        g_view_resources->m_arcball->drag(float(x), float(y));
        g_view_resources->m_cam->mouse_move(x, y);
        return 1; // indicate that we handled the mouse move event
        });
    REGISTER_MOUSE_CALLBACK(WM_LBUTTONDOWN, [](int x, int y, unsigned __int64 extra) {
        g_view_resources->m_arcball->beginDrag(float(x), float(y));
        return 1; // indicate that we handled the left button down event
        });
    REGISTER_MOUSE_CALLBACK(WM_LBUTTONUP, [](int x, int y, unsigned __int64 extra) {
        g_view_resources->m_arcball->endDrag();
        return 1; // indicate that we handled the left button up event
        });
    REGISTER_MOUSE_CALLBACK(WM_RBUTTONDOWN, [](int x, int y, unsigned __int64 extra) {
        g_view_resources->m_cam->begin_drag(x, y);
        return 1; // indicate that we handled the right button down event
        });
    REGISTER_MOUSE_CALLBACK(WM_RBUTTONUP, [](int x, int y, unsigned __int64 extra) {
        g_view_resources->m_cam->end_drag();
        return 1; // indicate that we handled the right button up event
        });
    REGISTER_MOUSE_CALLBACK(WM_MOUSEWHEEL, [](int delta, int ignore, unsigned __int64 extra) {
        g_view_resources->m_cam->zoom(float(delta));
        return 1; // indicate that we handled the mouse wheel event
        });
}

void theMeshView::create_model_view() {
    g_view_resources->m_draw_parts.clear();
    if (m_model) {
        const auto& parts = m_model->m_parts;
        for (btm::MeshExplicit<double>* part : parts) {
            std::unique_ptr<gl_prim> prim(create_prim(part));
            g_view_resources->m_draw_parts.push_back(std::move(prim));
        }
    }
}


// m_model_curvatures_k1
void theMeshView::createVertexK1View() {
    if (!m_model) {
        return; // No model loaded, cannot create face normals view
    }
    if (g_view_resources->m_model_curvatures_k1 == nullptr) {
        g_view_resources->m_model_curvatures_k1 = std::make_unique<visual_object>();
    }

    g_view_resources->m_model_curvatures_k1->clear_vao();
    g_view_resources->m_model_curvatures_k1->clear_mesh_data();
    const auto& parts = m_model->m_parts;
    for (btm::MeshExplicit<double>* part : parts) {
        float average_edge_length = static_cast<float>(part->average_edge_length());
        auto vertices = part->get_vertices();
        auto curvatures = part->get_vertex_curvatures();


        size_t num_vertices = vertices.size();
        for (std::uint32_t vertex_index = 0; vertex_index < num_vertices; ++vertex_index) {
            VertexExplicit<double> vertex = vertices[vertex_index];
            VertexCurvature<double> curvature = curvatures[vertex_index];

            fvec3 vertex_k1 = to_fvec3(curvature.k1_dir);
            fvec3 vertex_position = to_fvec3(vertex.position);
            fvec3 k1_end = vertex_position + vertex_k1 * (average_edge_length * 0.5f);
            g_view_resources->m_model_curvatures_k1->add_vector(vertex_position, k1_end);
        }
    }
    g_view_resources->m_model_curvatures_k1->create_prim();
    g_view_resources->m_model_curvatures_k1->set_color(fvec3(0.f, 0.f, 1.f)); // blue color for k1 directions
}
void theMeshView::createVertexK2View() {
    if (!m_model) {
        return; // No model loaded, cannot create face normals view
    }
    if (g_view_resources->m_model_curvatures_k2 == nullptr) {
        g_view_resources->m_model_curvatures_k2 = std::make_unique<visual_object>();
    }

    g_view_resources->m_model_curvatures_k2->clear_vao();
    g_view_resources->m_model_curvatures_k2->clear_mesh_data();
    const auto& parts = m_model->m_parts;
    for (btm::MeshExplicit<double>* part : parts) {
        float average_edge_length = static_cast<float>(part->average_edge_length());
        auto vertices = part->get_vertices();
        auto curvatures = part->get_vertex_curvatures();


        size_t num_vertices = vertices.size();
        for (std::uint32_t vertex_index = 0; vertex_index < num_vertices; ++vertex_index) {
            VertexExplicit<double> vertex = vertices[vertex_index];
            VertexCurvature<double> curvature = curvatures[vertex_index];

            fvec3 vertex_k2 = to_fvec3(curvature.k2_dir);
            fvec3 vertex_position = to_fvec3(vertex.position);
            fvec3 k2_end = vertex_position + vertex_k2 * (average_edge_length * 0.5f);
            g_view_resources->m_model_curvatures_k2->add_vector(vertex_position, k2_end);
        }
    }
    g_view_resources->m_model_curvatures_k2->create_prim();
    g_view_resources->m_model_curvatures_k2->set_color(fvec3(0.f, 1.f, 0.f)); // green color for k2 directions
}

void theMeshView::createFaceNormalsView(){
    if (!m_model) {
        return; // No model loaded, cannot create face normals view
    }
    if (g_view_resources->m_face_normals == nullptr) {
        g_view_resources->m_face_normals = std::make_unique<visual_object>();
    }

    g_view_resources->m_face_normals->clear_vao();
    g_view_resources->m_face_normals->clear_mesh_data();
    const auto& parts = m_model->m_parts;
    for (btm::MeshExplicit<double>* part : parts) {
        float average_edge_length = static_cast<float>(part->average_edge_length());
        for (std::uint32_t face_index = 0; face_index < part->faces.size(); ++face_index) {
            fvec3 face_normal = to_fvec3(part->faces[face_index].normal);
            fvec3 face_center = to_fvec3(part->faces[face_index].center);
            fvec3 normal_end = face_center + face_normal * (average_edge_length * 0.5f);
            g_view_resources->m_face_normals->add_vector(face_center, normal_end);
        }
    }
    g_view_resources->m_face_normals->create_prim();
    g_view_resources->m_face_normals->set_color(fvec3(1.f, 0.f, 0.f)); // red color for face normals
}

void theMeshView::createVertexNormalsView() {
    if (!m_model) {
        return; // No model loaded, cannot create vertex normals view
    }
    if (g_view_resources->m_vertex_normals == nullptr) {
        g_view_resources->m_vertex_normals = std::make_unique<visual_object>();
    }
    g_view_resources->m_vertex_normals->clear_vao();
    g_view_resources->m_vertex_normals->clear_mesh_data();
    const auto& parts = m_model->m_parts;
    for (btm::MeshExplicit<double>* part : parts) {
        float average_edge_length = static_cast<float>(part->average_edge_length());
        for (std::uint32_t vertex_index = 0; vertex_index < part->vertices.size(); ++vertex_index) {
            fvec3 vertex_position = to_fvec3(part->vertices[vertex_index].position);
            fvec3 vertex_normal = to_fvec3(part->vertices[vertex_index].normal);
            fvec3 normal_end = vertex_position + vertex_normal * (average_edge_length * 0.5f);
            g_view_resources->m_vertex_normals->add_vector(vertex_position, normal_end);
        }
    }
    g_view_resources->m_vertex_normals->create_prim();
    g_view_resources->m_vertex_normals->set_color(fvec3(0.f, 1.f, 0.f)); // green color for vertex normals
}

void theMeshView::reset_view() {
    g_view_resources->m_cam->set_position(btm::fvec3(0, 0, 25));
    g_view_resources->m_cam->set_target(btm::fvec3(0, 0, 0));
    g_view_resources->m_cam->set_up(btm::fvec3(0, 1, 0));
    g_view_resources->m_cam->set_fov(btm::dtr(20.f));

    double dx = m_model->m_bbox.max_p.x() - m_model->m_bbox.min_p.x();
    double dy = m_model->m_bbox.max_p.y() - m_model->m_bbox.min_p.y();
    float distance = 1.2f * g_view_resources->m_cam->computeCameraDistance(20, dx, dy);
    g_view_resources->m_cam->set_position(btm::fvec3(0, 0, distance));
    g_view_resources->m_cam->set_depth_range(0.1f, distance * 20.0f);

    g_view_resources->m_arcball->reset();
}

void theMeshView::init() {
    g_view_resources = new mesh_view_resources();
    g_view_resources->init();
    set_callbacks();
}

void theMeshView::render() {
    btm::GLContext* context = btm::get_current_gl_context();
    if (!context)
        return;

    context->begin_render();

    // OpenGL initialization
    // glEnable(GL_CULL_FACE);
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LEQUAL);
    glEnable(GL_MULTISAMPLE);

    glClearColor(0.2f, 0.4f, 0.6f, 1.f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    int width = context->width();
    int height = context->height();
    if (height <= 0)
        height = 1;

    g_view_resources->m_cam->set_aspect(width, height);
    g_view_resources->m_cam->set_viewport();

    fmat4 rot_mat(g_view_resources->m_arcball->rotation());
    // now it is safe to call resize on the arcball to update its internal viewport size,
    // which is used for mouse coordinate normalization
    g_view_resources->m_arcball->resize((float)width, (float)height);

    if (m_model) {
        g_view_resources->m_shader->use();

        g_view_resources->m_light->set_position(fvec3(-10, 0, 20));
        g_view_resources->m_light->set_color(fvec3(1.0f, 1.f, 1.f));

        g_view_resources->m_light->apply(g_view_resources->m_shader.get());
        g_view_resources->m_cam->apply(g_view_resources->m_shader.get());

        // render the mesh parts with the current rotation applied
        if (g_view_resources->m_draw_parts.size() > 0)
        {
            // render filled polygons first
            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
            glEnable(GL_DEPTH_TEST);

            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
            glEnable(GL_DEPTH_TEST);
            if (m_view_state.show_gaussian_curvature || m_view_state.show_mean_curvature || m_view_state.show_k1_k2_preview) {
                g_view_resources->m_shader->set_int("useVertexColor", 1);
            }
            for (const auto& part : g_view_resources->m_draw_parts) {
                part->set_draw_mode(GL_FILL);
                part->force_black = false;
                part->view_matrix = rot_mat;   // apply the current arcball rotation to the mesh parts
                part->render(g_view_resources->m_shader.get());
            }
            g_view_resources->m_shader->set_int("useVertexColor", 0);

            if (m_view_state.show_wireframe) {
                // then render wireframe on top
                glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
                glEnable(GL_POLYGON_OFFSET_LINE);
                glPolygonOffset(-1.0f, -1.0f); // pull lines toward camera
                for (const auto& part : g_view_resources->m_draw_parts) {
                    part->set_draw_mode(GL_LINE);
                    part->force_black = true; // render wireframe in black
                    part->set_use_vertex_color(0); // ensure vertex color is disabled for wireframe
                    part->render(g_view_resources->m_shader.get());
                }
            }
        }

        if (m_view_state.show_face_normals) {
            g_view_resources->m_face_normals->view_matrix = rot_mat;  // apply the same rotation to the face normals visualization
            g_view_resources->m_face_normals->render(g_view_resources->m_shader.get());
        }

        if (m_view_state.show_vertex_normals) {
            g_view_resources->m_vertex_normals->view_matrix = rot_mat;  // apply the same rotation to the vertex normals visualization
            g_view_resources->m_vertex_normals->render(g_view_resources->m_shader.get());
        }

        if (m_view_state.show_principal_curvatures) {
            g_view_resources->m_model_curvatures_k1->view_matrix = rot_mat;  // apply the same rotation to the principal curvatures visualization
            g_view_resources->m_model_curvatures_k1->render(g_view_resources->m_shader.get());
            g_view_resources->m_model_curvatures_k2->view_matrix = rot_mat;  // apply the same rotation to the principal curvatures visualization
            g_view_resources->m_model_curvatures_k2->render(g_view_resources->m_shader.get());
        }

        g_view_resources->m_shader->end();
        print_info();
    }
    // render coordinate system arrows
    if (g_view_resources->m_ucs_view.get() != nullptr) {
        g_view_resources->m_ucs_view->resize_window(width, height);
        g_view_resources->m_ucs_view->set_user_rotation(rot_mat); // apply the same rotation to the UCS view
        g_view_resources->m_ucs_view->render();
    }

    context->end_render();
}

void theMeshView::print_info() {
    g_view_resources->font2D->set_color(fvec4(1));
    g_view_resources->font2D->set_position(5, 5);
    g_view_resources->font2D->render((!m_model || m_model->get_name().empty()) ? "No model loaded" : m_model->get_name().c_str());
    if (!m_model) return;

    int ypos = 20;
    if (m_view_state.show_principal_curvatures) {
        g_view_resources->font2D->set_position(5, ypos);
        g_view_resources->font2D->render("Red: k1 direction, Green: k2 direction");
        ypos += 15;
    }
    if (m_view_state.show_curvture_map()) {
        g_view_resources->font2D->set_position(5, ypos);
        if (m_view_state.show_gaussian_curvature) {
            g_view_resources->font2D->render("Gaussian curvature map (red=negative(saddle), blue=near zero, green=positive(sphere-like))");
        }
        else if (m_view_state.show_mean_curvature) {
            g_view_resources->font2D->render("Mean curvature map (red=convex, blue=flat, green=concave)");
        }
        else if (m_view_state.show_k1_k2_preview) {
            g_view_resources->font2D->render("K1/K2 preview (red=concave ell, blue=concave cyl, green=saddle, yellow=flat, magenta=convex cyl, cyan=convex ell)");
        }
    }
}


void theMeshView::load_model(const std::string& filename) {
    if (m_model) {
        delete m_model;
        m_model = nullptr;
    }
    m_model = load_mesh_model(filename);
    create_model_view();
    g_view_resources->invalidate_visualizations();
    m_view_state.reset();
    reset_view();
}

void theMeshView::flip_mesh() {
    if (m_model)
    {
        for (btm::MeshExplicit<double>* part : m_model->m_parts) {
            part->flip_all_faces();
        }
        g_view_resources->m_draw_parts.clear();
        create_model_view();
    }
}

void theMeshView::calculate_curvatures() {
    if (m_model) {
        for (btm::MeshExplicit<double>* part : m_model->m_parts) {
            curvature::compute_vertex_curvatures<double>(*part);
        }
    }
}

void theMeshView::toggle_show_mesh() {
    m_view_state.show_wireframe = !m_view_state.show_wireframe;
}

void theMeshView::toggle_show_face_normals() {
    m_view_state.show_face_normals = !m_view_state.show_face_normals;
    if (m_view_state.show_face_normals) {
        createFaceNormalsView();
    }
}

void theMeshView::toggle_show_vertex_normals() {
    m_view_state.show_vertex_normals = !m_view_state.show_vertex_normals;
    if (m_view_state.show_vertex_normals) {
        createVertexNormalsView();
    }
}

void theMeshView::toggle_show_principal_curvatures(){
    m_view_state.show_principal_curvatures = !m_view_state.show_principal_curvatures;
    if (m_view_state.show_principal_curvatures) {
        createVertexK1View();
        createVertexK2View();
    }
}

void theMeshView::toggle_show_gaussian_curvatures(){
    toggle_map_view(MAP_TYPE_GAUSSIAN_CURVATURE);
}

void theMeshView::toggle_show_mean_curvatures(){
    toggle_map_view(MAP_TYPE_MEAN_CURVATURE);
}

void theMeshView::toggle_show_k1_k2_preview() {
    toggle_map_view(MAP_TYPE_K1_K2_PREVIEW);
}

bool theMeshView::toggle_map_view(int map_type) {
    if (!m_model) {
        return false; // no model loaded
    }

    bool map_view_flags[] = { m_view_state.show_gaussian_curvature, m_view_state.show_mean_curvature, m_view_state.show_k1_k2_preview };

    if (map_type < 1 || map_type > 3) {
        return false; // invalid map type
    }

    if (map_view_flags[map_type - 1]) {
        // check if map view is currently active, so disable it
        switch (map_type) {
        case MAP_TYPE_GAUSSIAN_CURVATURE:
            m_view_state.show_gaussian_curvature = false;
            break;
        case MAP_TYPE_MEAN_CURVATURE:
            m_view_state.show_mean_curvature = false;
            break;
        case MAP_TYPE_K1_K2_PREVIEW:
            m_view_state.show_k1_k2_preview = false;
            break;
        }
    }
    else {
        // map view is currently inactive, so enable it
        switch (map_type) {
        case MAP_TYPE_GAUSSIAN_CURVATURE:
            m_model->prepare_gaussian_curvature_preview();
            m_view_state.show_gaussian_curvature = true;
            m_view_state.show_mean_curvature = false;
            m_view_state.show_k1_k2_preview = false;
            break;
        case MAP_TYPE_MEAN_CURVATURE:
            m_model->prepare_mean_curvature_preview();
            m_view_state.show_gaussian_curvature = false;
            m_view_state.show_mean_curvature = true;
            m_view_state.show_k1_k2_preview = false;
            break;
        case MAP_TYPE_K1_K2_PREVIEW:
            m_model->prepare_k1_k2_preview();
            m_view_state.show_gaussian_curvature = false;
            m_view_state.show_mean_curvature = false;
            m_view_state.show_k1_k2_preview = true;
            break;
        }
        // rebuild model view to reflect the new curvature map visualization
        create_model_view(); 
    }
    return true;
}
