#define _CRT_SECURE_NO_WARNINGS
#include "gl_graphics.h"
#include "themeshview.h"
#include "themeshmodel.h"

#include "gl_camera.h"
#include "gl_shaders.h"
#include "gl_light.h"
#include "gl_prim.h"

#include "visual_objects.h"
#include "arcball.h"

#include "resource.h"

using namespace base_math;
using namespace base_opengl;

struct mesh_view_private {
    std::unique_ptr<gl_camera> m_cam;                      ///< gl_camera used to view the scene and compute view/projection.
    std::unique_ptr<gl_shader> m_shader;                ///< Shader program used for mesh and helper rendering.
    std::unique_ptr<gl_light> m_light;                  ///< Primary scene light affecting shading.

    std::unique_ptr<UCS_view> m_ucs_view;               ///< UCS (user coordinate system) view used to render the axes widget.

    // {{ the loaded model and the acompanying visualization objects
    std::vector<std::unique_ptr<gl_prim>> m_draw_parts; ///< Drawable mesh parts converted to OpenGL primitives.
    // }}

    std::vector<std::unique_ptr<visual_objects>> m_face_normals;        ///< Visualization of per-face normals as line segments.
    std::vector<std::unique_ptr<visual_objects>> m_model_curvatures_min;///< Visualization of minimum principal curvature directions per vertex.
    std::vector<std::unique_ptr<visual_objects>> m_model_curvatures_max;///< Visualization of maximum principal curvature directions per vertex.

    std::unique_ptr<arcball> m_arcball;                                 ///< Arcball for mouse interaction
};

meshViewWindow::meshViewWindow() : glViewWindow() {
    m_private = new mesh_view_private();
    m_private->m_ucs_view.reset(new UCS_view());
    m_private->m_arcball.reset(new arcball(800, 600));
    m_private->m_cam.reset(new gl_camera(fvec3(0, 0, 50), fvec3(0, 0, 0), fvec3(0, 1, 0)));
    m_private->m_cam->set_fov(dtr(30.f));
}

meshViewWindow::~meshViewWindow() {
    delete m_private;
}

void meshViewWindow::init() {
    m_private->m_ucs_view->initialize();

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
    init_done = true;
}

void meshViewWindow::reset_view() {
    m_private->m_cam->setup(fvec3(0, 0, 50), fvec3(0, 0, 0), fvec3(0, 1, 0));
    m_private->m_arcball->reset();
    m_private->m_model_curvatures_min.clear();
    m_private->m_model_curvatures_max.clear();
    m_view_state.reset();
    build_model_representation();
    build_normals_representation();
}

void meshViewWindow::build_model_representation() {
    m_private->m_draw_parts.clear();

    if (m_model) {
        const auto& parts = m_model->m_parts;
        m_private->m_draw_parts.reserve(parts.size());

        for (const base_math::mesh<double>* part : parts) {
            mesh_data md;
            collect_mesh_data(part, md);
            auto prim = std::make_unique<gl_prim>();
            prim->create_from_mesh(&md, GL_FILL);
            m_private->m_draw_parts.push_back(std::move(prim));
        }
    }
}

void meshViewWindow::build_normals_representation() {
    m_private->m_face_normals.clear();
    if (m_model) {
        const auto& parts = m_model->m_parts; // get_parts();
        m_private->m_face_normals.reserve(parts.size());

        for (base_math::mesh<double>* part : parts) {
            float vl = float(part->getAverageEdgeLength()) * 0.5f;
            part->computeFaceNormals();
            part->computeFaceProperties();
            auto mo = std::make_unique<visual_objects>();
            for (const auto& face : part->faces) {
                fvec3 face_normal = to_fvec3(face.second->normal);
                fvec3 face_center = to_fvec3(face.second->center);
                fvec3 normal_end = face_center + face_normal * vl;
                mo->add_vector(face_center, normal_end);
            }
            mo->create_prim();
            mo->get_prim()->set_color(fvec3(0.f, 0.f, 1.f));
            m_private->m_face_normals.push_back(std::move(mo));
        }
    }
}

LRESULT meshViewWindow::OnSize(int width, int height) {
    m_private->m_arcball->resize((float)width, (float)height);
    m_private->m_cam->set_aspect(width, height);
    if (m_private->m_ucs_view.get() != nullptr)
        m_private->m_ucs_view->resize_window(width, height);
    return 0;
}

void meshViewWindow::generate_model_curvature_view() {
    m_view_state.show_curvature = false; // ensure curvature visualization is disbled while we rebuild it

    m_private->m_model_curvatures_min.clear();
    m_private->m_model_curvatures_max.clear();
    if (m_model) {
        const auto& parts = m_model->m_parts; // get_parts();
        m_private->m_model_curvatures_min.reserve(parts.size());
        m_private->m_model_curvatures_max.reserve(parts.size());

        for (base_math::mesh<double>* part : parts) {
            double vl = part->getAverageEdgeLength() * 0.5f;
            auto mo_min = std::make_unique<visual_objects>();
            auto mo_max = std::make_unique<visual_objects>();
            // for each vertex, add curvature vectors
            for (const auto& v : part->vertices) {
                dvec3 p = v.second->position;
                dvec3 k1_end = p + v.second->curvature_info.k_min_dir * vl;
                dvec3 k2_end = p + v.second->curvature_info.k_max_dir * vl;
                // add min curvature vector in red
                mo_min->add_vector(to_fvec3(p), to_fvec3(k1_end));
                // add max curvature vector in green
                mo_max->add_vector(to_fvec3(p), to_fvec3(k2_end));
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

void meshViewWindow::generate_model_curvature_map_view() {
    m_model->prepare_curvature_preview();
    build_model_representation();
}

int meshViewWindow::onCommand(int cmd) {
    if (!init_done) return 0;
    if (m_model == nullptr) return 0;

    switch (cmd) {
    case ID_CALCULATE_MESHCURVATURES:
        m_model->compute_curvatures();
        generate_model_curvature_view();
        break;
    case ID_VIEW_MESH:
        m_view_state.show_wireframe = !m_view_state.show_wireframe;
        break;
    case ID_VIEW_FACENORMALS:
        m_view_state.show_normals = !m_view_state.show_normals;
        break;
    case ID_VIEW_PRINCIPALCURVATURES:
        if (m_model && m_model->curvatures_computed()) {
            m_view_state.show_curvature = !m_view_state.show_curvature;
        }
        break;
    case ID_VIEW_MAXIMUMCURVATUREMAP:
        if (m_model && m_model->curvatures_computed()) {
            generate_model_curvature_map_view();
            m_view_state.show_curvature_map = !m_view_state.show_curvature_map;
        }
        break;
    case ID_EDIT_FLIPMESH:
        m_model->flip_all_faces();
        build_model_representation(); // rebuild the OpenGL primitives to reflect the flipped geometry
        build_normals_representation();

        break;
    }

    return 0;
}

void meshViewWindow::onMouseWheel(int delta, unsigned __int64 extra_btn) {
    m_private->m_cam->zoom(float(delta));
}

void meshViewWindow::onMouseMove(int x, int y, unsigned __int64 extra) {
    m_private->m_arcball->drag(float(x), float(y));
    if (dragging) {
        int deltax = x - last_mouse_x;
        int deltay = y - last_mouse_y;

        m_private->m_cam->pan(float(deltax), float(deltay)); // adjust the camera position based on mouse movement

        last_mouse_x = x;
        last_mouse_y = y;
    }
}
void meshViewWindow::onLMouseDown(int x, int y, unsigned __int64 extra) {
    m_private->m_arcball->beginDrag(float(x), float(y));
}
void meshViewWindow::onLMouseUp(int x, int y, unsigned __int64 extra) {
    m_private->m_arcball->endDrag();
}
void meshViewWindow::onRMouseDown(int x, int y, unsigned __int64 extra) {
    dragging = true;
    last_mouse_x = x;
    last_mouse_y = y;
}
void meshViewWindow::onRMouseUp(int x, int y, unsigned __int64 extra) {
    dragging = false;
}

void meshViewWindow::render_scene(fmat4& cam_matrix, fmat4& rot_mat, gl_shader* shdr)
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
            //part->set_use_vertex_color(0); // ensure vertex color is disabled by default
            part->set_use_vertex_color(m_view_state.show_curvature_map ? 1 : 0); // ensure vertex color is disabled by default
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
                part->set_use_vertex_color(0); // ensure vertex color is disabled for wireframe
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

void meshViewWindow::render() {
    if (!hWnd || !hDC || !hGLRC) return;
    begin_render();

    m_private->m_cam->set_viewport();
    glClearColor(.25f, .25f, .25f, 1.f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    if (init_done) {
        fmat4 cam_matrix = m_private->m_cam->perspective();
        fmat4 rot_mat(m_private->m_arcball->rotation());

        render_scene(cam_matrix, rot_mat, m_private->m_shader.get());

        // render coordinate system arrows
        if (m_private->m_ucs_view.get() != nullptr) {
            m_private->m_ucs_view->set_user_rotation(rot_mat); // apply the same rotation to the UCS view
            m_private->m_ucs_view->render();
        }
    }

    end_render();
}