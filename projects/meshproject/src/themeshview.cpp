#include <iostream>
#include "themeshview.h"
#include "gl_context.h"
#include "glew.h"
#include "application.h"
#include "cmd_target.h"

#include "view_resources.h"

#include "themeshmodel.h"
#include "curvature.h"
#include "feature_line_extraction.h"

// cannot be a member of theMeshView because it is used in static callbacks, so we use a global instance
static mesh_view_resources* g_view_resources;

theMeshView::theMeshView(cModel* model/*=nullptr*/) : m_model(model) {
    m_view_resources = new mesh_view_resources();
    m_view_resources->init();
    g_view_resources = m_view_resources;
    set_callbacks();
}

void theMeshView::set_callbacks() {
    REGISTER_MOUSE_CALLBACK(WM_MOUSEMOVE, [](int x, int y, unsigned __int64 extra) {
        g_view_resources->m_arcball.drag(float(x), float(y));
        g_view_resources->m_cam.mouse_move(x, y);
        return 1; // indicate that we handled the mouse move event
        });
    REGISTER_MOUSE_CALLBACK(WM_LBUTTONDOWN, [](int x, int y, unsigned __int64 extra) {
        g_view_resources->m_arcball.beginDrag(float(x), float(y));
        return 1; // indicate that we handled the left button down event
        });
    REGISTER_MOUSE_CALLBACK(WM_LBUTTONUP, [](int x, int y, unsigned __int64 extra) {
        g_view_resources->m_arcball.endDrag();
        return 1; // indicate that we handled the left button up event
        });
    REGISTER_MOUSE_CALLBACK(WM_RBUTTONDOWN, [](int x, int y, unsigned __int64 extra) {
        g_view_resources->m_cam.begin_drag(x, y);
        return 1; // indicate that we handled the right button down event
        });
    REGISTER_MOUSE_CALLBACK(WM_RBUTTONUP, [](int x, int y, unsigned __int64 extra) {
        g_view_resources->m_cam.end_drag();
        return 1; // indicate that we handled the right button up event
        });
    REGISTER_MOUSE_CALLBACK(WM_MOUSEWHEEL, [](int delta, int ignore, unsigned __int64 extra) {
        g_view_resources->m_cam.zoom(float(delta));
        return 1; // indicate that we handled the mouse wheel event
        });
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

    m_view_resources->m_cam.set_aspect(width, height);
    m_view_resources->m_cam.set_viewport();

    fmat4 rot_mat(m_view_resources->m_arcball.rotation());
    // now it is safe to call resize on the arcball to update its internal viewport size,
    // which is used for mouse coordinate normalization
    m_view_resources->m_arcball.resize((float)width, (float)height);

    if (m_model) {
        m_view_resources->m_shader.use();

        m_view_resources->m_light.set_position(fvec3(-10, 0, 20));
        m_view_resources->m_light.set_color(fvec3(1.0f, 1.f, 1.f));

        m_view_resources->m_light.apply(&m_view_resources->m_shader);
        m_view_resources->m_cam.apply(&m_view_resources->m_shader);

        glEnable(GL_DEPTH_TEST);

        // render the mesh parts with the current rotation applied
        if (m_view_resources->m_draw_parts.size() > 0)
        {
            if (m_view_state.show_mesh) {
                // render filled polygons first
                glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

                if (m_view_state.show_curvture_map()) {
                    m_view_resources->m_shader.set_int("useVertexColor", 1);
                }
                for (const auto& part : m_view_resources->m_draw_parts) {
                    part->set_draw_mode(GL_FILL);
                    part->force_black = false;
                    part->view_matrix = rot_mat;   // apply the current arcball rotation to the mesh parts
                    part->render(&m_view_resources->m_shader);
                }
                m_view_resources->m_shader.set_int("useVertexColor", 0);
            }

            if (m_view_state.show_wireframe) {
                // then render wireframe on top
                glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
                glEnable(GL_POLYGON_OFFSET_LINE);
                glPolygonOffset(-1.0f, -1.0f); // pull lines toward camera
                for (const auto& part : m_view_resources->m_draw_parts) {
                    part->set_draw_mode(GL_LINE);
                    part->force_black = true; // render wireframe in black
                    part->view_matrix = rot_mat;   // apply the current arcball rotation to the mesh parts
                    part->set_use_vertex_color(0); // ensure vertex color is disabled for wireframe
                    part->render(&m_view_resources->m_shader);
                }
                glDisable(GL_POLYGON_OFFSET_LINE);
            }
        }

        //glDisable(GL_DEPTH_TEST);

        if (m_view_state.show_face_normals) {
            m_view_resources->m_face_normals.view_matrix = rot_mat;  // apply the same rotation to the face normals visualization
            m_view_resources->m_face_normals.render(&m_view_resources->m_shader);
        }

        if (m_view_state.show_vertex_normals) {
            m_view_resources->m_vertex_normals.view_matrix = rot_mat;  // apply the same rotation to the vertex normals visualization
            m_view_resources->m_vertex_normals.render(&m_view_resources->m_shader);
        }

        if (m_view_state.show_principal_k1) {
            m_view_resources->m_model_curvatures_k1.view_matrix = rot_mat;  // apply the same rotation to the principal curvatures visualization
            m_view_resources->m_model_curvatures_k1.render(&m_view_resources->m_shader);
        }
        if (m_view_state.show_principal_k2) {
            m_view_resources->m_model_curvatures_k2.view_matrix = rot_mat;  // apply the same rotation to the principal curvatures visualization
            m_view_resources->m_model_curvatures_k2.render(&m_view_resources->m_shader);
        }

        if (m_view_state.show_ridges || m_view_state.show_creases || m_view_state.show_valleys)
        {
            glEnable(GL_POLYGON_OFFSET_LINE);
            glDisable(GL_LIGHTING);
            glLineWidth(3.0f);
            glPolygonOffset(-3.0f, -3.0f); // pull lines toward camera
            if (m_view_state.show_ridges) {
                // apply the model rotation to the feature lines visualization
                m_view_resources->m_feature_lines[0].view_matrix = rot_mat;
                // render the feature lines with the current shader,
                // the graphics engine uses the GL_LINES primitive type to render the feature lines
                m_view_resources->m_feature_lines[0].render(&m_view_resources->m_shader);
            }
            if (m_view_state.show_valleys) {
                // apply the model rotation to the feature lines visualization
                m_view_resources->m_feature_lines[1].view_matrix = rot_mat;
                // render the feature lines with the current shader,
                // the graphics engine uses the GL_LINES primitive type to render the feature lines
                m_view_resources->m_feature_lines[1].render(&m_view_resources->m_shader);
            }
            if (m_view_state.show_creases) {
                // apply the model rotation to the feature lines visualization
                m_view_resources->m_feature_lines[2].view_matrix = rot_mat;
                // render the feature lines with the current shader,
                // the graphics engine uses the GL_LINES primitive type to render the feature lines
                m_view_resources->m_feature_lines[2].render(&m_view_resources->m_shader);
            }
            glDisable(GL_POLYGON_OFFSET_LINE);
            glLineWidth(1.0f);
        }
        glLineWidth(3.0f);
        m_view_resources->m_potential_ridges.view_matrix = rot_mat;
        m_view_resources->m_potential_ridges.render(&m_view_resources->m_shader);
        
        m_view_resources->m_potential_valleys.view_matrix = rot_mat;
        m_view_resources->m_potential_valleys.render(&m_view_resources->m_shader);
        
        m_view_resources->m_potential_creases.view_matrix = rot_mat;
        m_view_resources->m_potential_creases.render(&m_view_resources->m_shader);

        m_view_resources->m_boundary_edges.view_matrix = rot_mat;
        m_view_resources->m_boundary_edges.render(&m_view_resources->m_shader);
        glLineWidth(1.0f);

        glEnable(GL_DEPTH_TEST);

        m_view_resources->m_shader.end();
        print_info();
    }
    // render coordinate system arrows
    m_view_resources->m_ucs_view.resize_window(width, height);
    m_view_resources->m_ucs_view.set_user_rotation(rot_mat); // apply the same rotation to the UCS view
    m_view_resources->m_ucs_view.render();

    context->end_render();
}

void theMeshView::load_model(const std::string& filename) {
    if (m_model) {
        delete m_model;
        m_model = nullptr;
    }
    m_model = load_mesh_model(filename);
    create_model_view();
    m_view_resources->invalidate_visualizations();
    m_view_state.reset();
    reset_view();
}

void theMeshView::flip_mesh() {
    if (m_model)
    {
        for (btm::MeshExplicit<double>* part : m_model->m_parts) {
            part->flip_all_faces();
        }
        m_view_resources->m_draw_parts.clear();
        create_model_view();
    }
}

void theMeshView::reset_view_state() {
    m_view_state.reset();
}

void theMeshView::calculate_curvatures() {
    if (m_model) {
        for (btm::MeshExplicit<double>* part : m_model->m_parts) {
            curvature::compute_vertex_curvatures<double>(*part);
        }
    }
}

void theMeshView::calculate_feature_lines() {
    if (m_model) {
        for (int i = 0; i < 3; ++i) {
            m_view_resources->m_feature_lines[i].clear();
        }
        for (btm::MeshExplicit<double>* part : m_model->m_parts) {
            curvature::compute_vertex_curvatures<double>(*part);
            FeatureLineParameters<double> flp;
            FeatureLineExtractor<double> extractor(*part, part->vertex_curvatures, flp);

            std::vector<FeatureLine<double>> ridges = extractor.extractRidges();
            std::vector<FeatureLine<double>> valleys = extractor.extractValleys();
            std::vector<FeatureLine<double>> creases = extractor.extractCreases();

            // add feature lines to the draw object
            for (const auto& ridge : ridges) {
                for (size_t i = 0; i < ridge.points.size() - 1; ++i) {
                    fvec3 v1 = dvec_to_fvec<3>(ridge.points[i]);
                    fvec3 v2 = dvec_to_fvec<3>(ridge.points[i + 1]);
                    m_view_resources->m_feature_lines[0].add_vector(v1, v2);
                }
            }
            for (const auto& valley : valleys) {
                for (size_t i = 0; i < valley.points.size() - 1; ++i) {
                    fvec3 v1 = dvec_to_fvec<3>(valley.points[i]);
                    fvec3 v2 = dvec_to_fvec<3>(valley.points[i + 1]);
                    m_view_resources->m_feature_lines[1].add_vector(v1, v2);
                }
            }
            for (const auto& crease : creases) {
                for (size_t i = 0; i < crease.points.size() - 1; ++i) {
                    fvec3 v1 = dvec_to_fvec<3>(crease.points[i]);
                    fvec3 v2 = dvec_to_fvec<3>(crease.points[i + 1]);
                    m_view_resources->m_feature_lines[2].add_vector(v1, v2);
                }
            }
        }
        for (int i = 0; i < 3; ++i) {
            m_view_resources->m_feature_lines[i].create_prim();
        }
        m_view_resources->m_feature_lines[0].set_color(fvec3(1.f, 0.f, 0.f)); // red for ridges
        m_view_resources->m_feature_lines[1].set_color(fvec3(0.f, 0.f, 1.f)); // blue for valleys
        m_view_resources->m_feature_lines[2].set_color(fvec3(0.f, 1.f, 0.f)); // green for creases
        for (int i = 0; i < 3; ++i)
            m_view_resources->m_feature_lines[i].clear_mesh_data();
    }
}

#ifdef _DEBUG

#include "feature_line_params.h"
static FeatureLineParameters<double> g_feature_line_params;

// Check if vertex satisfies ridge/valley criteria
bool _isRidgeVertex(CurvatureField<double>& field, int vertex) {
    double k = field[vertex].kmax;
    return std::abs(k) > g_feature_line_params.ridge_threshold && k > 0.0;
}
bool _isValleyVertex(CurvatureField<double>& field, int vertex) {
    double k = field[vertex].kmax;
    return std::abs(k) > g_feature_line_params.valley_threshold && k < 0.0;
}

double _computeDihedralAngle(btm::MeshExplicit<double>* part, int edgeIndex) {
    const auto& e = part->edge(edgeIndex);
    const auto& ef = e.incident_faces;

    int f0 = ef[0]; // face 0
    int f1 = ef[1]; // face 1

    // Boundary edge: no dihedral (or treat as 0)
    if (f0 < 0 || f1 < 0) {
        //odprintf("Warning: Edge %d is a boundary edge; dihedral angle undefined. Returning 0.\n", edgeIndex);
        return 0.0;
    }

    const basevec3<double>& n0 = part->faces[f0].normal;
    const basevec3<double>& n1 = part->faces[f1].normal;

    double dihedral = btm::dihedralAngle<double>(n0, n1);
    return rtd<double>(dihedral);
}

void theMeshView::test_function() {
    if (!m_model) {
        return; // no model loaded
    }
    m_view_resources->m_potential_ridges.clear();
    m_view_resources->m_potential_valleys.clear();
    m_view_resources->m_potential_creases.clear();
    m_view_resources->m_boundary_edges.clear();
    m_view_resources->m_potential_ridges.clear_do();
    m_view_resources->m_potential_valleys.clear_do();
    m_view_resources->m_potential_creases.clear_do();
    m_view_resources->m_boundary_edges.clear_do();

    for (btm::MeshExplicit<double>* part : m_model->m_parts) {
        curvature::compute_vertex_curvatures<double>(*part);
        for (std::uint32_t v = 0; v < part->num_vertices(); ++v) {
            if (_isRidgeVertex(part->vertex_curvatures, v)) {
                const auto& vertex = part->vertices[v];
                fvec3 pos = dvec_to_fvec<3>(vertex.position);
                m_view_resources->m_potential_ridges.add_vertex(pos);
            }
            if (_isValleyVertex(part->vertex_curvatures, v)) {
                const auto& vertex = part->vertices[v];
                fvec3 pos = dvec_to_fvec<3>(vertex.position);
                m_view_resources->m_potential_valleys.add_vertex(pos);
            }
        }
    }

    int num_edges = 0;
    for (btm::MeshExplicit<double>* part : m_model->m_parts) {
        for (std::uint32_t v = 0; v < part->num_edges(); ++v) {
            const auto& e = part->edge(v);
            ++num_edges;
            if (e.is_boundary()) {
                const auto& vertex0 = part->vertices[e.v0()];
                const auto& vertex1 = part->vertices[e.v1()];
                fvec3 pos0 = dvec_to_fvec<3>(vertex0.position);
                fvec3 pos1 = dvec_to_fvec<3>(vertex1.position);
                m_view_resources->m_boundary_edges.add_vector(pos0, pos1);
            }

            double da = _computeDihedralAngle(part, v);
            if (da > g_feature_line_params.dihedralAngleThreshold) {
                const auto& vertex0 = part->vertices[e.v0()];
                const auto& vertex1 = part->vertices[e.v1()];
                fvec3 pos0 = dvec_to_fvec<3>(vertex0.position);
                fvec3 pos1 = dvec_to_fvec<3>(vertex1.position);
                m_view_resources->m_potential_creases.add_vector(pos0, pos1);
            }
        }
    }
    m_view_resources->m_potential_ridges.create_prim(GL_POINTS);
    m_view_resources->m_potential_ridges.set_color(fvec3(1.f, 0.f, 0.f));
    m_view_resources->m_potential_ridges.set_draw_mode(GL_POINT);
    m_view_resources->m_potential_ridges.set_draw_type(GL_POINTS);
    m_view_resources->m_potential_ridges.clear_mesh_data();

    m_view_resources->m_potential_valleys.create_prim(GL_POINTS);
    m_view_resources->m_potential_valleys.set_color(fvec3(0.f, 0.f, 1.f));
    m_view_resources->m_potential_valleys.set_draw_mode(GL_POINT);
    m_view_resources->m_potential_valleys.set_draw_type(GL_POINTS);
    m_view_resources->m_potential_valleys.clear_mesh_data();

    m_view_resources->m_potential_creases.create_prim(GL_LINES);
    m_view_resources->m_potential_creases.set_color(fvec3(0.f, 1.f, 0.f));
    m_view_resources->m_potential_creases.clear_mesh_data();

    m_view_resources->m_boundary_edges.create_prim(GL_LINES);
    m_view_resources->m_boundary_edges.set_color(fvec3(1.f, 1.f, 0.f));
    m_view_resources->m_boundary_edges.clear_mesh_data();
}
#endif