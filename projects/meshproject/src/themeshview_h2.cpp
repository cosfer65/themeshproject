#include "themeshview.h"
#include "themeshmodel.h"
#include "draw_objects.h"
#include "view_resources.h"

#include "prim.h"
// #include "draw_objects.h"




void theMeshView::toggle_show_wireframe() {
    m_view_state.show_wireframe = !m_view_state.show_wireframe;
}
void theMeshView::toggle_show_mesh() {
    m_view_state.show_mesh = !m_view_state.show_mesh;
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

void theMeshView::toggle_show_principal_k1() {
    if (!m_model->curvatures_calculated()) {
        m_view_state.show_principal_k1 = false;
        return;
    }

    m_view_state.show_principal_k1 = !m_view_state.show_principal_k1;
    if (m_view_state.show_principal_k1) {
        createVertexK1View();
    }
}

void theMeshView::toggle_show_principal_k2() {
    if (!m_model->curvatures_calculated()) {
        m_view_state.show_principal_k2 = false;
        return;
    }
    m_view_state.show_principal_k2 = !m_view_state.show_principal_k2;
    if (m_view_state.show_principal_k2) {
        createVertexK2View();
    }
}

void theMeshView::toggle_show_gaussian_curvatures() {
    toggle_map_view(MAP_TYPE_GAUSSIAN_CURVATURE);
}

void theMeshView::toggle_show_mean_curvatures() {
    toggle_map_view(MAP_TYPE_MEAN_CURVATURE);
}

void theMeshView::toggle_show_k1_k2_preview() {
    toggle_map_view(MAP_TYPE_K1_K2_PREVIEW);
}

bool theMeshView::toggle_map_view(int map_type) {
    if (!m_model) {
        return false; // no model loaded
    }
    if (!m_model->curvatures_calculated()) {
        m_view_state.show_gaussian_curvature = false;
        m_view_state.show_mean_curvature = false;
        m_view_state.show_k1_k2_preview = false;
        return false;
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

void theMeshView::toggle_show_ridges() {
    if (!m_model->curvatures_calculated()) {
        m_view_state.show_ridges = false;
        return;
    }
    m_view_state.show_ridges = !m_view_state.show_ridges;
}

void theMeshView::toggle_show_valleys() {
    if (!m_model->curvatures_calculated()) {
        m_view_state.show_valleys = false;
        return;
    }
    m_view_state.show_valleys = !m_view_state.show_valleys;
}

void theMeshView::toggle_show_creases() {
    if (!m_model->curvatures_calculated()) {
        m_view_state.show_creases = false;
        return;
    }
    m_view_state.show_creases = !m_view_state.show_creases;
}

void theMeshView::toggle_show_boundaries() {
    m_view_state.show_boundaries = !m_view_state.show_boundaries;
    if (m_view_state.show_boundaries) {
        createBoundariesView();
    }
}

void theMeshView::create_model_view() {
    m_view_resources->m_draw_parts.clear();
    if (m_model) {
        const auto& parts = m_model->m_parts;
        for (btm::MeshExplicit<double>* part : parts) {
            std::unique_ptr<gl_prim> prim(create_prim(part));
            m_view_resources->m_draw_parts.push_back(std::move(prim));
        }
    }
}

// m_model_curvatures_k1
void theMeshView::createVertexK1View() {
    if (!m_model) {
        return; // No model loaded, cannot create face normals view
    }

    m_view_resources->m_model_curvatures_k1.clear();

    if (!m_model->curvatures_calculated()) {
        return;
    }

    const auto& parts = m_model->m_parts;
    for (btm::MeshExplicit<double>* part : parts) {
        float average_edge_length = static_cast<float>(part->average_edge_length());
        auto vertices = part->get_vertices();
        auto curvatures = part->get_vertex_curvatures();

        size_t num_vertices = vertices.size();
        for (std::uint32_t vertex_index = 0; vertex_index < num_vertices; ++vertex_index) {
            VertexExplicit<double> vertex = vertices[vertex_index];
            VertexCurvature<double> curvature = curvatures[vertex_index];

            fvec3 vertex_k1 = dvec_to_fvec<3>(curvature.kmin_dir);
            fvec3 vertex_position = dvec_to_fvec<3>(vertex.position);
            fvec3 k1_end = vertex_position + vertex_k1 * (average_edge_length * 0.5f);
            m_view_resources->m_model_curvatures_k1.add_vector(vertex_position, k1_end);
        }
    }
    m_view_resources->m_model_curvatures_k1.create_prim();
    m_view_resources->m_model_curvatures_k1.set_color(fvec3(0.f, 0.f, 1.f)); // blue color for k1 directions
    m_view_resources->m_model_curvatures_k1.clear_mesh_data();
}
void theMeshView::createVertexK2View() {
    if (!m_model) {
        return; // No model loaded, cannot create face normals view
    }

    m_view_resources->m_model_curvatures_k2.clear();

    if (!m_model->curvatures_calculated())
        return;

    const auto& parts = m_model->m_parts;
    for (btm::MeshExplicit<double>* part : parts) {
        float average_edge_length = static_cast<float>(part->average_edge_length());
        auto vertices = part->get_vertices();
        auto curvatures = part->get_vertex_curvatures();

        size_t num_vertices = vertices.size();
        for (std::uint32_t vertex_index = 0; vertex_index < num_vertices; ++vertex_index) {
            VertexExplicit<double> vertex = vertices[vertex_index];
            VertexCurvature<double> curvature = curvatures[vertex_index];

            fvec3 vertex_k2 = dvec_to_fvec<3>(curvature.kmax_dir);
            fvec3 vertex_position = dvec_to_fvec<3>(vertex.position);
            fvec3 k2_end = vertex_position + vertex_k2 * (average_edge_length * 0.5f);
            m_view_resources->m_model_curvatures_k2.add_vector(vertex_position, k2_end);
        }
    }
    m_view_resources->m_model_curvatures_k2.create_prim();
    m_view_resources->m_model_curvatures_k2.set_color(fvec3(0.f, 1.f, 0.f)); // green color for k2 directions
    m_view_resources->m_model_curvatures_k2.clear_mesh_data();
}

void theMeshView::createFaceNormalsView() {
    if (!m_model) {
        return; // No model loaded, cannot create face normals view
    }

    m_view_resources->m_face_normals.clear();
    const auto& parts = m_model->m_parts;
    for (btm::MeshExplicit<double>* part : parts) {
        float average_edge_length = static_cast<float>(part->average_edge_length());
        for (std::uint32_t face_index = 0; face_index < part->faces.size(); ++face_index) {
            fvec3 face_normal = dvec_to_fvec<3>(part->faces[face_index].normal);
            fvec3 face_center = dvec_to_fvec<3>(part->faces[face_index].center);
            fvec3 normal_end = face_center + face_normal * (average_edge_length * 0.5f);
            m_view_resources->m_face_normals.add_vector(face_center, normal_end);
        }
    }
    m_view_resources->m_face_normals.create_prim();
    m_view_resources->m_face_normals.set_color(fvec3(1.f, 0.f, 0.f)); // red color for face normals
    m_view_resources->m_face_normals.clear_mesh_data();
}

void theMeshView::createVertexNormalsView() {
    if (!m_model) {
        return; // No model loaded, cannot create vertex normals view
    }
    m_view_resources->m_vertex_normals.clear();
    const auto& parts = m_model->m_parts;
    for (btm::MeshExplicit<double>* part : parts) {
        float average_edge_length = static_cast<float>(part->average_edge_length());
        for (std::uint32_t vertex_index = 0; vertex_index < part->vertices.size(); ++vertex_index) {
            fvec3 vertex_position = dvec_to_fvec<3>(part->vertices[vertex_index].position);
            fvec3 vertex_normal = dvec_to_fvec<3>(part->vertices[vertex_index].normal);
            fvec3 normal_end = vertex_position + vertex_normal * (average_edge_length * 0.5f);
            m_view_resources->m_vertex_normals.add_vector(vertex_position, normal_end);
        }
    }
    m_view_resources->m_vertex_normals.create_prim();
    m_view_resources->m_vertex_normals.set_color(fvec3(0.f, 1.f, 0.f)); // green color for vertex normals
    m_view_resources->m_vertex_normals.clear_mesh_data();
}

void theMeshView::createBoundariesView() {
    if (!m_model) {
        return; // No model loaded, cannot create boundaries view
    }

    m_view_resources->m_boundary_edges.clear();
    m_view_resources->m_boundary_edges.clear_do();
    for (btm::MeshExplicit<double>* part : m_model->m_parts) {
        for (std::uint32_t v = 0; v < part->num_edges(); ++v) {
            const auto& e = part->edge(v);
            if (e.is_boundary()) {
                const auto& vertex0 = part->vertices[e.v0()];
                const auto& vertex1 = part->vertices[e.v1()];
                fvec3 pos0 = dvec_to_fvec<3>(vertex0.position);
                fvec3 pos1 = dvec_to_fvec<3>(vertex1.position);
                m_view_resources->m_boundary_edges.add_vector(pos0, pos1);
            }
        }
    }
    m_view_resources->m_boundary_edges.create_prim();
    m_view_resources->m_boundary_edges.set_color(fvec3(1.f, 1.f, 0.f)); // yellow color for boundary edges
    m_view_resources->m_boundary_edges.clear_mesh_data();
}   

void theMeshView::reset_view() {
    m_view_resources->m_cam.set_position(btm::fvec3(0, 0, 25));
    m_view_resources->m_cam.set_target(btm::fvec3(0, 0, 0));
    m_view_resources->m_cam.set_up(btm::fvec3(0, 1, 0));
    m_view_resources->m_cam.set_fov(btm::dtr(20.f));

    if (m_model) {
        double dx = m_model->m_bbox.max_p.x() - m_model->m_bbox.min_p.x();
        double dy = m_model->m_bbox.max_p.y() - m_model->m_bbox.min_p.y();
        float distance = 1.2f * m_view_resources->m_cam.computeCameraDistance(20, dx, dy);
        m_view_resources->m_cam.set_position(btm::fvec3(0, 0, distance));
        m_view_resources->m_cam.set_depth_range(0.1f, distance * 20.0f);
    }

    m_view_resources->m_arcball.reset();
}

void theMeshView::print_info() {
    m_view_resources->font2D->set_color(fvec4(1));
    m_view_resources->font2D->set_position(5, 5);
    m_view_resources->font2D->render((!m_model || m_model->get_name().empty()) ? "No model loaded" : m_model->get_name().c_str());
    if (!m_model) return;

    int ypos = 20;
    if (m_view_state.show_principal_k1 || m_view_state.show_principal_k2) {
        m_view_resources->font2D->set_position(5, ypos);
        m_view_resources->font2D->render("Blue: k1 direction, Green: k2 direction");
        ypos += 15;
    }
    if (m_view_state.show_curvture_map()) {
        m_view_resources->font2D->set_position(5, ypos);
        if (m_view_state.show_gaussian_curvature) {
            m_view_resources->font2D->render("Gaussian curvature map (red=negative(saddle), blue=near zero, green=positive(sphere-like))");
        }
        else if (m_view_state.show_mean_curvature) {
            m_view_resources->font2D->render("Mean curvature map (red=convex, blue=flat, green=concave)");
        }
        else if (m_view_state.show_k1_k2_preview) {
            m_view_resources->font2D->render("K1/K2 preview (red=concave ell, blue=concave cyl, green=saddle, yellow=flat, magenta=convex cyl, cyan=convex ell)");
        }
    }
}
