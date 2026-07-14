#include <fstream>

#include "themeshmodel.h"
#include "aabb.h"
#include "string_utils.h"
#include "mesh_loader.h"

using namespace btm;

static void recalculate_model(cModel* mdl) {
    double min_x = std::numeric_limits<double>::max();
    double min_y = std::numeric_limits<double>::max();
    double min_z = std::numeric_limits<double>::max();
    double max_x = std::numeric_limits<double>::lowest();
    double max_y = std::numeric_limits<double>::lowest();
    double max_z = std::numeric_limits<double>::lowest();
    AABB<double> bbox({ min_x, min_y, min_z }, { max_x, max_y, max_z });
    for (auto part : mdl->m_parts) {
        part->recalculateMesh();
        AABB<double> part_bbox = part->getBoundingBox();
        bbox = merge(bbox, part_bbox);
    }
    mdl->m_bbox = bbox;
    dvec3 bbcenter = center(bbox);
    for (auto part : mdl->m_parts) {
        part->generate_edges();
        part->build_adjacency();

        part->orient();
        double vol = part->component_signed_volume();
        if (vol < 0.0) {
            part->flip_all_faces();
        }

        part->build_adjacency();
        part->build_attributes();
        part->recalculateMesh();

        part->translate(dvec3(-bbcenter.x(), -bbcenter.y(), -bbcenter.z()));
    }
}

#ifdef _DEBUG
void log_model(cModel* mdl) {
    std::fstream log_file("c:\\temp\\model_log.txt", std::ios::out);
    log_file << "Model Name: " << mdl->get_name() << std::endl;
    log_file << "Bounding Box Min: (" << mdl->m_bbox.min_p.x() << ", " << mdl->m_bbox.min_p.y() << ", " << mdl->m_bbox.min_p.z() << ")" << std::endl;
    log_file << "Bounding Box Max: (" << mdl->m_bbox.max_p.x() << ", " << mdl->m_bbox.max_p.y() << ", " << mdl->m_bbox.max_p.z() << ")" << std::endl;
    for (size_t i = 0; i < mdl->m_parts.size(); ++i) {
        auto part = mdl->m_parts[i];
        log_file << "Part " << i + 1 << ":" << std::endl;
        log_file << "  Number of Vertices: " << part->num_vertices() << std::endl;
        log_file << "  Number of Faces: " << part->num_faces() << std::endl;
        log_file << "  Number of Edges: " << part->num_edges() << std::endl;
        log_file << "  Average Edge Length: " << part->m_average_edge_length << std::endl;
        for (auto& adj : part->vertex_adjacency) {
            log_file << "  Vertex Adjacency: " << adj.incident_faces.size() << " incident faces, " << adj.neighbor_vertices.size() << " neighbor vertices" << std::endl;
        }
    }
}
#endif

cModel* load_mesh_model(const std::string& fnm) {
    if (file_extension(fnm) == "obj") {
        // load obj cModel
        cModel* mdl = new cModel();
        if (load_obj<btm::MeshExplicit<double>, double>(fnm, mdl->m_parts)) {
            recalculate_model(mdl);
#ifdef _DEBUG
            // log_model(mdl);
#endif
            mdl->set_name(fnm);
            return mdl;
        }
        return mdl;
    }
    return nullptr;
}

void cModel::prepare_gaussian_curvature_preview() {
    if (m_parts.empty())
        return;
    // if (!curvatures_computed()) {
    //     return; // curvatures not computed, cannot prepare preview
    // }

    // For this preview, we will just set the curvature_map_value to:
    // 1 for positive Gaussian curvature (sphere like),
    // 0 for negative (saddle like),
    // and 0.5 for flat areas.
    std::fstream log_file("c:\\temp\\gaussian_curvature.txt", std::ios::out);
    
    for (auto& part : m_parts) {
        for (auto& v : part->vertex_curvatures) {
            log_file << "Vertex Gaussian Curvature: " << v.gaussian << std::endl;
            if (v.gaussian < -1.e-2) { // negative Gaussian curvature (saddle) in red
                v.curvature_map_value = 0;
            }
            else if (v.signGauss < 1.e-2) { // near zero Gaussian curvature (flat) in blue
                v.curvature_map_value = 1;
            }
            else
                v.curvature_map_value = 2;  // positive Gaussian curvature (sphere-like) in green
        }
    }
}

void cModel::prepare_k1_k2_preview() {
    if (m_parts.empty())
        return;
    // if (!curvatures_computed()) {
    //     return; // curvatures not computed, cannot prepare preview
    // }
    std::fstream log_file("c:\\temp\\k1k2_curvature.txt", std::ios::out);
    for (auto& part : m_parts) {
        for (auto& v : part->vertex_curvatures) {
            double k1 = v.k1;
            double k2 = v.k2;
            log_file << "Vertex k1: " << k1 << ", k2: " << k2 << std::endl;

            if (k1 < -1.e-2) { // k_1 is negative
                if (k2 < -1.e-2) { // k_2 is negative, Concave ellipsoid 0 red
                    v.curvature_map_value = 0;
                }
                else if (k2 < 1.e-2) { // k_2 is near zero, Concave cylinder 1 blue
                    v.curvature_map_value = 1;
                }
                else { // k_2 is positive, Hyperboloid 2 green
                    v.curvature_map_value = 2;
                }
            }
            else if (k1 < 1.e-2) { // k_1 is near zero
                if (k2 < -1.e-2) { // k_2 is negative, Concave cylinder 1 blue
                    v.curvature_map_value = 1;
                }
                else if (k2 < 1.e-2) { // k_2 is near zero, Plane 3 yellow
                    v.curvature_map_value = 3;
                }
                else { // k_2 is positive, Convex cylinder 4 magenta
                    v.curvature_map_value = 4;
                }
            }
            else { // k_1 is positive
                if (k2 < -1.e-2) { // k_2 is negative, Hyperboloid 2 green
                    v.curvature_map_value = 2;
                }
                else if (k2 < 1.e-2) { // k_2 is near zero, Convex cylinder 4 magenta
                    v.curvature_map_value = 4;
                }
                else { // k_2 is positive, Convex ellipsoid 5 cyan
                    v.curvature_map_value = 5;
                }
            }
        }
    }
}

void cModel::prepare_mean_curvature_preview() {
    if (m_parts.empty())
        return;
    // if (!curvatures_computed()) {
    //     return; // curvatures not computed, cannot prepare preview
    // }
    std::fstream log_file("c:\\temp\\mean_curvature.txt", std::ios::out);
    for (auto& part : m_parts) {
        for (auto& v : part->vertex_curvatures) {
            double H = v.mean;
            log_file << "Vertex Mean Curvature: " << H << std::endl;

            if (H < -1.e-2) { // negative concave areas in red
                v.curvature_map_value = 0;
            }
            else if (H < 1.e-2) { // near zero mean curvature flat areas in blue
                v.curvature_map_value = 1;
            }
            else { // positive convex areas in green
                v.curvature_map_value = 2;
            }
        }
    }
}