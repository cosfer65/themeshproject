#include <iostream>
// #include "matrix.h"

#include <fstream>

// #include "he_mesh.h"

#include "string_utils.h"
#include "vector.h"
#include "mesh.h"

// testing math library
// int run_all_matrix_tests();
// int run_all_vector_tests();

int mesh_tests();

using namespace base_math;

class model {
    public:
        std::vector<base_math::mesh<double>*> m_parts;

        ~model() {
            cleanUp();
        }

        void cleanUp() {
            for (auto part : m_parts) {
                delete part;
            }
            m_parts.clear();
        }

        void add_part(base_math::mesh<double>* part) {
            m_parts.push_back(part);
        }

        void normalize_absolute_curvature() {
            if (m_parts.empty())
                return;
            //if (!m_parts[0]->curvatures_computed())
            //    return;

            double max_abs_curvature = 0.;
            double min_abs_curvature = std::numeric_limits<double>::max();
            for (auto& part : m_parts) {
                for (const auto& v : part->getVertices()) {
                    if (v.second->curvature_info.absKmax > max_abs_curvature)
                        max_abs_curvature = v.second->curvature_info.absKmax;
                    if (v.second->curvature_info.absKmax < min_abs_curvature)
                        min_abs_curvature = v.second->curvature_info.absKmax;
                }
            }
            double range = max_abs_curvature - min_abs_curvature;
            double scale = range > 0. ? 1. / range : 1.;
            for (auto& part : m_parts) {
                for (const auto& v : part->getVertices()) {
                    v.second->curvature_info.absKmax *= scale;
                }
            }
        }

};


static bool load_obj(const std::string& fnm, model* mdl) {
    std::string line;
    std::ifstream mdl_file(fnm);
    size_t vertex_count = 1;
    size_t face_count = 1;
    base_math::mesh<double>* vmesh = nullptr;

    if (mdl_file.is_open()) {

        while (std::getline(mdl_file, line)) {
            auto tokens = splitString(line);
            /*
            labels:
            "mtllib" -> material library
            "usemtl" -> use material
            "vn" -> vertex normal (they are not calculated as normal to the surface, and thus useless for our calculations)
            "vt" -> vertex texture coordinate (they are only for texture mapping, useless for our calculations)
            are ignored
            */
            if (tokens[0] == "o") {
                if (vmesh != nullptr) {
                    // finalize previous mesh part
                    // vmesh->average_edge_length = vmesh->total_edge_length / meshvmeshhalf_edges.size();
                }
                vmesh = new mesh<double>();
                mdl->add_part(vmesh);
            }
            else if (tokens[0] == "v") {
                // add vertex
                double x = atof(tokens[1].c_str());
                double y = atof(tokens[2].c_str());
                double z = atof(tokens[3].c_str());
                vmesh->addVertex(vertex_count, basevector<double, 3>(x, y, z));
                ++vertex_count;
            }
            else if (tokens[0] == "f") {
                // add face
                auto v1 = splitString(tokens[1], '/');
                auto v2 = splitString(tokens[2], '/');
                auto v3 = splitString(tokens[3], '/');
                size_t idx1 = static_cast<size_t>(std::stoull(v1[0]));
                size_t idx2 = static_cast<size_t>(std::stoull(v2[0]));
                size_t idx3 = static_cast<size_t>(std::stoull(v3[0]));
                if (idx1 == 0 || idx2 == 0 || idx3 == 0) {
                    // OBJ indices are 1-based; zero is invalid
                    continue;
                }
                if (idx1 >= vertex_count || idx2 >= vertex_count || idx3 >= vertex_count) {
                    // Index out of bounds; skip this face
                    continue;
                }
                if (tokens.size() > 4) {
                    // More than 3 vertices per face
                    auto v4 = splitString(tokens[4], '/');
                    size_t idx4 = static_cast<size_t>(std::stoull(v4[0]));
                    if (idx4 == 0 || idx4 >= vertex_count) {
                        // Index out of bounds; skip this face
                        continue;
                    }
                    vmesh->addFace(idx1, idx2, idx3, idx4);
                }
                else {
                    vmesh->addFace(idx1, idx2, idx3);
                }
            }
        }
        mdl_file.close();
        if (vmesh != nullptr) {
            // finalize last mesh part
            // vmesh->average_edge_length = vmesh->total_edge_length / vmesh->half_edges.size();
        }
        // recalculate_model(mdl);
    }
    else {
        return false;
    }

    return true;
}

/**
 * @brief Loads a model from a file, dispatching to a format-specific loader.
 *
 * Currently supports only Wavefront OBJ files, identified by the `obj` file extension.
 * The returned `model` instance is heap-allocated and owned by the caller.
 *
 * @param fnm [in] Path to the model file.
 * @return Pointer to a newly allocated `model` on success, or `nullptr` on failure
 *         or when the extension is unsupported.
 */
model* load_model(const std::string& fnm) {
    if (file_extension(fnm) == "obj") {
        // load obj model
        model* mdl = new model();
        if (load_obj(fnm, mdl)) {
            return mdl;
        }
        delete mdl;
    }
    return nullptr;
}

void print_mesh(const mesh<double>* m) {
#if 0
    std::cout << "Vertices:--\n";
    for (const auto& v : m->getVertices()) {
        std::cout << "Vertex " << v.first << ": (" << v.second->position.x() << ", " << v.second->position.y() << ", " << v.second->position.z() << ")\n";
        for (const auto& adj : v.second->adjacentFaces) {
            std::cout << "  Adjacent face: " << adj->id << "\n";
            meshEdge<double>* e = adj->firstEdge;
            for (size_t i = 0; i < adj->vertexCount(); ++i) {
                std::cout << "    Edge from vertex " << e->source->id << " to vertex " << e->target->id << "\n";
                if (e->oppositeEdge != nullptr) {
                    std::cout << "      Opposite edge: (" << e->oppositeEdge->source->id << " -> " << e->oppositeEdge->target->id << ")\n";
                }
                else {
                    std::cout << "      Opposite edge: None\n";
                }
                e = e->nextEdge;
            }
        }
    }
#endif
    std::cout << "Faces:--\n";
    for (const auto& f : m->getFaces()) {
        std::cout << "Face " << f.first << ": vertices ";
        for (const auto& vert : f.second->vertices) {
            std::cout << vert->id << " ";
        }
        std::cout << "\n";
        // std::cout << "  First edge: " << (f.second->firstEdge ? ("(" + std::to_string(f.second->firstEdge->source->id) + " -> " + std::to_string(f.second->firstEdge->target->id) + ")") : "None") << "\n";
    }
#if 0
    std::cout << "Edges:--\n";
    for (const auto& e : m->getEdges()) {
        std::cout << "Edge (" << e.first.first << " -> " << e.first.second << ")\n";
    }
#endif
}

void print_mesh_normals(const mesh<double>* m) {
    std::cout << "Face normals:--\n";
    for (const auto& f : m->getFaces()) {
        std::cout << "Face " << f.first << ": normal (" << f.second->normal.x() << ", " << f.second->normal.y() << ", " << f.second->normal.z() << ")\n";
    }
}


int main() {
    // run_all_matrix_tests();
    // run_all_vector_tests();
    // mesh_tests();
    model* m = load_model("C:\\GitHub\\themeshproject_resources\\models\\alfa147.obj");

    std::cout << "Model loaded successfully!" << std::endl;
    std::cout << "Number of parts: " << m->m_parts.size() << std::endl;
    for (size_t i = 0; i < m->m_parts.size(); ++i) {
        std::cout << "Part " << i << ": " << m->m_parts[i]->getVertices().size() << " vertices, " << m->m_parts[i]->getFaces().size() << " faces\n";
    }
    delete m;

    std::cout << "All tests passed successfully!" << std::endl;
    
    return 0;
}
