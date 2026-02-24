// Loads and normalizes 3D mesh models from supported file formats (currently Wavefront OBJ).
#include <fstream>

//#include "he_mesh.h"

#include "model.h"
#include "string_utils.h"
#include "geometry.h"

using namespace base_math;

/**
 * @brief Expands an aggregate axis-aligned bounding box to include a part's bounding box.
 *
 * Given the current model bounding box `bbox` and a part bounding box `part_bbox`,
 * this function updates `bbox` so that it becomes the minimum AABB enclosing both.
 *
 * @param bbox      [in,out] Aggregate bounding box of the model, as a 2x3 matrix.
 *                  Row 0 contains the minimum (x, y, z), row 1 the maximum (x, y, z).
 * @param part_bbox [in]     Bounding box of a single mesh part in the same format.
 */
static void expand_bbox(basematrix<double, 2, 3>& bbox, const basematrix<double, 2, 3>& part_bbox) {
    for (size_t i = 0; i < 3; ++i) {
        if (part_bbox(0, i) < bbox(0, i)) {
            bbox(0, i) = part_bbox(0, i);
        }
        if (part_bbox(1, i) > bbox(1, i)) {
            bbox(1, i) = part_bbox(1, i);
        }
    }
}

/**
 * @brief Computes the center point of an axis-aligned bounding box.
 *
 * The center is computed as the midpoint between the min and max
 * coordinates along each axis.
 *
 * @param bbox [in] Bounding box as a 2x3 matrix (row 0 = min, row 1 = max).
 * @return Center of the bounding box as an `dvec3`.
 */
static dvec3 bbox_center(const basematrix<double, 2, 3>& bbox) {
    return dvec3(
        (bbox(0, 0) + bbox(1, 0)) * 0.5f,
        (bbox(0, 1) + bbox(1, 1)) * 0.5f,
        (bbox(0, 2) + bbox(1, 2)) * 0.5f
    );
}

/**
 * @brief Recenters and recomputes derived properties for all parts of a model.
 *
 * This function:
 *  - Computes the combined bounding box of all mesh parts in the model.
 *  - Recenters the model so that the bounding box center lies at the origin.
 *  - Recomputes face normals, face properties, and vertex normals for each part.
 *
 * @param mdl [in,out] Pointer to the model to be updated. Must not be null.
 */
static void recalculate_model(model* mdl) {
    basematrix<double, 2, 3> bbox(0);
    for (auto part : mdl->m_parts) {
        basematrix<double, 2, 3> part_bbox = part->getBoundingBox();
        expand_bbox(bbox, part_bbox);
    }
    dvec3 bbcenter = bbox_center(bbox);
    for (auto part : mdl->m_parts) {
        part->breakQuads();
        part->translate(dvec3(-bbcenter.x(), -bbcenter.y(), -bbcenter.z()));
        part->computeFaceNormals();
        part->computeFaceProperties();
        part->computeVertexNormals();
    }
}

/**
 * @brief Loads a Wavefront OBJ file into a `model` instance.
 *
 * This loader:
 *  - Parses lines of the OBJ file using `splitString`.
 *  - Creates a new `half_edge_mesh` for each `o` object declaration and adds it as a part.
 *  - Adds vertices for `v` lines (positions only; normals/texcoords are ignored).
 *  - Adds triangular faces for `f` lines, parsing only the vertex indices and ignoring
 *    texture/normal indices.
 *  - Computes `average_edge_length` for each mesh part.
 *  - Calls `recalculate_model` to recenter and update derived properties.
 *
 * Notes:
 *  - Material (`mtllib`, `usemtl`), vertex normals (`vn`), and texture coordinates (`vt`)
 *    are intentionally ignored as they are not needed for geometric calculations here.
 *  - Vertex and face indices are assumed to be 1-based, consistent with OBJ format.
 *
 * @param fnm [in]  Path to the OBJ file.
 * @param mdl [out] Pre-allocated model that will receive the mesh parts.
 * @return `true` on successful load and processing, `false` if the file could not be opened.
 */
static bool load_obj(const std::string& fnm, model* mdl) {
    std::string line;
    std::ifstream mdl_file(fnm);
    size_t vertex_count = 1;
    size_t face_count = 1;
    base_math::mesh<double>* mesh = nullptr;

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
                if (mesh != nullptr) {
                    // finalize previous mesh part
                    // mesh->average_edge_length = mesh->total_edge_length / mesh->getEdges().size();
                }
                mesh = new base_math::mesh<double>();
                mdl->add_part(mesh);
            }
            else if (tokens[0] == "v") {
                // add vertex
                double x = atof(tokens[1].c_str());
                double y = atof(tokens[2].c_str());
                double z = atof(tokens[3].c_str());
                mesh->addVertex(vertex_count, base_math::dvec3(x, y, z));
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
                    mesh->addFace(idx1, idx2, idx3, idx4);
                }
                else {
                    mesh->addFace(idx1, idx2, idx3);
                }

                ++face_count;
            }
        }
        mdl_file.close();
        // finalize last mesh part
        recalculate_model(mdl);
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