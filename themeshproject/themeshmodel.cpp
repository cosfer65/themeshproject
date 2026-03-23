#include <Windows.h>
#include "themeshmodel.h"

#include <fstream>
#include "base_definitions.h"
#include "nurbs_core.h"

using namespace base_math;

void compute_vertex_curvatures(mesh<double>* fmesh);
void cModel::compute_curvatures() {
    if (state.curvatures_computed) {
        return; // curvatures already calculated, no need to recompute
    }
    for (auto part : m_parts) {
        compute_vertex_curvatures(part);
    }
    state.curvatures_computed = true; // mark curvatures as calculated
}

void cModel::prepare_curvature_preview() {
    if (m_parts.empty())
        return;
    if (!curvatures_computed()) {
        return; // curvatures not computed, cannot prepare preview
    }

    // For this preview, we will just set the curvature_map_value to:
    // 1 for positive Gaussian curvature (sphere like),
    // 0 for negative (saddle like),
    // and 0.5 for flat areas.
    int negatives = 0;
    int positives = 0;
    for (auto& part : m_parts) {
        for (const auto& v : part->vertices) {
            v.second->curvature_map_value = 0.5;
            if (v.second->curvature_info.signGauss == 1) {
                v.second->curvature_map_value = 1;
                positives++;
            }
            else if (v.second->curvature_info.signGauss == -1) {
                v.second->curvature_map_value = 0;
                negatives++;
            }
        }
    }
    OutputDebugString(("Curvature preview prepared: " + std::to_string(positives) + " positive curvature vertices, " + 
        std::to_string(negatives) + " negative curvature vertices.\n").c_str());
}


//////////////////////////////////////////////////////////////////////////////////////////////////////


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

static dvec3 bbox_center(const basematrix<double, 2, 3>& bbox) {
    return dvec3(
        (bbox(0, 0) + bbox(1, 0)) * 0.5f,
        (bbox(0, 1) + bbox(1, 1)) * 0.5f,
        (bbox(0, 2) + bbox(1, 2)) * 0.5f
    );
}

static void recalculate_model(cModel* mdl) {
    basematrix<double, 2, 3> bbox(0);
    for (auto part : mdl->m_parts) {
        part->recalculateMesh();
        basematrix<double, 2, 3> part_bbox = part->getBoundingBox();
        expand_bbox(bbox, part_bbox);
    }
    dvec3 bbcenter = bbox_center(bbox);
    for (auto part : mdl->m_parts) {
        part->breakQuads();
        part->translate(dvec3(-bbcenter.x(), -bbcenter.y(), -bbcenter.z()));
        // part->computeFaceNormals();
        // part->computeFaceProperties();
        // part->computeVertexNormals();
    }
}

static bool load_obj(const std::string& fnm, cModel* mdl) {
    std::string line;
    std::ifstream mdl_file(fnm);
    size_t vertex_count = 1;
    size_t face_count = 1;
    base_math::mesh<double>* mesh = nullptr;

    if (mdl_file.is_open()) {
        while (std::getline(mdl_file, line)) {
            auto tokens = splitString(line);
            if (tokens.empty()) {
                continue; // skip empty lines
            }
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

static void file_dump_model(const std::string& fnm, const cModel* mdl) {
    std::ofstream out(fnm);
    if (!out.is_open()) {
        return;
    }
    for (const auto& part : mdl->m_parts) {
        out << "o Part" << std::endl;
        for (const auto& v_pair : part->vertices) {
            const basevector<double, 3>& pos = v_pair.second->position;
            out << "v " << pos.x() << " " << pos.y() << " " << pos.z() << std::endl;
        }
        for (const auto& f_pair : part->faces) {
            const meshFace<double>* face = f_pair.second;
            out << "f";
            for (const auto& vert : face->vertices) {
                out << " " << vert->id; // OBJ format is 1-based
            }
            out << std::endl;
        }
    }
    out.close();
}

static bool load_nurbs(const std::string& fnm, cModel* mdl) {
    std::vector<base_math::mesh<double>*> meshes = create_from_nurbs_file<double>(fnm, 15, 15);
    if (!meshes.empty()) {
        for (auto mesh : meshes) {
            mdl->add_part(mesh);
        }
        recalculate_model(mdl);
#ifdef _DEBUG
        file_dump_model("c:\\temp\\tessellated_output.obj", mdl);
#endif
        return true;
    }

    return false;
}

cModel* load_mesh_model(const std::string& fnm) {
    if (file_extension(fnm) == "obj") {
        // load obj cModel
        cModel* mdl = new cModel();
        if (load_obj(fnm, mdl)) {
            return mdl;
        }
        return mdl;
    }
    else if (file_extension(fnm) == "nurbs") {
        cModel* mdl = new cModel();
        if (load_nurbs(fnm, mdl)) {
            return mdl;
        }
        return mdl;
    }
    return nullptr;
}

base_math::mesh<double>* create_sphere_mesh(double radius) {
    base_math::mesh<double>* ms = new base_math::mesh<double>;

    int sectorCount = 48;
    int stackCount = 48;

    double x, y, z, xz;                              // vertex position
    double lengthInv = 1.0f / radius;    // normal
    double sectorStep = 2 * PI<double> / sectorCount;
    double stackStep = PI<double> / stackCount;
    double sectorAngle, stackAngle;
    int vertexCount = 0;
    for (int i = 0; i <= stackCount; ++i)
    {
        stackAngle = -(PI<double> / 2 - i * stackStep);        // starting from -pi/2 to pi/2
        xz = radius * cos(stackAngle);
        y = radius * sin(stackAngle);

        // add (sectorCount+1) vertices per stack
        // the first and last vertices have same position and normal, but different tex coords
        for (int j = 0; j <= sectorCount; ++j)
        {
            sectorAngle = j * sectorStep;           // starting from 0 to 2pi

            // vertex position
            x = xz * cos(sectorAngle);             // r * cos(u) * cos(v)
            z = xz * sin(sectorAngle);             // r * cos(u) * sin(v)
            ms->addVertex(vertexCount, dvec3(x, y, z));

            vertexCount++;
        }
    }

    // indices
    //  k1--k1+1
    //  |  / |
    //  | /  |
    //  k2--k2+1
    int k1, k2;
    for (int i = 0; i < stackCount; ++i)
    {
        k1 = i * (sectorCount + 1);     // beginning of current stack
        k2 = k1 + sectorCount + 1;      // beginning of next stack

        for (int j = 0; j < sectorCount; ++j, ++k1, ++k2)
        {
            // 2 triangles per sector excluding 1st and last stacks
            if (i != 0)
                ms->addFace(k1, k2, k1 + 1);   // k1---k2---k1+1

            if (i != (stackCount - 1))
                ms->addFace(k1 + 1, k2, k2 + 1); // k1+1---k2---k2+1
        }
    }

    return ms;
}


cModel* generate_sphere_model(){
    cModel* mdl = new cModel();
    base_math::mesh<double>* sphere_mesh = create_sphere_mesh(1.0);
    mdl->add_part(sphere_mesh);
    recalculate_model(mdl);
    return mdl;
}