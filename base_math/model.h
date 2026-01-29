#ifndef __model_h__
#define __model_h__

#include <string>
#include <map>

#include "mesh.h"

namespace base_math {
    // template <typename T>
    // typedef std::map<std::string, half_edge_mesh<T>*> model;

    inline std::vector<std::string> splitString(const std::string& str, char delimiter = ' ') {
        std::vector<std::string> tokens;
        std::stringstream ss(str);
        std::string token;
        while (getline(ss, token, delimiter)) {
            tokens.push_back(token);
        }
        return tokens;
    }


#if 0
    template <typename T>
    half_edge_mesh<T>* parse_static_mesh(const static_mesh_def<T>& mesh_def) {
        half_edge_mesh<T>* mesh = new half_edge_mesh<T>;
        // add vertices
        size_t vertex_id = 1;
        for (const std::tuple<float, float, float>& vert : mesh_def.vertices) {
            mesh->add_vertex(vertex_id, std::get<0>(vert), std::get<1>(vert), std::get<2>(vert));
            ++vertex_id;
        }

        // add faces
        size_t face_id = 1;
        for (const std::tuple<size_t, size_t, size_t>& fc : mesh_def.faces) {
            mesh->add_face(face_id, std::get<0>(fc), std::get<1>(fc), std::get<2>(fc));
            ++face_id;
        }

        //mesh->orient_mesh();
        mesh->compute_face_normals();
        mesh->compute_face_centers();
        mesh->compute_vertex_normals();

        return mesh;
    }
#endif

    std::map<std::string, half_edge_mesh<float>*>* parse_obj_model(const std::string& fnm) {
        std::string line;
        std::ifstream mdl(fnm);

        size_t vertex_count = 0;
        size_t face_count = 0;

        std::map<std::string, half_edge_mesh<float>*>* model = new std::map<std::string, half_edge_mesh<float>*>;
        half_edge_mesh<float>* current_object = nullptr;
        if (mdl.is_open())
        {
            while (std::getline(mdl, line))
            {
                auto tokens = splitString(line);

                if (tokens[0] == "o") {
                    if (current_object) {
                        current_object->compute_face_normals();
                        current_object->compute_face_centers();
                        current_object->compute_vertex_normals();
                    }
                    // create new mesh
                    current_object = new half_edge_mesh<float>;
                    model[tokens[1]] = current_object;
                }
                else if (tokens[0] == "v") {
                    // add vertex
                    float x = (float)atof(tokens[1].c_str());
                    float y = (float)atof(tokens[2].c_str());
                    float z = (float)atof(tokens[3].c_str());
                    current_object->add_vertex((vertex_count, x, y, z));
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
                    current_object->add_face(face_count, idx1, idx2, idx3);
                    ++face_count;
                }
            }
            mdl.close();
            if (current_object) {
                current_object->compute_face_normals();
                current_object->compute_face_centers();
                current_object->compute_vertex_normals();
            }
        }
        else
        {
            return nullptr;
        }


        return model;
    }



}

#endif // __model_h__

