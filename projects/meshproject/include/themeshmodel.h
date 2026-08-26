#pragma once

#include <string>
#include <vector>

#include "aabb.h"
#include "mesh_explicit.h"

class cModel {
    std::string m_name; ///< Name of the model, used for identification and display purposes.
public:
    btm::AABB<double> m_bbox;
    std::vector<btm::MeshExplicit<double>*> m_parts;

    ~cModel() {
        clean_up();
    }
    void clean_up() {
        for (auto part : m_parts) {
            delete part;
        }
        m_parts.clear();
    }
    void add_part(btm::MeshExplicit<double>* part) {
        m_parts.push_back(part);
    }
    const std::string& get_name() const {
        return m_name;
    }
    void set_name(const std::string& name) {
        m_name = name;
    }
    bool curvatures_calculated() const {
        if (m_parts.empty()) return false;
        return m_parts[0]->curvatures_calculated();
    }

    void prepare_gaussian_curvature_preview();
    void prepare_k1_k2_preview();
    void prepare_mean_curvature_preview();
};

cModel* load_mesh_model(const std::string& fnm);
