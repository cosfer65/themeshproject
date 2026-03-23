#pragma once

#include "mesh.h"

struct model_state {
    bool curvatures_computed = false; ///< Flag indicating whether curvature has been calculated for the current model.
    void reset() {
        curvatures_computed = false; // reset curvature calculation state when the model changes
    }
};


class cModel {
    model_state state; ///< State information about the model, such as whether curvature has been calculated.
    public:
        std::vector<base_math::mesh<double>*> m_parts;
        ~cModel() {
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
        void compute_curvatures();
        void flip_all_faces() {
            for (auto part : m_parts) {
                part->flip_all_faces();
            }
        }
        bool curvatures_computed() const {
            return state.curvatures_computed;
        }

        void prepare_curvature_preview();
};

cModel* load_mesh_model(const std::string& fnm);
cModel* generate_sphere_model();