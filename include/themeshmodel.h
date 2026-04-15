#pragma once

#include "mesh.h"

class cModel {
    bool m_curvatures_computed = false; ///< Flag indicating whether curvature has been calculated for the current model.
	std::string m_name; ///< Name of the model, used for identification and display purposes.
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
            return m_curvatures_computed;
        }
        bool& curvatures_computed(){
            return m_curvatures_computed;
        }

        void prepare_gaussian_curvature_preview();
        void prepare_k1_k2_preview();
        void prepare_mean_curvature_preview();
        const std::string& get_name() const {
            return m_name;
		}
        void set_name(const std::string& name) {
            m_name = name;
		}
};

cModel* load_mesh_model(const std::string& fnm);
cModel* generate_sphere_model();