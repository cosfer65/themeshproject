/// \file model.h
/// \brief Declares the `model` class and the `load_model` function used to
///        represent and load multi-part 3D meshes.
///
/// This header defines a lightweight container for multiple mesh parts based on
/// `base_math::half_edge_mesh<double>` and a factory function that loads such
/// models from an external file.

#ifndef __model_loader_h__
#define __model_loader_h__

#include <vector>
#include <string>

#include "mesh.h"

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
        if (!m_parts[0]->curvatures_computed())
            return;

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

/// \brief Loads a `model` instance from the specified file.
///
/// This function parses the file identified by \p fnm and constructs a
/// `model` composed of one or more `base_math::half_edge_mesh<double>` parts,
/// depending on the file contents.
///
/// \param fnm Path to the model file to be loaded.
/// \return A pointer to a newly created `model` instance on success; behavior
///         on failure depends on the implementation (may return `nullptr` or
///         throw an exception).
model* load_model(const std::string& fnm);

#endif // __model_loader_h__
