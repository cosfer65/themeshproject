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

/// \class model
/// \brief Aggregates multiple half-edge mesh parts into a single logical model.
///
/// The `model` class holds a collection of pointers to `base_math::half_edge_mesh<double>`
/// instances, each representing a separate part of a 3D model (for example,
/// different sub-meshes or components). The class provides methods to add parts
/// and to retrieve the complete list of stored parts.
///
/// \note The `model` class does not manage the lifetime of the mesh pointers
///       it stores. Callers are responsible for ensuring that the pointed-to
///       meshes remain valid for the lifetime of the `model` instance.
class model {
    /// \brief Collection of mesh parts that make up the model.
    ///
    /// Each entry points to a `base_math::half_edge_mesh<double>` instance.
    std::vector<base_math::half_edge_mesh<double>*> m_parts;
public:
    /// \brief Constructs an empty model with no mesh parts.
    model() {}

    /// \brief Destroys the model instance.
    ///
    /// This destructor deletes the stored mesh pointers.
    ~model() {
        for (auto part : m_parts) {
            delete part;
        }
    }

    /// \brief Adds a mesh part to the model.
    ///
    /// The provided pointer is stored as-is; ownership is not transferred.
    ///
    /// \param part Pointer to a `base_math::half_edge_mesh<double>` instance
    ///             representing a model part. Must remain valid for as long
    ///             as the model uses it.
    void add_part(base_math::half_edge_mesh<double>* part) {
        m_parts.push_back(part);
    }

    /// \brief Provides read-only access to the collection of mesh parts.
    ///
    /// \return A constant reference to the internal vector of mesh pointers.
    const std::vector<base_math::half_edge_mesh<double>*>& get_parts() const {
        return m_parts;
    }

    /// \brief Normalizes the absolute Gaussian curvature values across all mesh parts.
    ///
    /// This method rescales the `absGaussCurvature` attribute stored in each vertex's
    /// `curvature_data` so that the values are linearly mapped into the range
    /// \f$[0,1]\f$ over the entire model.
    ///
    /// The procedure is:
    ///  - If there are no parts in the model, the function returns immediately.
    ///  - If curvature data has not been computed for the first part
    ///    (as indicated by `curvatures_computed()`), the function also returns
    ///    without modifying any data.
    ///  - It then iterates over all vertices in all parts to determine the global
    ///    minimum and maximum of `absGaussCurvature`.
    ///  - Using these extrema, it computes a scale factor \f$\frac{1}{\text{max} - \text{min}}\f$.
    ///    If the range is zero (all vertices have identical curvature), a scale factor
    ///    of `1.0f` is used to avoid division by zero.
    ///  - Finally, it multiplies each vertex's `absGaussCurvature` by the computed
    ///    scale factor, effectively normalizing the values relative to the global
    ///    range of the model.
    ///
    /// \note This operation is destructive: the original absolute Gaussian curvature
    ///       values are overwritten with their normalized counterparts.
    ///       Callers that require the original values should back them up beforehand.
    void normalize_absolute_curvature() {
        if (m_parts.empty())
            return;
        if (!m_parts[0]->curvatures_computed())
            return;
        double max_abs_gauss = 0.;
        double min_abs_gauss = std::numeric_limits<double>::max();
        for (auto& part : m_parts) {
            for (const auto& v : part->get_vertices()) {
                if (v.second->curvature_data.absGaussCurvature > max_abs_gauss)
                    max_abs_gauss = v.second->curvature_data.absGaussCurvature;
                if (v.second->curvature_data.absGaussCurvature < min_abs_gauss)
                    min_abs_gauss = v.second->curvature_data.absGaussCurvature;
            }
        }
        double range = max_abs_gauss - min_abs_gauss;
        double scale = range > 0. ? 1. / range : 1.;
        for (auto& part : m_parts) {
            for (const auto& v : part->get_vertices()) {
                v.second->curvature_data.absGaussCurvature *= scale;
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
