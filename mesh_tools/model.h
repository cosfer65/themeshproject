/// \file model.h
/// \brief Declares the `model` class and the `load_model` function used to
///        represent and load multi-part 3D meshes.
///
/// This header defines a lightweight container for multiple mesh parts based on
/// `base_math::half_edge_mesh<float>` and a factory function that loads such
/// models from an external file.

#ifndef __model_loader_h__
#define __model_loader_h__

#include <vector>
#include <string>

#include "mesh.h"

/// \class model
/// \brief Aggregates multiple half-edge mesh parts into a single logical model.
///
/// The `model` class holds a collection of pointers to `base_math::half_edge_mesh<float>`
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
    /// Each entry points to a `base_math::half_edge_mesh<float>` instance.
    std::vector<base_math::half_edge_mesh<float>*> m_parts;
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
    /// \param part Pointer to a `base_math::half_edge_mesh<float>` instance
    ///             representing a model part. Must remain valid for as long
    ///             as the model uses it.
    void add_part(base_math::half_edge_mesh<float>* part) {
        m_parts.push_back(part);
    }

    /// \brief Provides read-only access to the collection of mesh parts.
    ///
    /// \return A constant reference to the internal vector of mesh pointers.
    const std::vector<base_math::half_edge_mesh<float>*>& get_parts() const {
        return m_parts;
    }
};

/// \brief Loads a `model` instance from the specified file.
///
/// This function parses the file identified by \p fnm and constructs a
/// `model` composed of one or more `base_math::half_edge_mesh<float>` parts,
/// depending on the file contents.
///
/// \param fnm Path to the model file to be loaded.
/// \return A pointer to a newly created `model` instance on success; behavior
///         on failure depends on the implementation (may return `nullptr` or
///         throw an exception).
model* load_model(const std::string& fnm);

#endif // __model_loader_h__
