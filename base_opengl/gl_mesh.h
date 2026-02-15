#ifndef __gl_mesh_h__
#define __gl_mesh_h__

#include <vector>
#include "vector.h"

#include "he_mesh.h"

// namespace base_math {
//     /**
//     * forward declarations
//     */
//     template <typename T>
//     class half_edge_mesh;
// }


using namespace base_math;

namespace base_opengl {

    /**
     * @struct mesh_data
     * @brief Stores raw mesh data for rendering or processing.
     *
     * Contains vertex, normal, and index counts, as well as the corresponding data arrays.
     */
    struct mesh_data {
        size_t num_vertices=0;                ///< Number of vertices in the mesh.
        size_t num_normals=0;                 ///< Number of normals in the mesh.
        size_t num_indices=0;                 ///< Number of indices in the mesh.
        size_t num_curvatures=0;              ///< Number of curvature values (if available).
        std::vector<float> vertices;        ///< Flat array of vertex positions (x, y, z).
        std::vector<float> normals;         ///< Flat array of normal vectors (x, y, z).
        std::vector<unsigned int> indices;  ///< Indices defining mesh faces.
        std::vector<float> curvatures;      ///< Optional array of curvature values per vertex (if available).
    };

    /**
     * @class gl_mesh
     * @brief Represents a simple mesh structure for drawing.
     *
     * Stores vertices, normals, indices, and texture coordinates.
     */
    class gl_mesh {
    public:
        std::vector<fvec3> vertices;         ///< List of vertex positions.
        std::vector<fvec3> normals;          ///< List of normal vectors.
        std::vector<unsigned int> indices;   ///< Indices for mesh faces.

        /**
         * @brief Default constructor.
         */
        gl_mesh() {}

        /**
         * @brief Destructor.
         */
        ~gl_mesh() {}

        /**
         * @brief Returns the number of vertices in the mesh.
         * @return Number of vertices.
         */
        size_t n_vertices() {
            return vertices.size();
        }

        /**
         * @brief Adds a vertex to the mesh.
         * @param v The vertex position to add.
         */
        void addVertex(const fvec3& v) {
            vertices.push_back(v);
        }

        /**
         * @brief Adds a normal vector to the mesh.
         * @param n The normal vector to add.
         */
        void addNormal(const fvec3& n) {
            normals.push_back(n);
        }

        /**
         * @brief Adds a triangle face to the mesh using three indices.
         * @param i1 First index.
         * @param i2 Second index.
         * @param i3 Third index.
         */
        void addIndices(unsigned int i1, unsigned int i2, unsigned int i3) {
            indices.push_back(i1);
            indices.push_back(i2);
            indices.push_back(i3);
        }

        /**
         * @brief Adds an edge or line to the mesh using two indices.
         * @param i1 First index.
         * @param i2 Second index.
         */
        void addIndices(unsigned int i1, unsigned int i2) {
            indices.push_back(i1);
            indices.push_back(i2);
        }
    };

    /**
     * @brief Collects mesh data from a gl_mesh object into a mesh_data structure.
     * @param mesh Pointer to the gl_mesh object.
     * @param mdata Reference to the mesh_data structure to fill.
     * @return True if successful, false otherwise.
     */
    bool collect_mesh_data(gl_mesh* mesh, mesh_data& mdata);

    /**
     * @brief Collects mesh data from a half_edge_mesh object into a mesh_data structure.
     * @param mesh Pointer to the half_edge_mesh object.
     * @param mdata Reference to the mesh_data structure to fill.
     * @return True if successful, false otherwise.
     */
    template <typename T>
    bool collect_mesh_data(const half_edge_mesh<T>* mesh, mesh_data& mdata) {
        // build data from half-edge mesh based on its faces
        // this allows for proper duplication of vertices/normals per face
        // it is better for flat shading
        size_t index = mdata.vertices.size() / 3;
        const std::map<size_t, vertex<T>*>& vertices = mesh->get_vertices();
        for (auto f_pair : mesh->faces) {
            face<T>* f = f_pair.second;
            size_t v_ids[3] = { f->v1, f->v2, f->v3 };
            fvec3 face_normal = fvec3(float(f->normal.x()), float(f->normal.y()), float(f->normal.z()));
            for (int i = 0; i < 3; ++i) {
                dvec3 vert = vertices.at(v_ids[i])->coords;
                size_t i1 = index;
                mdata.vertices.push_back(float(vert.x()));
                mdata.vertices.push_back(float(vert.y()));
                mdata.vertices.push_back(float(vert.z()));
                mdata.normals.push_back(face_normal.x());
                mdata.normals.push_back(face_normal.y());
                mdata.normals.push_back(face_normal.z());
                mdata.indices.push_back(static_cast<unsigned int>(i1));
                if (mesh->curvatures_computed()) {
                    float curvature = float(vertices.at(v_ids[i])->curvature_data.absKmax);//  absGaussCurvature);
                    mdata.curvatures.push_back(curvature);
                    mdata.curvatures.push_back(curvature);
                    mdata.curvatures.push_back(0.75f);
                }
                ++index;
            }
        }

        mdata.num_vertices = mdata.vertices.size();
        mdata.num_normals = mdata.normals.size();
        mdata.num_indices = mdata.indices.size();
        mdata.num_curvatures = mdata.curvatures.size();

        return true;
    }


    /**
     * @brief Creates a simple mesh representing the UCS (Universal Coordinate System).
     * @return Pointer to the created gl_mesh object.
     */
    gl_mesh* create_UCS_mesh();

}

#endif // __gl_mesh_h__
