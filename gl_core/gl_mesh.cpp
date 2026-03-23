#include "mesh.h"
#include "gl_mesh.h"

namespace base_opengl {

    /**
     * @brief Collects vertex, normal, and index data from a gl_mesh and stores it in a mesh_data structure.
     * 
     * This function extracts all vertex positions, normals, and indices from the given gl_mesh and
     * stores them in the provided mesh_data object.
     * 
     * @param mesh Pointer to the gl_mesh containing the mesh data.
     * @param mdata Reference to the mesh_data structure to populate.
     * @return true if data collection is successful.
     */
    bool collect_mesh_data(gl_mesh* mesh, mesh_data& mdata) {
        for (auto v : mesh->vertices) {
            mdata.vertices.push_back(v.x());
            mdata.vertices.push_back(v.y());
            mdata.vertices.push_back(v.z());
        }
        for (auto n : mesh->normals) {
            mdata.normals.push_back(n.x());
            mdata.normals.push_back(n.y());
            mdata.normals.push_back(n.z());
        }
        for (auto idx : mesh->indices) {
            mdata.indices.push_back(static_cast<unsigned int>(idx));
        }

        mdata.num_vertices = mdata.vertices.size();
        mdata.num_normals = mdata.normals.size();
        mdata.num_indices = mdata.indices.size();
        mdata.num_curvatures = 0; // gl_mesh does not have curvature data

        return true;
    }

    /**
     * @brief Creates a gl_mesh representing the UCS (Universal Coordinate System) axes.
     * 
     * The mesh consists of three colored axes: X (red), Y (green), and Z (blue), each represented by lines.
     * Additional vertices are added to create arrowheads for each axis.
     * 
     * @return Pointer to the newly created gl_mesh representing the UCS.
     */
    gl_mesh* create_UCS_mesh() {
        gl_mesh* ms = new gl_mesh;
        ms->addVertex(fvec3(0, 0, 0));
        // x->red
        ms->addVertex(fvec3(1, 0, 0));
        ms->addVertex(fvec3(0.8f, 0.2f, 0));
        ms->addVertex(fvec3(0.8f, -0.2f, 0));
        ms->addIndices(0, 1);
        ms->addIndices(1, 2);
        ms->addIndices(1, 3);

        // y->green
        ms->addVertex(fvec3(0, 1, 0));
        ms->addVertex(fvec3(0.2f, 0.8f, 0));
        ms->addVertex(fvec3(-0.2f, 0.8f, 0));
        ms->addIndices(0, 4);
        ms->addIndices(4, 5);
        ms->addIndices(4, 6);

        // z->blue
        ms->addVertex(fvec3(0, 0, 1));
        ms->addVertex(fvec3(0, 0.2f, 0.8f));
        ms->addVertex(fvec3(0, -0.2f, 0.8f));
        ms->addIndices(0, 7);
        ms->addIndices(7, 8);
        ms->addIndices(7, 9);

        return ms;
    }

}