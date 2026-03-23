#include "gl_graphics.h"
#include "mesh.h"
#include "gl_prim.h"

namespace base_opengl {
    /**
     * @brief Creates a gl_prim object from the provided mesh data.
     * Initializes OpenGL buffers for vertex positions, normals, and indices.
     * @param mesh Pointer to mesh_data containing geometry information.
     * @param drmode OpenGL draw mode (default: GL_FILL).
     * @param dr_el Whether to use element drawing (default: true).
     */
    void gl_prim::create_from_mesh(mesh_data* mesh, GLenum drmode /*= GL_FILL*/, bool dr_el /*= true*/) {
        m_mesh_data = *mesh;

        if (m_mesh_data.num_vertices == 0)
            return;

        draw_elements = dr_el;
        draw_mode = drmode;

        int idx = 0;

        glGenVertexArrays(1, &vao);
        glBindVertexArray(vao);

        GLuint position_buffer;
        glGenBuffers(1, &position_buffer);
        glBindBuffer(GL_ARRAY_BUFFER, position_buffer);
        glBufferData(GL_ARRAY_BUFFER, m_mesh_data.num_vertices * sizeof(float), &m_mesh_data.vertices[0], GL_STATIC_DRAW);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, (void*)0);
        glEnableVertexAttribArray(idx);

        if (m_mesh_data.num_normals) {
            GLuint normals_buffer;
            glGenBuffers(1, &normals_buffer);
            glBindBuffer(GL_ARRAY_BUFFER, normals_buffer);
            glBufferData(GL_ARRAY_BUFFER, m_mesh_data.num_normals * sizeof(float), &m_mesh_data.normals[0], GL_STATIC_DRAW);
            ++idx;
            glVertexAttribPointer(idx, 3, GL_FLOAT, GL_FALSE, 0, (void*)0);
            glEnableVertexAttribArray(idx);
        }

        if (m_mesh_data.num_curvatures > 0)
        {
            GLuint color_buffer;
            glGenBuffers(1, &color_buffer);
            glBindBuffer(GL_ARRAY_BUFFER, color_buffer);
            glBufferData(GL_ARRAY_BUFFER, m_mesh_data.num_curvatures * sizeof(float), &m_mesh_data.curvatures[0], GL_STATIC_DRAW);
            ++idx;
            glVertexAttribPointer(idx, 3, GL_FLOAT, GL_FALSE, 0, (void*)0);
            glEnableVertexAttribArray(idx);
            use_vertex_color = 1;
        }


        GLuint index_buffer;
        glGenBuffers(1, &index_buffer);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, index_buffer);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, m_mesh_data.num_indices * sizeof(unsigned int), &m_mesh_data.indices[0], GL_STATIC_DRAW);
        glBindVertexArray(0);
    }

    /**
     * @brief Creates a gl_prim from a half_edge_mesh.
     * Converts the half_edge_mesh to mesh_data and initializes a gl_prim object.
     * @param ms Pointer to half_edge_mesh.
     * @param drmode OpenGL draw mode (default: GL_LINE).
     * @param dr_el Whether to use element drawing (default: true).
     * @return Pointer to the created gl_prim.
     */
    // gl_prim* create_prim(half_edge_mesh<double>* ms, GLenum drmode/*=GL_LINE*/, bool dr_el/*=true*/) {
    //     if (!ms) return nullptr;
    //     mesh_data mdata;
    //     collect_mesh_data<double>(ms, mdata);
    //     gl_prim* prim = new gl_prim;
    //     prim->create_from_mesh(&mdata, drmode);
    //     prim->set_draw_mode(drmode);
    //     return prim;
    // }

    /**
     * @brief Creates a gl_prim from a gl_mesh.
     * Converts the gl_mesh to mesh_data and initializes a gl_prim object.
     * @param ms Pointer to gl_mesh.
     * @param drmode OpenGL draw mode (default: GL_LINE).
     * @param dr_el Whether to use element drawing (default: true).
     * @return Pointer to the created gl_prim.
     */
    gl_prim* create_prim(gl_mesh* ms, GLenum drmode/*=GL_LINE*/, bool dr_el/*=true*/) {
        if (!ms) return nullptr;
        mesh_data mdata;
        collect_mesh_data(ms, mdata);
        gl_prim* prim = new gl_prim;
        prim->create_from_mesh(&mdata, drmode);
        prim->set_draw_mode(drmode);
        return prim;
    }
#if 0
    /**
     * @brief Creates a unit cube primitive.
     * Generates a half_edge_mesh representing a unit cube and returns a gl_prim.
     * @param drmode OpenGL draw mode (default: GL_LINE).
     * @param dr_el Whether to use element drawing (default: true).
     * @return Pointer to the created gl_prim.
     */
    gl_prim* create_cube(GLenum drmode/*=GL_LINE*/, bool dr_el/*=true*/) {
        std::unique_ptr<half_edge_mesh<double>> ms(create_unit_cube<double, double>());
        gl_prim* p = create_prim(ms.get(), drmode, dr_el);
        return p;
    }

    /**
     * @brief Creates a unit sphere primitive.
     * Generates a half_edge_mesh representing a unit sphere and returns a gl_prim.
     * @param drmode OpenGL draw mode (default: GL_LINE).
     * @param dr_el Whether to use element drawing (default: true).
     * @return Pointer to the created gl_prim.
     */
    gl_prim* create_sphere(GLenum drmode/*=GL_LINE*/, bool dr_el/*=true*/) {
        std::unique_ptr<half_edge_mesh<double>> ms(create_unit_sphere<double, double>());
        gl_prim* p = create_prim(ms.get(), drmode, dr_el);
        return p;
    }

    /**
     * @brief Creates a unit cylinder primitive.
     * Generates a half_edge_mesh representing a unit cylinder and returns a gl_prim.
     * @param drmode OpenGL draw mode (default: GL_LINE).
     * @param dr_el Whether to use element drawing (default: true).
     * @return Pointer to the created gl_prim.
     */
    gl_prim* create_cylinder(GLenum drmode/*=GL_LINE*/, bool dr_el/*=true*/) {
        std::unique_ptr<half_edge_mesh<double>> ms(create_unit_cylinder<double, double>());
        gl_prim* p = create_prim(ms.get(), drmode, dr_el);
        return p;
    }

    /**
     * @brief Creates a unit cone primitive.
     * Generates a half_edge_mesh representing a unit cone and returns a gl_prim.
     * @param drmode OpenGL draw mode (default: GL_LINE).
     * @param dr_el Whether to use element drawing (default: true).
     * @return Pointer to the created gl_prim.
     */
    gl_prim* create_cone(GLenum drmode/*=GL_LINE*/, bool dr_el/*=true*/) {
        std::unique_ptr<half_edge_mesh<double>> ms(create_unit_cone<double, double>());
        gl_prim* p = create_prim(ms.get(), drmode, dr_el);
        return p;
    }

    /**
     * @brief Creates a unit dodecahedron primitive.
     * Generates a half_edge_mesh representing a unit dodecahedron and returns a gl_prim.
     * @param drmode OpenGL draw mode (default: GL_LINE).
     * @param dr_el Whether to use element drawing (default: true).
     * @return Pointer to the created gl_prim.
     */
    gl_prim* create_dodecahedron(GLenum drmode/*=GL_LINE*/, bool dr_el/*=true*/) {
        std::unique_ptr<half_edge_mesh<double>> ms(create_unit_dodecahedron<double, double>());
        gl_prim* p = create_prim(ms.get(), drmode, dr_el);
        return p;
    }

    /**
     * @brief Creates a unit icosahedron primitive.
     * Generates a half_edge_mesh representing a unit icosahedron and returns a gl_prim.
     * @param drmode OpenGL draw mode (default: GL_LINE).
     * @param dr_el Whether to use element drawing (default: true).
     * @return Pointer to the created gl_prim.
     */
    gl_prim* create_icosahedron(GLenum drmode/*=GL_LINE*/, bool dr_el/*=true*/) {
        std::unique_ptr<half_edge_mesh<double>> ms(create_unit_icosahedron<double, double>());
        gl_prim* p = create_prim(ms.get(), drmode, dr_el);
        return p;
    }

    /**
     * @brief Creates a unit octahedron primitive.
     * Generates a half_edge_mesh representing a unit octahedron and returns a gl_prim.
     * @param drmode OpenGL draw mode (default: GL_LINE).
     * @param dr_el Whether to use element drawing (default: true).
     * @return Pointer to the created gl_prim.
     */
    gl_prim* create_octa(GLenum drmode/*=GL_LINE*/, bool dr_el/*=true*/) {
        std::unique_ptr<half_edge_mesh<double>> ms(create_unit_octa<double, double>());
        gl_prim* p = create_prim(ms.get(), drmode, dr_el);
        return p;
    }

    /**
     * @brief Creates a unit pentahedron primitive.
     * Generates a half_edge_mesh representing a unit pentahedron and returns a gl_prim.
     * @param drmode OpenGL draw mode (default: GL_LINE).
     * @param dr_el Whether to use element drawing (default: true).
     * @return Pointer to the created gl_prim.
     */
    gl_prim* create_penta(GLenum drmode/*=GL_LINE*/, bool dr_el/*=true*/) {
        std::unique_ptr<half_edge_mesh<double>> ms(create_unit_penta<double, double>());
        gl_prim* p = create_prim(ms.get(), drmode, dr_el);
        return p;
    }

    /**
     * @brief Creates a unit plane primitive.
     * Generates a half_edge_mesh representing a unit plane and returns a gl_prim.
     * @param drmode OpenGL draw mode (default: GL_LINE).
     * @param dr_el Whether to use element drawing (default: true).
     * @return Pointer to the created gl_prim.
     */
    gl_prim* create_plane(GLenum drmode/*=GL_LINE*/, bool dr_el/*=true*/) {
        std::unique_ptr<half_edge_mesh<double>> ms(create_unit_plane<double, double>());
        gl_prim* p = create_prim(ms.get(), drmode, dr_el);
        return p;
    }

    /**
     * @brief Creates a unit tetrahedron primitive.
     * Generates a half_edge_mesh representing a unit tetrahedron and returns a gl_prim.
     * @param drmode OpenGL draw mode (default: GL_LINE).
     * @param dr_el Whether to use element drawing (default: true).
     * @return Pointer to the created gl_prim.
     */
    gl_prim* create_tetra(GLenum drmode/*=GL_LINE*/, bool dr_el/*=true*/) {
        std::unique_ptr<half_edge_mesh<double>> ms(create_unit_tetra<double, double>());
        gl_prim* p = create_prim(ms.get(), drmode, dr_el);
        return p;
    }

    /**
     * @brief Creates a unit torus primitive.
     * Generates a half_edge_mesh representing a unit torus and returns a gl_prim.
     * @param drmode OpenGL draw mode (default: GL_LINE).
     * @param dr_el Whether to use element drawing (default: true).
     * @return Pointer to the created gl_prim.
     */
    gl_prim* create_torus(GLenum drmode/*=GL_LINE*/, bool dr_el/*=true*/) {
        std::unique_ptr<half_edge_mesh<double>> ms(create_unit_torus<double, double>());
        gl_prim* p = create_prim(ms.get(), drmode, dr_el);
        return p;
    }  
#endif

    /**
     * @class gl_ucs
     * @brief Represents a UCS (Universal Coordinate System) primitive for rendering axes.
     * Inherits from gl_prim.
     */
    class gl_ucs :public gl_prim {
    protected:
    public:
        /**
         * @brief Default constructor for gl_ucs.
         */
        gl_ucs() {
        }
        /**
         * @brief Destructor for gl_ucs.
         */
        virtual ~gl_ucs() {
        }
        /**
         * @brief Renders the UCS axes using the provided shader.
         * Draws three colored axes (X: red, Y: green, Z: blue) using GL_LINES.
         * @param _shader Pointer to the shader used for rendering.
         */
        virtual void render(gl_shader* _shader) {
            if (!vao) return;
            _shader->set_int("object_or_vertex_color", 0);

            // position object
            fmat4 ob_matrix = tmat * rmat * smat;
            ob_matrix = ob_matrix.transpose();    // convert for OpenGL!
            ob_matrix = view_matrix * ob_matrix;

            // pass transformation to shader
            _shader->set_mat4("model", ob_matrix);

            glBindVertexArray(vao);
            unsigned int point_count = (unsigned int)m_mesh_data.num_indices / 3;
            _shader->set_vec4("object_color", fvec4(1, 0, 0, 1));
            glDrawElements(GL_LINES, point_count, GL_UNSIGNED_INT, 0);
            _shader->set_vec4("object_color", fvec4(0, 1, 0, 1));
            glDrawElements(GL_LINES, point_count, GL_UNSIGNED_INT, (const void*)(point_count * sizeof(unsigned int)));
            _shader->set_vec4("object_color", fvec4(0, 0, 1, 1));
            glDrawElements(GL_LINES, point_count, GL_UNSIGNED_INT, (const void*)(2 * point_count * sizeof(unsigned int)));
            glBindVertexArray(0);
        }
    };

    /**
     * @brief Creates a UCS (Universal Coordinate System) primitive.
     * Generates a gl_mesh representing the UCS and returns a gl_ucs object.
     * @param drmode OpenGL draw mode (default: GL_FILL).
     * @param dr_el Whether to use element drawing (default: true).
     * @return Pointer to the created gl_prim (as gl_ucs).
     */
    gl_prim* create_UCS(GLenum drmode /*= GL_FILL*/, bool dr_el /*= true*/) {
        std::unique_ptr<gl_mesh> ms(create_UCS_mesh());
        mesh_data mdata;
        collect_mesh_data(ms.get(), mdata);
        gl_prim* p = new gl_ucs();
        p->create_from_mesh(&mdata, drmode, dr_el);
        p->set_draw_mode(drmode);
        return p;
    }
}
