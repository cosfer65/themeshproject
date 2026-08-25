#pragma once

#include "prim.h"

class draw_object : public btm::gl_prim {
    int do_vertices = 0;
    int do_vectors = 0;
public:
    draw_object() {
        draw_type = GL_LINES;
        draw_mode = GL_LINE;
        draw_elements = true;
    }

    void add_vertex(const btm::fvec3& v) {
        int idx = m_mesh_data.add_vertex(v);
        ++do_vertices;
        m_mesh_data.add_index(idx);
    }   
    void add_vector(const btm::fvec3& vstart, const btm::fvec3& vend) {
        unsigned int idx_start = (unsigned int)m_mesh_data.num_vertices / 3;
        m_mesh_data.add_vertex(vstart);
        m_mesh_data.add_vertex(vend);
        m_mesh_data.add_indices(idx_start, idx_start + 1);
        ++do_vectors;
    }
    void create_prim(GLenum drmode=GL_LINE) {
        create_from_mesh(&m_mesh_data, drmode);
    }

    int num_do_vertices() const { return do_vertices; }
    int num_do_vectors() const { return do_vectors; }
    void clear_do() { do_vertices = 0; do_vectors = 0; }
};
