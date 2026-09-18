#pragma once

#include "prim.h"

class draw_object : public btm::gl_prim {
    btm::mesh_data m_mesh_data;
public:
    draw_object() {
        draw_type = GL_LINES;
        draw_mode = GL_LINE;
        draw_elements = true;
    }

    void add_vertex(const btm::fvec3& v) {
        int idx = m_mesh_data.add_vertex(v);
        m_mesh_data.add_index(idx);
    }   
    void add_vector(const btm::fvec3& vstart, const btm::fvec3& vend) {
        unsigned int idx_start = (unsigned int)m_mesh_data.num_vertices;// / 3;
        m_mesh_data.add_vertex(vstart);
        m_mesh_data.add_vertex(vend);
        m_mesh_data.add_indices(idx_start, idx_start + 1);
    }
    void create_prim(GLenum drmode=GL_LINE) {
        create_from_mesh(&m_mesh_data, drmode);
    }
};
