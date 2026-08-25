#pragma once

#include "camera.h"
#include "shaders.h"
#include "light.h"
#include "prim.h"
#include "arcball.h"
#include "draw_objects.h"
#include "ucs_view.h"
#include "font.h"

struct mesh_view_resources {
    mesh_view_resources() = default;
    gl_camera m_cam;                   ///< gl_camera used to view the scene and compute view/projection.
    gl_shader m_shader;                ///< Shader program used for mesh and helper rendering.
    gl_light m_light;                  ///< Primary scene light affecting shading.

    UCS_view m_ucs_view;               ///< UCS (user coordinate system) view used to render the axes widget.

    arcball m_arcball;                                 ///< Arcball for mouse interaction
    btm::gl_font* font2D;

    // {{ the loaded model and the acompanying visualization objects
    std::vector<std::unique_ptr<btm::gl_prim>> m_draw_parts; ///< Drawable mesh parts converted to OpenGL primitives.
    // }}

    draw_object m_face_normals;
    draw_object m_vertex_normals;
    draw_object m_model_curvatures_k1;                  ///< Visualization of minimum principal curvature directions per vertex.
    draw_object m_model_curvatures_k2;                  ///< Visualization of maximum principal curvature directions per vertex.

    draw_object m_feature_lines[3];            ///< Visualization of feature lines (e.g., sharp edges) in the mesh. 0-ridges, 1-valleys, 2-creases

    draw_object m_potential_ridges;           ///< Visualization of potential ridge feature points in the mesh.
    draw_object m_potential_valleys;          ///< Visualization of potential valley feature points in the mesh.    
    draw_object m_potential_creases;          ///< Visualization of potential crease feature edges in the mesh.
    draw_object m_boundary_edges;             ///< Visualization of boundary edges in the mesh.

    void invalidate_visualizations() {
        m_face_normals.clear();
        m_vertex_normals.clear();
        m_model_curvatures_k1.clear();
        m_model_curvatures_k2.clear();
        for (int i = 0; i < 3; ++i) {
            m_feature_lines[i].clear();
        }
        m_potential_ridges.clear();
        m_potential_valleys.clear();
        m_potential_creases.clear();
        m_boundary_edges.clear();
    }

    void init() {
        m_cam = gl_camera(fvec3(0, 0, 50), fvec3(0, 0, 0), fvec3(0, 1, 0));

        m_shader.add_file(GL_VERTEX_SHADER, "resources/shaders/meshprojectVertexShader.glsl");
        m_shader.add_file(GL_FRAGMENT_SHADER, "resources/shaders/meshprojectFragmentShader.glsl");
        m_shader.load();

        m_light.set_position(fvec3(-20, 20, 50));

        m_ucs_view.initialize();

        font2D = get_font_manager().create_font("Consolas", "Consolas", 12);
    }
};
