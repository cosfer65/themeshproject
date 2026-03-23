#pragma once

#include "matrix.h"
#include "gl_shaders.h"
#include "gl_window.h"

class cModel;
struct mesh_view_private;

struct view_state {
    bool show_wireframe = false;  ///< Whether to render the mesh in wireframe mode.
    bool show_curvature = false; ///< Whether to visualize curvature directions.
    bool show_normals = false;   ///< Whether to visualize face normals.
    bool show_curvature_map = false;   ///< Whether to visualize the curvature map.

    void reset() {
        show_curvature = false; // disable curvature visualization until it is regenerated for the new model
        show_normals = false;   // disable normals visualization until it is regenerated for the new model
        show_wireframe = false;  // disable wireframe visualization until the user toggles it
        show_curvature_map = false;  // disable map visualization until the user toggles it
    }
};

class meshViewWindow : public glViewWindow {
    cModel* m_model = nullptr;
    mesh_view_private* m_private;
    view_state m_view_state;                         ///< Current visualization state (e.g., whether to show curvature or normals)

    void build_model_representation();
    void build_normals_representation();
    void generate_model_curvature_view();
    void generate_model_curvature_map_view();
    bool init_done = false;
    bool dragging = false;                           ///< Indicates whether the user is currently dragging with the mouse (for panning)
    int last_mouse_x = 0;                            ///< Last recorded mouse X position (used for calculating deltas during dragging)
    int last_mouse_y = 0;                            ///< Last recorded mouse Y position (used for calculating deltas during dragging)
public:
    meshViewWindow();
    virtual ~meshViewWindow();

    void init();

    void set_model(cModel* model) {
        m_model = model;
    }
    void reset_view();

    void step_simulation(float fElapsed) {};
    void render_scene(base_math::fmat4& cam_matrix, base_math::fmat4& rot_mat, base_opengl::gl_shader* shdr);
    virtual void render();
    virtual int onCommand(int cmd);
    virtual LRESULT OnSize(int cx, int cy);

    virtual void onMouseWheel(int delta, unsigned __int64 extra_btn);
    virtual void onMouseMove(int x, int y, unsigned __int64 extra);
    virtual void onLMouseDown(int x, int y, unsigned __int64 extra);
    virtual void onLMouseUp(int x, int y, unsigned __int64 extra);
    virtual void onRMouseDown(int x, int y, unsigned __int64 extra);
    virtual void onRMouseUp(int x, int y, unsigned __int64 extra);
};
