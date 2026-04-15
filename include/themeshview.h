#pragma once

#include "matrix.h"
#include "gl_shaders.h"
#include "gl_window.h"

class cModel;
struct mesh_view_private;

struct view_state {
    bool show_wireframe = false;  ///< Whether to render the mesh in wireframe mode.
    bool show_face_normals = false;   ///< Whether to visualize face normals.
    bool show_vertex_normals = false; ///< Whether to visualize vertex normals.
	bool show_principal_curvatures = false; ///< Whether to visualize principal curvature directions.
	bool show_gaussian_curvature = false; ///< Whether to visualize Gaussian curvature using a color map.
	bool show_mean_curvature = false; ///< Whether to visualize mean curvature using a color map.
	bool show_k1_k2_preview = false; ///< Whether to show the convex ellipsoid preview based on k1 and k2 directions.

    void reset() {
		// disable all visualization options by default
        show_wireframe = false;
        show_face_normals = false;
        show_vertex_normals = false;
        show_principal_curvatures = false;
        show_gaussian_curvature = false;
        show_mean_curvature = false;
        show_k1_k2_preview = false;
    }

    bool show_curvture_map() {
        return show_gaussian_curvature | show_mean_curvature | show_k1_k2_preview;
    }
};

class meshViewWindow : public glViewWindow {
    bool init_done = false;
    bool dragging = false;                           ///< Indicates whether the user is currently dragging with the mouse (for panning)
    int last_mouse_x = 0;                            ///< Last recorded mouse X position (used for calculating deltas during dragging)
    int last_mouse_y = 0;                            ///< Last recorded mouse Y position (used for calculating deltas during dragging)

    cModel* m_model = nullptr;
    mesh_view_private* m_private;
    view_state m_view_state;                         ///< Current visualization state (e.g., whether to show curvature or normals)

    void create_model_view();
    bool generate_face_normals_view();
    bool generate_vertex_normals_view();
    bool generate_principal_curvature_view();
    bool generate_gaussian_curvature_map_view();
    bool generate_k1_k2_map_view();
    bool generate_mean_curvature_map_view();

#define MAP_TYPE_GAUSSIAN_CURVATURE 1
#define MAP_TYPE_MEAN_CURVATURE 2
#define MAP_TYPE_K1_K2_PREVIEW 3
	bool toggle_map_view(int map_type);
    void print_info();
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
