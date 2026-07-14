#pragma once

#include <string>

class cModel;

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

class theMeshView {
    cModel* m_model = nullptr;
    view_state m_view_state;                         ///< Current visualization state (e.g., whether to show curvature or normals)
    void init();
    void create_model_view();
    void createVertexK1View();
    void createVertexK2View();
    void createFaceNormalsView();
    void createVertexNormalsView();
    void set_callbacks();
    void reset_view();

#define MAP_TYPE_GAUSSIAN_CURVATURE 1
#define MAP_TYPE_MEAN_CURVATURE 2
#define MAP_TYPE_K1_K2_PREVIEW 3
    bool toggle_map_view(int map_type);
    void print_info();

public:
    theMeshView(cModel* model = nullptr);

    void set_model(cModel* model) { m_model = model; }
    void load_model(const std::string& filename);
    void render();

    void flip_mesh();
    void calculate_curvatures();
    void toggle_show_mesh();
    void toggle_show_face_normals();
    void toggle_show_vertex_normals();
    void toggle_show_principal_curvatures();
    void toggle_show_gaussian_curvatures();
    void toggle_show_mean_curvatures();
    void toggle_show_k1_k2_preview();
};
