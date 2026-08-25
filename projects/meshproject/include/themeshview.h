#pragma once

#include <string>

class cModel;

struct view_state {
    bool show_mesh = true;          ///< Whether to render the mesh surface.
    bool show_wireframe = false;  ///< Whether to render the mesh in wireframe mode.
    bool show_face_normals = false;   ///< Whether to visualize face normals.
    bool show_vertex_normals = false; ///< Whether to visualize vertex normals.
    
    bool show_principal_k1 = false; ///< Whether to visualize principal curvature directions.
    bool show_principal_k2 = false; ///< Whether to visualize principal curvature directions.

    bool show_gaussian_curvature = false; ///< Whether to visualize Gaussian curvature using a color map.
    bool show_mean_curvature = false; ///< Whether to visualize mean curvature using a color map.
    bool show_k1_k2_preview = false; ///< Whether to show the convex ellipsoid preview based on k1 and k2 directions.
    bool show_ridges = false; ///< Whether to visualize ridges on the mesh.
    bool show_valleys = false; ///< Whether to visualize valleys on the mesh.
    bool show_creases = false; ///< Whether to visualize creases on the mesh.

    void reset() {
        // enable mesh rendering by default
        show_mesh = true;
        // disable all visualization options by default
        show_wireframe = false;
        show_face_normals = false;
        show_vertex_normals = false;

        show_principal_k1 = false;
        show_principal_k2 = false;

        show_gaussian_curvature = false;
        show_mean_curvature = false;
        show_k1_k2_preview = false;
        show_ridges = false;
        show_valleys = false;
        show_creases = false;
    }

    bool show_curvture_map() {
        return show_gaussian_curvature | show_mean_curvature | show_k1_k2_preview;
    }
};

struct mesh_view_resources;

class theMeshView {
    mesh_view_resources* m_view_resources;
    cModel* m_model = nullptr;
    view_state m_view_state;                         ///< Current visualization state (e.g., whether to show curvature or normals)

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

    void load_model(const std::string& filename);
    void render();

    void flip_mesh();
    void calculate_curvatures();
    void calculate_feature_lines();

    void toggle_show_mesh();
    void toggle_show_wireframe();
    void toggle_show_face_normals();
    void toggle_show_vertex_normals();

    void toggle_show_principal_k1();
    void toggle_show_principal_k2();

    void toggle_show_gaussian_curvatures();
    void toggle_show_mean_curvatures();
    void toggle_show_k1_k2_preview();

    void toggle_show_ridges();
    void toggle_show_valleys();
    void toggle_show_creases();

    void reset_view_state();
#ifdef _DEBUG
    void test_function();
#endif
};
