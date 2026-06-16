#include "application.h"
#include "glew.h"
#include "timer.h"

#include "gl_context.h"

#include "camera.h"
#include "shaders.h"
#include "prim.h"

#include "mesh_explicit.h"

using namespace btm;

// global application resources 
std::unique_ptr<gl_camera> g_cam;     ///< gl_camera used to view the scene and compute view/projection.
std::unique_ptr<gl_shader> g_shader;  ///< Shader program used for mesh and helper rendering.

// gl_prim is a wrapper around OpenGL calls to render the mesh data.
std::unique_ptr<gl_prim> mesh_renderer;

static bool create_mesh(MeshExplicit<float>& r_mesh) {
    // vertices of a cube
    static std::vector<basevec3<float>> cube_vertices = {
        {0.5f, 0.5f, -0.5f},   {0.5f, -0.5f, -0.5f},  {0.5f, 0.5f, 0.5f},
        {0.5f, -0.5f, 0.5f},   {-0.5f, 0.5f, -0.5f}, {-0.5f, -0.5f, -0.5f},
        {-0.5f, 0.5f, 0.5f},  {-0.5f, -0.5f, 0.5f},
    };
    // triangles of the cube (12 triangles, 2 per face)
    static std::vector<FaceExplicit> cube_triangles = {
        {4, 2, 0}, {2, 7, 3}, {6, 5, 7}, {1, 7, 5}, {0, 3, 1}, {4, 1, 5},
        {4, 6, 2}, {2, 6, 7}, {6, 4, 5}, {1, 3, 7}, {0, 2, 3}, {4, 0, 1}
    };

    // Add vertices and triangles to the mesh. The MeshExplicit class will store them in its internal data structures.
    for (const auto &v : cube_vertices) {
        r_mesh.add_vertex(0, v);
    }
    for (const auto &tri : cube_triangles) {
        r_mesh.add_face(tri);
    }

    return true;
}

void render() {
    // first we check OpenGL state
    btm::GLContext* context = btm::get_current_gl_context();
    if (!context)
        return;

    begin_render();

    glClearColor(0.2f, 0.4f, 0.6f, 1.f);

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // gets the window dimensions from the GL context and sets the camera aspect ratio and viewport accordingly. 
    // This ensures that the rendered scene is displayed correctly within the window.
    int width = context->width();
    int height = context->height();
    g_cam->set_aspect(width, height);
    g_cam->set_viewport();

    // bring the shader into use 
    g_shader->use();
    // call the camera's apply method to set the view and projection matrices in the shader.
    // The shader will use these matrices to transform the vertex positions from world space to clip space for rendering.
    g_cam->apply(g_shader.get());

    // now call the renderer to draw the mesh. The MeshRenderer will use the shader and the mesh data to issue 
    // OpenGL draw calls to render the triangles of the mesh on the screen.
    mesh_renderer->force_black = false;  // set to true to render the mesh in black (e.g., for wireframe)
    mesh_renderer->render(g_shader.get());

    g_shader->end();

    end_render();
}

void step_simulation(float fElapsed)
{
    if (mesh_renderer) {
        mesh_renderer->rotate_by(dtr(fElapsed * 10), dtr(fElapsed * 20), dtr(fElapsed * 30));
    }
}

int APIENTRY WinMain(HINSTANCE hInstance, HINSTANCE, LPTSTR, int nCmdShow) {
    application the_app;

    init_framework();

    FrameWindow *pFrame = create_main_window(false, 800, 600, "Mesh Data Structures");

    MeshExplicit<float> mesh;
    create_mesh(mesh);
    // The mesh_renderer/gl_prim will generate the vertex buffer and other necessary data for rendering based on the mesh topology and attributes.
    mesh_renderer.reset(create_prim<float>(&mesh));

    // Initialize application resources (camera, shader)
    /// {
    g_cam.reset(new btm::gl_camera(btm::fvec3(0, 0, 20), btm::fvec3(0, 0, 0), btm::fvec3(0, 1, 0)));
    g_cam->set_fov(btm::dtr(10.f));

    g_shader.reset(new gl_shader);
    g_shader->add_file(GL_VERTEX_SHADER, "resources/shaders/generic_VertexShader.glsl");
    g_shader->add_file(GL_FRAGMENT_SHADER, "resources/shaders/generic_FragmentShader.glsl");
    g_shader->load();
    /// }

    start_timer();
    get_elapsed_time();

    while (pollEvents()) {
        float fElapsed = (float)get_elapsed_time();
        step_simulation(fElapsed);
        render();
    }
    the_app.terminate();
    stop_timer();

    return 0;
}
