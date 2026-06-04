#include "application.h"
#include "glew.h"
#include "timer.h"

#include "window.h"
#include "gl_context.h"

#include "camera.h"
#include "shaders.h"
#include "light.h"

#include "mesh_explicit.h"
#include "mesh_renderer.h"

#include <vector>

using namespace btm;

struct render_resources {
    std::unique_ptr<gl_camera> g_cam;     ///< gl_camera used to view the scene and compute view/projection.
    std::unique_ptr<gl_shader> g_shader;  ///< Shader program used for mesh and helper rendering.
};

static bool create_mesh(MeshExplicit<float>& r_mesh) {
    // vertices of a cube
    static std::vector<basevec3<float>> cube_vertices = {
        {0.5f, 0.5f, -0.5f},   {0.5f, -0.5f, -0.5f},  {0.5f, 0.5f, 0.5f},
        {0.5f, -0.5f, 0.5f},   {-0.5f, 0.5f, -0.5f}, {-0.5f, -0.5f, -0.5f},
        {-0.5f, 0.5f, 0.5f},  {-0.5f, -0.5f, 0.5f},
    };
    // triangles of the cube (12 triangles, 2 per face)
    static std::vector<Triangle> cube_triangles = {
        {4, 2, 0}, {2, 7, 3}, {6, 5, 7}, {1, 7, 5}, {0, 3, 1}, {4, 1, 5},
        {4, 6, 2}, {2, 6, 7}, {6, 4, 5}, {1, 3, 7}, {0, 2, 3}, {4, 0, 1}
    };

    // Add vertices and triangles to the mesh. The MeshExplicit class will store them in its internal data structures.
    for (const auto &v : cube_vertices) {
        r_mesh.add_vertex(v);
    }
    for (const auto &tri : cube_triangles) {
        r_mesh.add_triangle(tri);
    }

    // build mesh adjacency and attributes. The adjacency will be used to compute the vertex normals and other attributes for rendering.
    r_mesh.build_adjacency();
    r_mesh.build_attributes();

    return true;
}

void render(render_resources& resources, MeshRenderer<float>* mesh_renderer) {
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
    resources.g_cam->set_aspect(width, height);
    resources.g_cam->set_viewport();
    // We then compute the camera matrix (combination of view and projection) and pass it to the shader. 
    fmat4 cam_matrix = resources.g_cam->perspective();
    // bring the shader into use and set the camera matrix uniform. 
    // The shader will use this matrix to transform the vertex positions from world space to clip space for rendering.
    resources.g_shader->use();
    resources.g_shader->set_mat4("camera", cam_matrix);

    // now call the renderer to draw the mesh. The MeshRenderer will use the shader and the mesh data to issue 
    // OpenGL draw calls to render the triangles of the mesh on the screen.
    mesh_renderer->render(resources.g_shader.get());

    resources.g_shader->end();

    end_render();
}

void step_simulation(float fElapsed, MeshRenderer<float>* mesh_renderer)
{
    if (mesh_renderer) {
        mesh_renderer->rotate_by(dtr(fElapsed * 10), dtr(fElapsed * 20), dtr(fElapsed * 30));
    }
}

int APIENTRY WinMain(HINSTANCE hInstance, HINSTANCE, LPTSTR, int nCmdShow) {
    application the_app;
    render_resources resources;

    init_framework();

    FrameWindow *pFrame = create_main_window(false, 800, 600, "Mesh Data Structures");

    MeshExplicit<float> mesh;
    create_mesh(mesh);
    // The MeshRenderer will generate the vertex buffer and other necessary data for rendering based on the mesh topology and attributes.
    std::unique_ptr<MeshRenderer<float>> mesh_renderer(new MeshRenderer<float>(mesh));

    // Initialize application resources (camera, shader)
    /// {
    resources.g_cam.reset(new btm::gl_camera(btm::fvec3(0, 0, 20), btm::fvec3(0, 0, 0), btm::fvec3(0, 1, 0)));
    resources.g_cam->set_fov(btm::dtr(10.f));

    resources.g_shader.reset(new gl_shader);
    resources.g_shader->add_file(GL_VERTEX_SHADER, "resources/shaders/generic_VertexShader.glsl");
    resources.g_shader->add_file(GL_FRAGMENT_SHADER, "resources/shaders/generic_FragmentShader.glsl");
    resources.g_shader->load();
    /// }

    start_timer();
    get_elapsed_time();

    while (pollEvents()) {
        float fElapsed = (float)get_elapsed_time();
        step_simulation(fElapsed, mesh_renderer.get());
        render(resources, mesh_renderer.get());

        {
            static int m_nFrames = 0; // frame Counter
            static float tot = 0;     // time couner
            tot += fElapsed;          // increment counter
            m_nFrames++;
            if (tot >= 1.f) // one second reached
            {
                char txt[200];
                sprintf_s(txt, "Mesh Data Structures, fps:%d", m_nFrames);
                HWND hWnd = pFrame->hWnd;
                SetWindowText(hWnd, txt);
                tot = 0; // reset counters
                m_nFrames = 0;
            }
        }
    }
    the_app.terminate();
    stop_timer();

    return 0;
}
