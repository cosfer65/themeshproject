#include "framework.h"
#include "resource.h"
#include "graphics_test.h"
#include "arcball.h"

#undef min
#undef max

#include "gl_camera.h"
#include "gl_prim.h"
#include "gl_light.h"

using namespace base_math;
using namespace base_opengl;
using namespace base_app;

#define MAX_LOADSTRING 100

/**
 * @brief Dialog procedure for the About dialog box.
 *
 * Handles initialization and command messages for the About dialog.
 *
 * @param hDlg Handle to the dialog box.
 * @param message Message identifier.
 * @param wParam Additional message information.
 * @param lParam Additional message information.
 * @return TRUE if message was processed, FALSE otherwise.
 */
INT_PTR CALLBACK About(HWND hDlg, UINT message, WPARAM wParam, LPARAM lParam)
{
    UNREFERENCED_PARAMETER(lParam);
    switch (message)
    {
    case WM_INITDIALOG:
        return (INT_PTR)TRUE;

    case WM_COMMAND:
        if (LOWORD(wParam) == IDOK || LOWORD(wParam) == IDCANCEL)
        {
            EndDialog(hDlg, LOWORD(wParam));
            return (INT_PTR)TRUE;
        }
        break;
    }
    return (INT_PTR)FALSE;
}

/**
 * @class graphics_test
 * @brief Main application class for the graphics test.
 *
 * Inherits from gl_application and manages the window, camera, viewport, shaders, and primitives.
 */
class graphics_test :public gl_application {
    char m_szTitle[MAX_LOADSTRING];                  ///< The title bar text

    // std::unique_ptr<gl_viewport> m_view;            ///< Viewport for rendering
    std::unique_ptr<gl_camera> m_cam;               ///< gl_camera for the scene
    std::unique_ptr<gl_shader> m_shader;            ///< Shader used for rendering

    // Primitives to render
    std::unique_ptr<gl_prim> m_ucs;
    std::unique_ptr<gl_prim> m_cone;
    std::unique_ptr<gl_prim> m_cube;
    std::unique_ptr<gl_prim> m_cylinder;
    std::unique_ptr<gl_prim> m_dodeca;
    std::unique_ptr<gl_prim> m_icosa;
    std::unique_ptr<gl_prim> m_octa;
    std::unique_ptr<gl_prim> m_penta;
    std::unique_ptr<gl_prim> m_plane;
    std::unique_ptr<gl_prim> m_sphere;
    std::unique_ptr<gl_prim> m_tetra;
    std::unique_ptr<gl_prim> m_torus;

    std::unique_ptr < gl_light> m_light;

    std::unique_ptr<arcball> m_arcball;           ///< Arcball for mouse interaction
public:
    /**
     * @brief Constructor. Initializes the viewport.
     */
    graphics_test() {
        m_cam.reset(new gl_camera(fvec3(0, 0, 50), fvec3(0, 0, 0), fvec3(0, 1, 0)));
        //m_view.reset(new gl_viewport());
        m_arcball.reset(new arcball(1200, 700));
    }

    /**
     * @brief Pre-creates the application window and sets window properties.
     *
     * @param hInstance Handle to the application instance.
     * @return 1 on success.
     */
    virtual int precreate_window(HINSTANCE hInstance) {
        char m_szWindowClass[MAX_LOADSTRING];            // the main window class name
        LoadString(hInstance, IDS_APP_TITLE, m_szTitle, MAX_LOADSTRING);
        LoadString(hInstance, IDC_GRAPHICSTEST, m_szWindowClass, MAX_LOADSTRING);
        szWindowClass = m_szWindowClass;

        m_window.szTitle = m_szTitle;
        m_window.prefered_width = 1200;
        m_window.prefered_height = 700;

        gl_application::precreate_window(hInstance);

        m_wcex.hIcon = LoadIcon(hInstance, MAKEINTRESOURCE(IDI_GRAPHICSTEST));
        m_wcex.hIconSm = LoadIcon(hInstance, MAKEINTRESOURCE(IDI_SMALL));
        m_wcex.hCursor = LoadCursor(NULL, IDC_ARROW);
        m_wcex.hbrBackground = (HBRUSH)(COLOR_WINDOW + 1);
        m_wcex.lpszMenuName = MAKEINTRESOURCE(IDC_GRAPHICSTEST);
        m_wcex.lpszClassName = szWindowClass.c_str();

        return 1;
    }

    /**
     * @brief Destructor.
     */
    virtual ~graphics_test() {
    }

    /**
     * @brief Initializes the application, camera, shader, and primitives.
     *
     * Sets up the field of view, camera, loads shaders, creates and positions primitives, and configures OpenGL state.
     * @return 1 on success.
     */
    virtual int init_application() {
        // zoom is the angle of the field of view
        m_cam->set_fov(dtr<float>(20));

        // create the basic shader
        m_shader.reset(new gl_shader);
        m_shader->add_file(GL_VERTEX_SHADER, "resources/shaders/mesh_tools_VertexShader.glsl");
        m_shader->add_file(GL_FRAGMENT_SHADER, "resources/shaders/mesh_tools_FragmentShader.glsl");
        m_shader->load();

#if 1
        m_light.reset(new gl_light(gl_light::SPOTLIGHT));
        //m_light->set_position(fvec3(-100, 0, 100));

        m_light->set_position(fvec3(-20, 20, 20));
        m_light->set_ambient(fvec3(0.75f));
        m_light->set_diffuse(fvec3(0.5f));
        m_light->set_specular(fvec3(0.1f));
#else
        m_light.reset(new gl_light(gl_light::DIRLIGHT));
        m_light->set_direction(fvec3(20, 30, 40));
#endif

        m_ucs.reset(create_UCS());
        m_ucs->move_to(-6, 4, 0);
        m_ucs->set_scale(fvec3(1.5));

        m_cone.reset(create_cone<double, double>(GL_FILL));
        m_cone->move_to(fvec3(-2, 4, 0));
        m_cone->set_scale(fvec3(1.5f));

        m_cube.reset(create_cube<double, double>(GL_FILL));
        m_cube->move_to(fvec3(2, 4, 0));
        m_cube->set_scale(fvec3(1.5f));

        m_cylinder.reset(create_cylinder<double, double>(GL_FILL));
        m_cylinder->move_to(fvec3(6, 4, 0));
        m_cylinder->set_scale(fvec3(1.5f));

        m_dodeca.reset(create_dodecahedron<double, double>(GL_FILL));
        m_dodeca->move_to(fvec3(-6, 0, 0));
        m_dodeca->set_scale(fvec3(1.5f));

        m_icosa.reset(create_icosahedron<double, double>(GL_FILL));
        m_icosa->move_to(fvec3(-2, 0, 0));
        m_icosa->set_scale(fvec3(1.5f));

        m_octa.reset(create_octa<double, double>(GL_FILL));
        m_octa->move_to(fvec3(2, 0, 0));
        m_octa->set_scale(fvec3(1.5f));

        m_penta.reset(create_penta<double, double>());
        m_penta->move_to(fvec3(6, 0, 0));
        m_penta->set_scale(fvec3(1.5f));

        m_plane.reset(create_plane<double, double>());
        m_plane->move_to(fvec3(-6, -4, 0));
        m_plane->set_scale(fvec3(1.5f));

        m_sphere.reset(create_sphere<double, double>(GL_FILL));
        m_sphere->move_to(fvec3(-2, -4, 0));
        m_sphere->set_scale(fvec3(1.5f));

        m_tetra.reset(create_tetra<double, double>());
        m_tetra->move_to(fvec3(2, -4, 0));
        m_tetra->set_scale(fvec3(1.5f));

        m_torus.reset(create_torus<double, double>());
        m_torus->move_to(fvec3(6, -4, 0));
        m_torus->set_scale(fvec3(1.5f));

        // OpenGL initialization we want for this sample
        //glEnable(GL_CULL_FACE);
        glEnable(GL_DEPTH_TEST);
        glDepthFunc(GL_LEQUAL);
        glEnable(GL_MULTISAMPLE);

        return 1;
    }

    /**
     * @brief Advances the simulation by the elapsed time.
     *
     * @param fElapsed Time elapsed since last step, in seconds.
     */
    virtual void step_simulation(float fElapsed) {
    }

    /**
     * @brief Renders the scene.
     *
     * Sets up the viewport, clears the screen, sets up the camera and shader, and renders all primitives.
     */
    virtual void render() {
        // set the viewport to the whole window
        m_cam->set_viewport();

        // clear screen
        glClearColor(.25f, .25f, .25f, 1.f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        // Convert the arcball's current rotation quaternion into a rotation matrix
        fmat4 rot_mat(m_arcball->rotation());

        // Update the view matrix of all primitives to match the arcball's rotation
        m_ucs->view_matrix = rot_mat;
        m_cone->view_matrix = rot_mat;
        m_cube->view_matrix = rot_mat;
        m_cylinder->view_matrix = rot_mat;
        m_dodeca->view_matrix = rot_mat;
        m_icosa->view_matrix = rot_mat;
        m_octa->view_matrix = rot_mat;
        m_penta->view_matrix = rot_mat;
        m_plane->view_matrix = rot_mat;
        m_sphere->view_matrix = rot_mat;
        m_tetra->view_matrix = rot_mat;
        m_torus->view_matrix = rot_mat;

        // combine the view and camera matrices into one
        // these matrices are COLUMN MAJOR!
        fmat4 cam_matrix = m_cam->perspective();

        // enable the shader
        m_shader->use();
        // set the combined view matrix
        m_shader->set_mat4("camera", cam_matrix);


        // UCS takes care of its own color
        m_ucs->render(m_shader.get());
        m_light->apply(m_shader.get());

        // render the objects using color
        m_shader->set_int("use_color_or_texture", 1);
        m_shader->set_vec4("object_color", fvec4(.7f, .8f, .9f, 1.f));

        m_cone->render(m_shader.get());
        m_cube->render(m_shader.get());
        m_cylinder->render(m_shader.get());
        m_dodeca->render(m_shader.get());
        m_icosa->render(m_shader.get());
        m_octa->render(m_shader.get());
        m_penta->render(m_shader.get());
        m_plane->render(m_shader.get());
        m_sphere->render(m_shader.get());
        m_tetra->render(m_shader.get());
        m_torus->render(m_shader.get());

        m_shader->end();
    }

    /**
     * @brief Handles window resize events.
     *
     * Updates the viewport aspect ratio.
     *
     * @param width New window width.
     * @param height New window height.
     */
    virtual void resize_window(int width, int height) {
        m_cam->set_aspect(width, height);
        // m_view->set_window_aspect(width, height);
        m_arcball->resize((float)width, (float)height);
    }

    /**
     * @brief Handles command messages from the window.
     *
     * @param cmd Command identifier.
     * @return 1 if handled.
     */
    virtual int onCommand(int cmd) {
        switch (cmd) {
        case IDM_EXIT:
            PostQuitMessage(0);
            break;
        case IDM_ABOUT:
            DialogBox(hInst, MAKEINTRESOURCE(IDD_ABOUTBOX), m_window.hWnd, About);
            break;
        }
        return 1;
    }

    /// @brief Handles mouse move events while interacting with the mesh view.
    ///
    /// Updates the arcball controller with the current cursor position to continuously
    /// update the view rotation during an active drag operation.
    ///
    /// @param x     Current mouse X position in window/client coordinates.
    /// @param y     Current mouse Y position in window/client coordinates.
    /// @param extra Additional mouse state flags (currently unused).
    void onMouseMove(int x, int y, unsigned __int64 extra) {
        m_arcball->drag(float(x), float(y));
    }

    /// @brief Handles left mouse button press to start an arcball rotation.
    ///
    /// Captures the initial mouse position and notifies the arcball controller that
    /// a drag operation has started, enabling interactive rotation of the mesh view.
    ///
    /// @param x     Mouse X position at the time of button press.
    /// @param y     Mouse Y position at the time of button press.
    /// @param extra Additional mouse state flags (currently unused).
    void onLMouseDown(int x, int y, unsigned __int64 extra) {
        m_arcball->beginDrag(float(x), float(y));
    }

    /// @brief Handles left mouse button release to end an arcball rotation.
    ///
    /// Signals the arcball controller to finalize the current drag operation, freezing
    /// the current view rotation until a new drag is started.
    ///
    /// @param x     Mouse X position at the time of button release (unused).
    /// @param y     Mouse Y position at the time of button release (unused).
    /// @param extra Additional mouse state flags (currently unused).
    void onLMouseUp(int x, int y, unsigned __int64 extra) {
        m_arcball->endDrag();
    }

    /**
     * @brief Handles mouse wheel events.
     *
     * @param delta Amount of wheel movement.
     * @param extra_btn Mouse button state.
     */
    virtual void onMouseWheel(int delta, WPARAM extra_btn) {
    }
};

/// Global instance of the graphics_test application.
graphics_test my_app;