#include "framework.h"
#include "resource.h"
#include "graphics_test.h"

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

    std::unique_ptr<gl_viewport> m_view;            ///< Viewport for rendering
    std::unique_ptr<gl_camera> m_cam;               ///< Camera for the scene
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
public:
    /**
     * @brief Constructor. Initializes the viewport.
     */
    graphics_test() {
        m_view.reset(new gl_viewport());
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
        m_view->set_fov(dtr<float>(20));
        // camera position, look at point, and up vector orientation
        m_cam.reset(new gl_camera(fvec3(0, 0, 50), fvec3(0, 0, 0), fvec3(0, 1, 0)));

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

        m_cone.reset(create_cone(GL_FILL));
        m_cone->move_to(fvec3(-2, 4, 0));
        m_cone->set_scale(fvec3(1.5f));

        m_cube.reset(create_cube(GL_FILL));
        m_cube->move_to(fvec3(2, 4, 0));
        m_cube->set_scale(fvec3(1.5f));

        m_cylinder.reset(create_cylinder(GL_FILL));
        m_cylinder->move_to(fvec3(6, 4, 0));
        m_cylinder->set_scale(fvec3(1.5f));

        m_dodeca.reset(create_dodecahedron(GL_FILL));
        m_dodeca->move_to(fvec3(-6, 0, 0));
        m_dodeca->set_scale(fvec3(1.5f));

        m_icosa.reset(create_icosahedron(GL_FILL));
        m_icosa->move_to(fvec3(-2, 0, 0));
        m_icosa->set_scale(fvec3(1.5f));

        m_octa.reset(create_octa(GL_FILL));
        m_octa->move_to(fvec3(2, 0, 0));
        m_octa->set_scale(fvec3(1.5f));

        m_penta.reset(create_penta());
        m_penta->move_to(fvec3(6, 0, 0));
        m_penta->set_scale(fvec3(1.5f));

        m_plane.reset(create_plane());
        m_plane->move_to(fvec3(-6, -4, 0));
        m_plane->set_scale(fvec3(1.5f));

        m_sphere.reset(create_sphere(GL_FILL));
        m_sphere->move_to(fvec3(-2, -4, 0));
        m_sphere->set_scale(fvec3(1.5f));

        m_tetra.reset(create_tetra());
        m_tetra->move_to(fvec3(2, -4, 0));
        m_tetra->set_scale(fvec3(1.5f));

        m_torus.reset(create_torus());
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
        m_view->set_viewport();

        // clear screen
        glClearColor(.25f, .25f, .25f, 1.f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        // combine the view and camera matrices into one
        // these matrices are COLUMN MAJOR!
        fmat4 cam_matrix = m_cam->perspective() * m_view->perspective();

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
        m_view->set_window_aspect(width, height);
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

    /**
     * @brief Handles mouse movement events.
     *
     * Rotates primitives based on mouse movement and button state.
     *
     * @param dx Change in x position.
     * @param dy Change in y position.
     * @param extra_btn Mouse button state.
     */
    virtual void onMouseMove(int dx, int dy, WPARAM extra_btn) {
        if (extra_btn & MK_LBUTTON) {
            m_ucs->rotate_by(dtr(dy / 5.f), dtr(dx / 5.f), 0);
            m_cone->rotate_by(dtr(dy / 5.f), dtr(dx / 5.f), 0);
            m_cube->rotate_by(dtr(dy / 5.f), dtr(dx / 5.f), 0);
            m_cylinder->rotate_by(dtr(dy / 5.f), dtr(dx / 5.f), 0);
            m_dodeca->rotate_by(dtr(dy / 5.f), dtr(dx / 5.f), 0);
            m_icosa->rotate_by(dtr(dy / 5.f), dtr(dx / 5.f), 0);
            m_octa->rotate_by(dtr(dy / 5.f), dtr(dx / 5.f), 0);
            m_penta->rotate_by(dtr(dy / 5.f), dtr(dx / 5.f), 0);
            m_plane->rotate_by(dtr(dy / 5.f), dtr(dx / 5.f), 0);
            m_sphere->rotate_by(dtr(dy / 5.f), dtr(dx / 5.f), 0);
            m_tetra->rotate_by(dtr(dy / 5.f), dtr(dx / 5.f), 0);
            m_torus->rotate_by(dtr(dy / 5.f), dtr(dx / 5.f), 0);
        }
        if (extra_btn & MK_RBUTTON) {
            m_ucs->rotate_by(0, 0, dtr(dy / 5.f));
            m_cone->rotate_by(0, 0, dtr(dy / 5.f));
            m_cube->rotate_by(0, 0, dtr(dy / 5.f));
            m_cylinder->rotate_by(0, 0, dtr(dy / 5.f));
            m_dodeca->rotate_by(0, 0, dtr(dy / 5.f));
            m_icosa->rotate_by(0, 0, dtr(dy / 5.f));
            m_octa->rotate_by(0, 0, dtr(dy / 5.f));
            m_penta->rotate_by(0, 0, dtr(dy / 5.f));
            m_plane->rotate_by(0, 0, dtr(dy / 5.f));
            m_sphere->rotate_by(0, 0, dtr(dy / 5.f));
            m_tetra->rotate_by(0, 0, dtr(dy / 5.f));
            m_torus->rotate_by(0, 0, dtr(dy / 5.f));
        }
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