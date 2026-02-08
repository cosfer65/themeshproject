/// @file mesh_tools.cpp
/// @brief Entry point and main application class for the Mesh Tools OpenGL-based viewer.
///
/// This file defines the `mesh_tools` application class, which derives from
/// `base_app::gl_application` and configures a Win32 window hosting an OpenGL
/// rendering context used by `mesh_view`. It also contains the About dialog
/// procedure and the global application instance.

#include "resource.h"
#include "mesh_tools.h"
#include "base_app.h"
#include "mesh_view.h"

/// Use the `base_app` namespace that provides the `gl_application` base class
/// and related windowing and application infrastructure.
using namespace base_app;

#define MAX_LOADSTRING 100

/// @brief Dialog procedure for the standard "About" dialog box.
///
/// Handles initialization and command notifications for the About dialog that
/// is displayed from the application menu.
///
/// @param hDlg Handle to the dialog box window.
/// @param message Current Windows message.
/// @param wParam Additional message-specific information (WPARAM).
/// @param lParam Additional message-specific information (LPARAM).
/// @return Returns `TRUE` if the message is processed, otherwise `FALSE`.
INT_PTR CALLBACK About(HWND hDlg, UINT message, WPARAM wParam, LPARAM lParam)
{
    UNREFERENCED_PARAMETER(lParam);
    switch (message)
    {
    case WM_INITDIALOG:
        // Perform any one-time initialization for the About dialog here.
        return (INT_PTR)TRUE;

    case WM_COMMAND:
        // Close the dialog when the user presses OK or Cancel.
        if (LOWORD(wParam) == IDOK || LOWORD(wParam) == IDCANCEL)
        {
            EndDialog(hDlg, LOWORD(wParam));
            return (INT_PTR)TRUE;
        }
        break;
    }
    return (INT_PTR)FALSE;
}

/// @brief Main application class for the Mesh Tools viewer.
///
/// `mesh_tools` configures the main window, initializes the mesh view,
/// and forwards input and render callbacks to the `mesh_view` instance.
class mesh_tools : public gl_application {
    /// Title bar text for the main window.
    char m_szTitle[MAX_LOADSTRING];

    /// Mesh view handler responsible for OpenGL rendering and user interaction.
    std::unique_ptr<mesh_view> m_mesh_view;

public:
    /// @brief Constructs the application and creates the mesh view instance.
    mesh_tools() {
        m_mesh_view.reset(new mesh_view());
    }

    /// @brief Pre-creation hook used to configure the main Win32 window.
    ///
    /// This method loads the window title and class name from resources,
    /// sets preferred window dimensions, and customizes the window class
    /// (icons, cursor, background, and menu).
    ///
    /// @param hInstance Handle to the current application instance.
    /// @return Returns 1 on success.
    virtual int precreate_window(HINSTANCE hInstance) {
        char m_szWindowClass[MAX_LOADSTRING];            // the main window class name
        LoadString(hInstance, IDS_APP_TITLE, m_szTitle, MAX_LOADSTRING);
        LoadString(hInstance, IDC_MESH_TOOLS, m_szWindowClass, MAX_LOADSTRING);
        szWindowClass = m_szWindowClass;

        m_window.szTitle = m_szTitle;
        m_window.prefered_width = 1200;
        m_window.prefered_height = 700;

        gl_application::precreate_window(hInstance);

        m_wcex.hIcon = LoadIcon(hInstance, MAKEINTRESOURCE(IDI_MESH_TOOLS));
        m_wcex.hIconSm = LoadIcon(hInstance, MAKEINTRESOURCE(IDI_SMALL));
        m_wcex.hCursor = LoadCursor(NULL, IDC_ARROW);
        m_wcex.hbrBackground = (HBRUSH)(COLOR_WINDOW + 1);
        m_wcex.lpszMenuName = MAKEINTRESOURCE(IDC_MESH_TOOLS);
        m_wcex.lpszClassName = szWindowClass.c_str();

        return 1;
    }

    /// @brief Virtual destructor for cleanup of application resources.
    virtual ~mesh_tools() {
    }

    /// @brief Initializes the application-specific components.
    ///
    /// This currently initializes the `mesh_view` which sets up any
    /// OpenGL state and mesh-related resources needed for rendering.
    ///
    /// @return Returns 1 on successful initialization.
    virtual int init_application() {
        // initialize the mesh view
        m_mesh_view->initialize();

        return 1;
    }

    /// @brief Renders a single frame.
    ///
    /// Delegates rendering to the underlying `mesh_view` instance.
    virtual void render() {
        m_mesh_view->render();
    }

    /// @brief Handles menu and command notifications.
    ///
    /// Processes high-level application commands such as exit and the
    /// About dialog, and forwards all other commands to `mesh_view`.
    ///
    /// @param cmd The command identifier (e.g., menu item ID).
    /// @return Returns 1 after the command is handled.
    virtual int onCommand(int cmd) {
        switch (cmd) {
        case IDM_EXIT:
            PostQuitMessage(0);
            break;
        case IDM_ABOUT:
            DialogBox(hInst, MAKEINTRESOURCE(IDD_ABOUTBOX), m_window.hWnd, About);
            break;
        default:
            m_mesh_view->onCommand(cmd);
            break;
        }
        return 1;
    }

    /// @brief Mouse-move event handler.
    ///
    /// Forwards the event to the `mesh_view` so it can update hover
    /// feedback, camera manipulation, or interaction state.
    ///
    /// @param x Current mouse X position in client coordinates.
    /// @param y Current mouse Y position in client coordinates.
    /// @param extra Additional flags or state encoded by the base framework.
    virtual void onMouseMove(int x, int y, unsigned __int64 extra) {
        m_mesh_view->onMouseMove(x, y, extra);
    }

    /// @brief Left mouse button down event handler.
    ///
    /// Forwards the event to `mesh_view` to begin selections or camera drags.
    ///
    /// @param x Mouse X position in client coordinates.
    /// @param y Mouse Y position in client coordinates.
    /// @param extra Additional flags or state encoded by the base framework.
    virtual void onLMouseDown(int x, int y, unsigned __int64 extra) {
        m_mesh_view->onLMouseDown(x, y, extra);
    }

    /// @brief Left mouse button up event handler.
    ///
    /// Forwards the event to `mesh_view` to finalize selections or drags.
    ///
    /// @param x Mouse X position in client coordinates.
    /// @param y Mouse Y position in client coordinates.
    /// @param extra Additional flags or state encoded by the base framework.
    virtual void onLMouseUp(int x, int y, unsigned __int64 extra) {
        m_mesh_view->onLMouseUp(x, y, extra);
    }

    /// @brief Mouse wheel event handler.
    ///
    /// Typically used by `mesh_view` for zooming or scrolling behavior.
    ///
    /// @param delta Wheel delta (positive or negative).
    /// @param extra_btn Additional button state (`WPARAM` from the message).
    virtual void onMouseWheel(int delta, WPARAM extra_btn) {
        m_mesh_view->onMouseWheel(delta, extra_btn);
    }

    /// @brief Window resize event handler.
    ///
    /// Notifies `mesh_view` of the new client size so that it can update
    /// the OpenGL viewport and projection.
    ///
    /// @param width New client area width in pixels.
    /// @param height New client area height in pixels.
    virtual void resize_window(int width, int height) {
        m_mesh_view->resize(width, height);
    }

    virtual void onRMouseDown(int x, int y, unsigned __int64 extra) {
        m_mesh_view->onRMouseDown(x, y, extra);
    }
    virtual void onRMouseUp(int x, int y, unsigned __int64 extra) {
        m_mesh_view->onRMouseUp(x, y, extra);
    }


};

/// @brief Global application instance.
///
/// The lifetime of this object controls initialization and shutdown
/// of the Mesh Tools application.
mesh_tools my_app;