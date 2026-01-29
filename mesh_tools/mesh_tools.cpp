#include "resource.h"
#include "mesh_tools.h"
#include "base_app.h"
#include "mesh_view.h"

using namespace base_app;

#define MAX_LOADSTRING 100

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

class mesh_tools :public gl_application {
    char m_szTitle[MAX_LOADSTRING];                     ///< The title bar text

    std::unique_ptr <mesh_view> m_mesh_view;               ///< Mesh view handler   
public:
    mesh_tools() {
        m_mesh_view.reset(new mesh_view());
    }

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

    virtual ~mesh_tools() {
    }

    virtual int init_application() {
        // initialize the mesh view
        m_mesh_view->initialize();

        return 1;
    }

    virtual void render() {
        m_mesh_view->render();
    }

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

    virtual void onMouseMove(int dx, int dy, WPARAM extra_btn) {
        m_mesh_view->mouse_move(dx, dy, extra_btn);
    }

    virtual void onMouseWheel(int delta, WPARAM extra_btn) {
        m_mesh_view->mouse_wheel(delta, extra_btn);
    }

    virtual void resize_window(int width, int height) {
        m_mesh_view->resize(width, height);
    }

};

/// Global application instance
mesh_tools my_app;