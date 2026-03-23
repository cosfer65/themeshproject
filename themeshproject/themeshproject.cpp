#include "framework.h"
#include "themeshproject.h"
#include "gl_window.h"
#include "core_application.h"

#include "themeshframe.h"
#include "themeshview.h"
#include "themeshmodel.h"

#include "app_core.h"
#include "common_dialogs.h"

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

class meshApplication : public winApplication {
    cModel* m_model = nullptr;
public:
    meshApplication() {}
    ~meshApplication() {}

    virtual void precreate_window(HINSTANCE hInstance, WNDCLASSEX* m_wcex) {
        m_wcex->hIcon = LoadIcon(hInstance, MAKEINTRESOURCE(IDI_MDIAPP));
        m_wcex->hIconSm = LoadIcon(hInstance, MAKEINTRESOURCE(IDI_SMALL));
        m_wcex->lpszMenuName = MAKEINTRESOURCE(IDC_MDIAPP);
    }

    FrameWindow* getMainWindow(HINSTANCE hInstance) {
        if (pFrame == nullptr) {
            pFrame = new meshFrameWindow(hInstance);
        }
        return pFrame;
    }

    virtual cWindow* get_active_view() {
        return getViewWindow();
    }

    meshViewWindow* getViewWindow() {
        if (pFrame) {
            return dynamic_cast<meshViewWindow*>(pFrame->get_view());
        }
        return nullptr;
    }

    bool load_model(const std::string& fnm) {
        if (m_model) {
            delete m_model;
            m_model = nullptr;
        }
        m_model = load_mesh_model(fnm);
        if (m_model) {
            meshViewWindow* pView = getViewWindow();
            if (pView) {
                pView->set_model(m_model);
                pView->reset_view();
            }
        }
        return m_model != nullptr;
    }

    bool generate_sphere() {
        if (m_model) {
            delete m_model;
            m_model = nullptr;
        }
        m_model = generate_sphere_model();
        if (m_model) {
            meshViewWindow* pView = getViewWindow();
            if (pView) {
                pView->set_model(m_model);
                pView->reset_view();
            }
        }
        return m_model != nullptr;
    }

    virtual int onCommand(int cmd) {
        switch (cmd) {
        case ID_FILE_LOADMODEL:
        {
            // generate_sphere();
            const char* cfname = OpenFileDialog("All Files\0*.*\0Obj Files\0*.obj\0NURBS Files\0*.nurbs\0");
            if (cfname) {
                load_model(cfname);
            }
        }
        break;
        case IDM_EXIT:
            PostQuitMessage(0);
            break;
        case IDM_ABOUT:
            DialogBox(get_hInstance(), MAKEINTRESOURCE(IDD_ABOUTBOX), pFrame->hWnd, About);
            break;
        default:
            if (pFrame) {
                pFrame->get_view()->onCommand(cmd);
            }
        }
        return 1;
    }

    virtual void step_simulation(float fElapsed) {
        meshViewWindow* pView = getViewWindow();
        if (pView) {
            pView->step_simulation(fElapsed);
        }
    }
    virtual void render() {
        meshViewWindow* pView = getViewWindow();
        if (pView) {
            pView->render();
        }
    }

    virtual void init_application() {
        meshViewWindow* pView = getViewWindow();
        pView->init();
    }
};

meshApplication gApp;