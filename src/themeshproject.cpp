#include <windows.h>
#include "themeshproject.h"
#include "timer.h"
#include "window.h"
#include "application.h"

#include "themeshframe.h"
#include "themeshview.h"
#include "themeshmodel.h"

#include "utils.h"
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

class meshApplication : public application {
    cModel* m_model = nullptr;
public:
    meshApplication() {
        // register menu handlers
        ON_COMMAND(this, ID_FILE_LOADMODEL, meshApplication::onFileOpen);
        ON_COMMAND(this, IDM_EXIT, meshApplication::onFileExit);
        ON_COMMAND(this, IDM_ABOUT, meshApplication::onHelpAbout);
    }
    virtual ~meshApplication() {}

    // menu handlers
    int onFileOpen(int cmd) {
        const char* cfname = OpenFileDialog("All Files\0*.*\0Obj Files\0*.obj\0NURBS Files\0*.nurbs\0");
        if (cfname) {
            load_model(cfname);
        }
        return 0;
    }

    int onFileExit(int cmd) {
        PostQuitMessage(0);
        return 0;
    }

    int onHelpAbout(int cmd) {
        FrameWindow* pFrame = getMainWindow(GetModuleHandle(nullptr));
        if (pFrame) {
            DialogBox(get_hInstance(), MAKEINTRESOURCE(IDD_ABOUTBOX), pFrame->hWnd, About);
        }
        return 0;
    }


    virtual void precreate_window(HINSTANCE hInstance, WNDCLASSEX* m_wcex) {
        m_wcex->hIcon = LoadIcon(hInstance, MAKEINTRESOURCE(IDI_MDIAPP));
        m_wcex->hIconSm = LoadIcon(hInstance, MAKEINTRESOURCE(IDI_SMALL));
        m_wcex->lpszMenuName = MAKEINTRESOURCE(IDC_MDIAPP);
    }

    FrameWindow* createMainWindow(HINSTANCE hInstance) {
        return new meshFrameWindow(hInstance);
    }
    virtual cWindow* get_active_view() {
        return getViewWindow();
    }
    //
    meshViewWindow* getViewWindow() {
        FrameWindow* pFrame = getMainWindow(GetModuleHandle(nullptr));
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
        if (pView) {
            pView->init();
        }
    }
};

// Windows application entry point
int APIENTRY WinMain(HINSTANCE hInstance, HINSTANCE, LPTSTR, int nCmdShow)
{
    meshApplication the_app;

    btm::init_framework();

    btm::FrameWindow* pFrame = btm::create_main_window(true, 1000, 750, "TheMeshProject");

    btm::start_timer();
    btm::get_elapsed_time();
    
    while (btm::pollEvents())
    {
        // Idle time → render a frame
        // calculate elapsed time since last frame
        float fElapsed = (float)btm::get_elapsed_time();
        // step the simulation and render
        the_app.step_simulation(fElapsed);
        the_app.render();
    }
    the_app.terminate();
    btm::stop_timer();

    return 0;
}
