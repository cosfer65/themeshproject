#include <string>
#include "gl_window.h"
#include "gl_graphics.h"
#include "gl_timer.h"

#include "core_application.h"

using namespace base_opengl;

class winApplication;
static winApplication* theApp = nullptr;
static HINSTANCE ghInstance;
HINSTANCE get_hInstance() {
    return ghInstance;
}

winApplication* GetApp() {
    return theApp;
}

winApplication::winApplication() {
    if (theApp != nullptr) {
        MessageBox(nullptr, "Application instance already exists!", "Error", MB_ICONERROR);
        exit(1);
    }
    theApp = this;
}
winApplication::~winApplication() {
    if (pFrame != nullptr) {
        delete pFrame;
        pFrame = nullptr;
    }
}

FrameWindow* winApplication::getMainWindow(HINSTANCE hInstance) {
    if (pFrame == nullptr) {
        pFrame = new FrameWindow(hInstance);
    }
    return pFrame;
}

static bool RegisterWindowClass(HINSTANCE hInst, LPCSTR className, UINT style, WNDPROC wndProc) {
    WNDCLASS wc = { 0 };
    wc.style = style;
    wc.lpfnWndProc = wndProc;
    wc.hInstance = hInst;
    wc.hCursor = LoadCursor(nullptr, IDC_ARROW);
    wc.hIcon = LoadIcon(nullptr, IDI_APPLICATION);
    wc.hbrBackground = (HBRUSH)(COLOR_BTNFACE + 1);
    wc.lpszClassName = className;
    return RegisterClass(&wc) != 0;
}

int APIENTRY WinMain(HINSTANCE hInstance, HINSTANCE, LPTSTR, int nCmdShow) {
    char szTitle[] = "TheMeshProject";
    char szWindowClass[] = "FrameWindowClass";
    char szViewClass[] = "ViewWindowClass";
    ghInstance = hInstance;

    WNDCLASSEX m_wcex;
    ZeroMemory(&m_wcex, sizeof(WNDCLASSEX));
    m_wcex.cbSize = sizeof(WNDCLASSEX);
    m_wcex.style = CS_HREDRAW | CS_VREDRAW;
    m_wcex.lpfnWndProc = cWindow::StaticWndProc;
    m_wcex.hInstance = hInstance;
    m_wcex.hCursor = LoadCursor(NULL, IDC_ARROW);
    m_wcex.hbrBackground = (HBRUSH)(COLOR_APPWORKSPACE);
    m_wcex.lpszClassName = szWindowClass;

    theApp->precreate_window(hInstance, &m_wcex);

    if (!RegisterClassEx(&m_wcex) ||
        !RegisterWindowClass(hInstance, szViewClass, CS_HREDRAW | CS_VREDRAW | CS_OWNDC | CS_DBLCLKS, cWindow::StaticWndProc)) {
        MessageBox(nullptr, "Failed to register window classes", "Error", MB_ICONERROR);
        return 0;
    }

    FrameWindow* pFrame = theApp->getMainWindow(hInstance);

    DWORD windowStyle = WS_OVERLAPPEDWINDOW;                // define our window style
    DWORD windowExtendedStyle = WS_EX_APPWINDOW;            // define the window's extended style

    RECT windowRect = { 0, 0, 800, 600 };              // define our window coordinates

    // adjust window, account for window borders
    AdjustWindowRectEx(&windowRect, windowStyle, 0, windowExtendedStyle);

    HWND hFrame = CreateWindowEx(
        0,
        szWindowClass,
        szTitle,
        WS_OVERLAPPEDWINDOW,
        CW_USEDEFAULT, CW_USEDEFAULT, windowRect.right - windowRect.left, windowRect.bottom - windowRect.top,
        nullptr,
        nullptr,
        hInstance,
        pFrame // lpCreateParams -> 'this' for FrameWindow
    );

    if (!hFrame) {
        MessageBox(nullptr, "Failed to create frame window", "Error", MB_ICONERROR);
        return 0;
    }

    ShowWindow(hFrame, nCmdShow);
    UpdateWindow(hFrame);

    // initialize wrangler library
    glewInit();
    // initialize application-specific resources
    theApp->init_application();
    // start the timer before we enter the main loop
    get_global_timer()->start();
    get_global_timer()->get_elapsed_time();

    MSG msg;
    bool running = true;

    while (running) {
        if (PeekMessage(&msg, NULL, 0, 0, PM_REMOVE) != 0)
        {
            if (msg.message == WM_QUIT)
                running = false;

            TranslateMessage(&msg);
            DispatchMessage(&msg);
        }
        else  // no messages, just loop for next frame
        {
            // Idle time → render a frame
            float fElapsed = (float)get_global_timer()->get_elapsed_time();
            theApp->step_simulation(fElapsed);
            theApp->render();
        }
    }
    theApp->terminate();
    get_global_timer()->stop();

    return (int)msg.wParam;
}

