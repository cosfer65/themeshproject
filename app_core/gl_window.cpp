#define _CRT_SECURE_NO_WARNINGS
#include <windows.h>
#include "gl_graphics.h"
#include "gl_window.h"
#include "core_application.h"

void cWindow::handle_mouse_message(UINT message, WPARAM wParam, LPARAM lParam) {
    int x, y;
    x = LOWORD(lParam);
    y = HIWORD(lParam);
    unsigned __int64 extra = wParam & 0x7f;

    switch (message)
    {
    case WM_MOUSEMOVE:
        onMouseMove(x, y, extra);
        break;
    case WM_LBUTTONDOWN:
        onLMouseDown(x, y, extra);
        break;
    case WM_LBUTTONUP:
        onLMouseUp(x, y, extra);
        break;
    case WM_LBUTTONDBLCLK:
        onLDblClick(x, y, extra);
        break;
    case WM_RBUTTONDOWN:
        onRMouseDown(x, y, extra);
        break;
    case WM_RBUTTONUP:
        onRMouseUp(x, y, extra);
        break;
    case WM_RBUTTONDBLCLK:
        onRDblClick(x, y, extra);
        break;
    case WM_MBUTTONDOWN:
        onMMouseDown(x, y, extra);
        break;
    case WM_MBUTTONUP:
        onMMouseUp(x, y, extra);
        break;
    case WM_MBUTTONDBLCLK:
        onMDblClick(x, y, extra);
        break;

    case WM_MOUSEWHEEL:
        onMouseWheel(GET_WHEEL_DELTA_WPARAM(wParam), extra);
        break;
    }
}

LRESULT cWindow::WndProc(UINT msg, WPARAM wParam, LPARAM lParam) {
    if (msg >= WM_MOUSEFIRST && msg <= WM_MOUSELAST) {
        handle_mouse_message(msg, wParam, lParam);
        return 0;
    }
    switch (msg) {
    case WM_SIZE:                                                     // window size change
        switch (wParam)                                               // evaluate
        {
        case SIZE_MINIMIZED:                                          // minimized?
            OnMinimize(LOWORD(lParam), HIWORD(lParam));
            return 0;

        case SIZE_MAXIMIZED:                                          // maximized?
            OnMinimized(LOWORD(lParam), HIWORD(lParam));
            return 0;

        case SIZE_RESTORED:                                           // restored?
            OnRestored(LOWORD(lParam), HIWORD(lParam));
            return 0;
        default:
            OnSize(LOWORD(lParam), HIWORD(lParam));
            return 0;
        }
        break;

    case WM_KEYDOWN:
        if ((wParam >= 0) && (wParam <= 255)) {                       // Is Key (wParam) In A Valid Range?
            onKeyDown((int)wParam);
            return 0;
        }
        break;

    case WM_KEYUP:
        if ((wParam >= 0) && (wParam <= 255)) {                       // Is Key (wParam) In A Valid Range?
            onKeyUp((int)wParam);                             // Set The Selected Key (wParam) To False
            return 0;                                                 // Return
        }
        break;

    case WM_CREATE:
        return OnCreate();
        break;
    case WM_PAINT:
        return OnPaint();
        break;
    case WM_COMMAND:
        return GetApp()->onCommand(LOWORD(wParam));
        break;
    case WM_DESTROY:
        return OnDestroy();
        break;
    }
    return DefWindowProc(hWnd, msg, wParam, lParam);
}

LRESULT CALLBACK cWindow::StaticWndProc(HWND hWnd, UINT msg, WPARAM wParam, LPARAM lParam) {
    cWindow* pThis = nullptr;

    if (msg == WM_NCCREATE) {
        CREATESTRUCT* cs = reinterpret_cast<CREATESTRUCT*>(lParam);
        pThis = reinterpret_cast<cWindow*>(cs->lpCreateParams);
        SetWindowLongPtr(hWnd, GWLP_USERDATA, (LONG_PTR)pThis);
        pThis->hWnd = hWnd;
    }
    else {
        pThis = reinterpret_cast<cWindow*>(GetWindowLongPtr(hWnd, GWLP_USERDATA));
    }

    if (pThis && pThis->hWnd==hWnd) {
        return pThis->WndProc(msg, wParam, lParam);
    }
    return DefWindowProc(hWnd, msg, wParam, lParam);
}

////////////////////////////////////////////////////////////////
ViewWindow::ViewWindow() {}

LRESULT ViewWindow::OnPaint() {
    PAINTSTRUCT ps;
    HDC hdc = BeginPaint(hWnd, &ps);

    RECT rc;
    GetClientRect(hWnd, &rc);

    char buf[128];
    strcpy(buf, "SDI View");

    DrawText(hdc, buf, -1, &rc, DT_SINGLELINE | DT_CENTER | DT_VCENTER);

    EndPaint(hWnd, &ps);
    return 0;
}

////////////////////////////////////////////////////////////////
glViewWindow::glViewWindow() : ViewWindow() {}

bool glViewWindow::SetupPixelFormat() {
    PIXELFORMATDESCRIPTOR pfd = {
        sizeof(PIXELFORMATDESCRIPTOR),
        1,
        PFD_DRAW_TO_WINDOW | PFD_SUPPORT_OPENGL | PFD_DOUBLEBUFFER,
        PFD_TYPE_RGBA,
        32,
        0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0,
        32,   // depth buffer
        32,    // stencil buffer
        0,
        PFD_MAIN_PLANE,
        0, 0, 0, 0
    };

    int pf = ChoosePixelFormat(hDC, &pfd);
    if (!pf) return false;

    if (!SetPixelFormat(hDC, pf, &pfd))
        return false;

    return true;
}

LRESULT glViewWindow::OnCreate() {
    hDC = GetDC(hWnd);
    if (!SetupPixelFormat()) return -1;

    hGLRC = wglCreateContext(hDC);
    wglMakeCurrent(hDC, hGLRC);

    return 0;
}
LRESULT glViewWindow::OnPaint() {
    PAINTSTRUCT ps;
    BeginPaint(hWnd, &ps);
    EndPaint(hWnd, &ps);
    return 0;
}
LRESULT glViewWindow::OnDestroy() {
    if (hGLRC) {
        wglMakeCurrent(nullptr, nullptr);
        wglDeleteContext(hGLRC);
        hGLRC = nullptr;
    }
    if (hDC) {
        ReleaseDC(hWnd, hDC);
        hDC = nullptr;
    }
    return 0;
}

void glViewWindow::begin_render() {
    wglMakeCurrent(hDC, hGLRC);
}
void glViewWindow::end_render() {
    SwapBuffers(hDC);
}

void glViewWindow::render() {
    if (!hWnd || !hDC || !hGLRC) return;

    begin_render();

    RECT rc;
    GetClientRect(hWnd, &rc);
    int width = rc.right - rc.left;
    int height = rc.bottom - rc.top;
    if (height <= 0) height = 1;

    glViewport(0, 0, width, height);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // Example: simple triangle
    glBegin(GL_TRIANGLES);
    glColor3f(1, 0, 0); glVertex2f(-0.5f, -0.5f);
    glColor3f(0, 1, 0); glVertex2f(0.5f, -0.5f);
    glColor3f(0, 0, 1); glVertex2f(0.0f, 0.5f);
    glEnd();

    end_render();
}

FrameWindow::FrameWindow(HINSTANCE hInstance) : hInst(hInstance) {}

LRESULT FrameWindow::OnDestroy() {
    if (pView) {
        delete pView;
        pView = nullptr;
    }
    PostQuitMessage(0);
    return 0;
}

LRESULT FrameWindow::OnCreate() {
    // Create the view object
    pView = get_view();

    HWND hView = CreateWindowEx(
        0, "ViewWindowClass", nullptr, WS_CHILD | WS_VISIBLE,
        0, 0, 0, 0,
        hWnd, nullptr, hInst,
        pView // lpCreateParams -> 'this' for ViewWindow
    );

    if (!hView) {
        MessageBox(hWnd, "Failed to create view window", "Error", MB_ICONERROR);
        return -1;
    }
    return 0;
}

LRESULT FrameWindow::OnSize(int cx, int cy) {
    if (pView && pView->hWnd) {
        MoveWindow(pView->hWnd, 0, 0, cx, cy, TRUE);
        pView->OnSize(cx, cy);
    }
    return 0;
}