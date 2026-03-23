#pragma once

#include <windows.h>
#include "app_core.h" // force link of library

class cWindow {
    void handle_mouse_message(UINT message, WPARAM wParam, LPARAM lParam);
public:
    HWND      hWnd = nullptr;
    cWindow() {};
    virtual ~cWindow() {}
    // Static window procedure
    static LRESULT CALLBACK StaticWndProc(HWND hWnd, UINT msg, WPARAM wParam, LPARAM lParam);
    // Instance window procedure
    virtual LRESULT WndProc(UINT msg, WPARAM wParam, LPARAM lParam);

    virtual LRESULT OnPaint() { 
        PAINTSTRUCT ps;
        HDC hdc = BeginPaint(hWnd, &ps);
        EndPaint(hWnd, &ps);

        return 0; 
    }
    virtual LRESULT OnSize(int cx, int cy) { return 0; }
    virtual LRESULT OnCreate() { return 0; }
    virtual LRESULT OnDestroy() { hWnd = nullptr; return 0; }
    virtual LRESULT OnMinimize(int wid, int hei) { return 0; }
    virtual LRESULT OnMinimized(int wid, int hei) { return 0; }
    virtual LRESULT OnRestored(int wid, int hei) { return 0; }
    virtual int onCommand(int cmd) { return 0; }

    virtual void render() {}

    virtual void onMouseMove(int x, int y, unsigned __int64 extra) {}
    virtual void onLMouseDown(int x, int y, unsigned __int64 extra) {}
    virtual void onLMouseUp(int x, int y, unsigned __int64 extra) {}
    virtual void onLDblClick(int x, int y, unsigned __int64 extra) {}
    virtual void onRMouseDown(int x, int y, unsigned __int64 extra) {}
    virtual void onRMouseUp(int x, int y, unsigned __int64 extra) {}
    virtual void onRDblClick(int x, int y, unsigned __int64 extra) {}
    virtual void onMMouseDown(int x, int y, unsigned __int64 extra) {}
    virtual void onMMouseUp(int x, int y, unsigned __int64 extra) {}
    virtual void onMDblClick(int x, int y, unsigned __int64 extra) {}
    virtual void onMouseWheel(int delta, unsigned __int64 extra) {}
    virtual void onKeyDown(int key) {}
    virtual void onKeyUp(int key) {}
};

class ViewWindow :public cWindow {
public:
    ViewWindow();
    virtual ~ViewWindow() {}

    virtual LRESULT OnPaint();
};

class glViewWindow : public ViewWindow {
    bool SetupPixelFormat();
public:
    HGLRC hGLRC = nullptr;
    HDC   hDC = nullptr;

    void begin_render();
    void end_render();

    glViewWindow();
    virtual ~glViewWindow() {}

    virtual LRESULT OnCreate();
    virtual LRESULT OnPaint();
    virtual LRESULT OnDestroy();

    virtual void render();
};

class FrameWindow :public cWindow {
public:
    HINSTANCE  hInst = nullptr;
    ViewWindow* pView = nullptr;

    FrameWindow(HINSTANCE hInstance);

    virtual LRESULT OnCreate();
    virtual ViewWindow* get_view() {
        if (pView == nullptr) {
            pView = new ViewWindow();
        }
        return pView;
    }

    virtual LRESULT OnSize(int cx, int cy);
    LRESULT OnDestroy();
};
