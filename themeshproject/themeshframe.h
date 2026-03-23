#pragma once

#include "gl_window.h"

class meshFrameWindow : public FrameWindow {
public:

    meshFrameWindow(HINSTANCE hInstance) : FrameWindow(hInstance) {}

    virtual ViewWindow* get_view();
    virtual LRESULT OnMinimize(int wid, int hei);
    virtual LRESULT OnMinimized(int wid, int hei);
    virtual LRESULT OnRestored(int wid, int hei);

};
