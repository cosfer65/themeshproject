#pragma once

#include "window.h"

class meshFrameWindow : public btm::FrameWindow {
public:

    meshFrameWindow(HINSTANCE hInstance) : btm::FrameWindow(hInstance) {}

    virtual btm::glView* create_view();

    virtual LRESULT OnMinimize(int wid, int hei);
    virtual LRESULT OnMaximize(int wid, int hei);
    virtual LRESULT OnRestored(int wid, int hei);

};
