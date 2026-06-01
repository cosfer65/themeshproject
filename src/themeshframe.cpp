#include "themeshframe.h"
#include "themeshview.h"

btm::glView* meshFrameWindow::create_view() {
    return new meshViewWindow();
}

LRESULT meshFrameWindow::OnMinimize(int wid, int hei) {
    OnSize(wid, hei);
    return 0;
}

LRESULT meshFrameWindow::OnMaximize(int wid, int hei) {
    OnSize(wid, hei);
    return 0;
}

LRESULT meshFrameWindow::OnRestored(int wid, int hei) {
    OnSize(wid, hei);
    return 0;
}