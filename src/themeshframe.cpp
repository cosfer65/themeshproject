#include "themeshframe.h"
#include "themeshview.h"

ViewWindow* meshFrameWindow::get_view() {
    if (pView == nullptr) {
        pView = new meshViewWindow();
    }
    return pView;
}

LRESULT meshFrameWindow::OnMinimize(int wid, int hei) {
    OnSize(wid, hei);
    return 0;
}

LRESULT meshFrameWindow::OnMinimized(int wid, int hei) {
    OnSize(wid, hei);
    return 0;
}

LRESULT meshFrameWindow::OnRestored(int wid, int hei) {
    OnSize(wid, hei);
    return 0;
}