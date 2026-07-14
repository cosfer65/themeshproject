#include "application.h"
#include "dynamic_menu.h"
#include "common_dialogs.h"
#include "window.h"
#include "themeshview.h"

static theMeshView* meshView = nullptr;

static void create_application_menu(btm::FrameWindow* pFrame, btm::Menu& menu)
{
    // Create top-level menus
    HMENU fileMenu = menu.create_submenu("File");

    // direct callback in add_item, which is simpler for small apps but less flexible for larger ones where you might want to separate command handling logic
    int openId = menu.add_item(fileMenu, "Open...", [pFrame]() {
        // open file dialog
        const char* fnm = OpenFileDialog("Obj Files\0*.obj\0NURBS Files\0*.nurbs\0All Files\0*.*\0");
        if (fnm) {
            if (meshView) {
                meshView->load_model(fnm);
            }
        }
        });

    menu.add_separator(fileMenu);

    int exitId = menu.add_item(fileMenu, "Exit", [pFrame]() {
        // quit, we can also post a message to the main window to trigger the close event
        PostQuitMessage(0);
        });

    HMENU editMenu = menu.create_submenu("Edit");

    int flipMeshId = menu.add_item(editMenu, "Flip mesh", []() {
        if (meshView) {
            meshView->flip_mesh();
        }
        return 1;
        });

    HMENU calculateMenu = menu.create_submenu("Calculate");
    int calcCurvsId = menu.add_item(calculateMenu, "Curvatures", []() {
        if (meshView) {
            meshView->calculate_curvatures();
        }
        return 1;
        });
    
    HMENU viewMenu = menu.create_submenu("View");
    int viewMeshId = menu.add_item(viewMenu, "Mesh", []() {
        if (meshView) {
            meshView->toggle_show_mesh();
        }
        return 1;
        });
    
    int viewFaceNormsId = menu.add_item(viewMenu, "Face Normals", []() {
        if (meshView) {
            meshView->toggle_show_face_normals();
        }
        return 1;
        });
    int viewVertexNormsId = menu.add_item(viewMenu, "Vertex Normals", []() {
        if (meshView) {
            meshView->toggle_show_vertex_normals();
        }
        return 1;
        });
    int viewPrincipalId = menu.add_item(viewMenu, "Principal Curvatures", []() {
        if (meshView) {
            meshView->toggle_show_principal_curvatures();
        }
        return 1;
        });
    int viewGaussId = menu.add_item(viewMenu, "Gaussian Curvatures", []() {
        if (meshView) {
            meshView->toggle_show_gaussian_curvatures();
        }
        return 1;
        });
    int viewMeanId = menu.add_item(viewMenu, "Mean Curvatures", []() {
        if (meshView) {
            meshView->toggle_show_mean_curvatures();
        }
        return 1;
        });
    int viewK1K2Id = menu.add_item(viewMenu, "K1/K2 Curvatures", []() {
        if (meshView) {
            meshView->toggle_show_k1_k2_preview();
        }
        return 1;
        });
    
    HMENU helpMenu = menu.create_submenu("Help");
    int aboutId = menu.add_item(helpMenu, "About", []() {
        return 1;
        });

    // Attach to window
    menu.attach_to_window(pFrame->hWnd);
}

int APIENTRY WinMain(HINSTANCE hInstance, HINSTANCE, LPTSTR, int nCmdShow)
{
    btm::application the_app;

    btm::init_framework();

    btm::FrameWindow* pFrame = btm::create_main_window(false, 1000, 750, "TheMeshProject");
    btm::Menu menu;
    create_application_menu(pFrame, menu);

    meshView = new theMeshView();

    while (btm::pollEvents())
    {
        meshView->render();
    }
    the_app.terminate();
 
    return 0;
}