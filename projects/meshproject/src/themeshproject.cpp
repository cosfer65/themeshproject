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
        const char* fnm = OpenFileDialog("Obj Files\0*.obj\0All Files\0*.*\0");
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
    // available only in DEBUG builds, for testing purposes
#ifdef _DEBUG
    int view_test_fun = menu.add_item(editMenu, "Test Function", []() {
        if (meshView) {
            meshView->test_function();
        }
        return 1;
        });
#endif


    HMENU calculateMenu = menu.create_submenu("Calculate");
    int calcCurvsId = menu.add_item(calculateMenu, "Curvatures", []() {
        if (meshView) {
            meshView->calculate_curvatures();
        }
        return 1;
        });
    int calcFeatLinesId = menu.add_item(calculateMenu, "Feature Lines", []() {
        if (meshView) {
            meshView->calculate_feature_lines();
        }
        return 1;
        });

    HMENU viewMenu = menu.create_submenu("View");
    int viewMeshId = menu.add_item(viewMenu, "Shaded Mesh", []() {
        if (meshView) {
            meshView->toggle_show_mesh();
        }
        return 1;
        });
    
    int viewWireframeId = menu.add_item(viewMenu, "Wireframe", []() {
        if (meshView) {
            meshView->toggle_show_wireframe();
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
    int viewPrincipalK1Id = menu.add_item(viewMenu, "Principal Curvatures Kmax", []() {
        if (meshView) {
            meshView->toggle_show_principal_k1();
        }
        return 1;
        });
    int viewPrincipalK2Id = menu.add_item(viewMenu, "Principal Curvatures Kmin", []() {
        if (meshView) {
            meshView->toggle_show_principal_k2();
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
    int viewK1K2Id = menu.add_item(viewMenu, "Kmax/Kmin Curvatures", []() {
        if (meshView) {
            meshView->toggle_show_k1_k2_preview();
        }
        return 1;
        });

    int viewRidges = menu.add_item(viewMenu, "Ridges", []() {
        if (meshView) {
            meshView->toggle_show_ridges();
        }
        return 1;
        });
    int viewValleys = menu.add_item(viewMenu, "Valleys", []() {
        if (meshView) {
            meshView->toggle_show_valleys();
        }
        return 1;
        });
    int viewCreases = menu.add_item(viewMenu, "Creases", []() {
        if (meshView) {
            meshView->toggle_show_creases();
        }
        return 1;
        });
    int viewBoundaries = menu.add_item(viewMenu, "Boundaries", []() {
        if (meshView) {
            meshView->toggle_show_boundaries();
        }
        return 1;
        });


    int viewResetId = menu.add_item(viewMenu, "Reset View", []() {
        if (meshView) {
            meshView->reset_view_state();
        }
        return 1;
        });


    // dialog boxes are not yet implemented, so we will not add a help menu for now
    // HMENU helpMenu = menu.create_submenu("Help");
    // int aboutId = menu.add_item(helpMenu, "About", []() {
    //     return 1;
    //     });

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