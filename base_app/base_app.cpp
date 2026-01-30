#include "gl_graphics.h"

#include "base_app.h"
#include "glew.h"
#include "gl_timer.h"
#include "gl_shaders.h"

using namespace base_opengl;

//-----------------------------------------------------------------------------
// Main application instance pointer
//-----------------------------------------------------------------------------
base_app::gl_application* theApp = NULL;

#ifdef _DEBUG
/**
 * @brief Outputs a debug message with a prompt and float value.
 * @param prompt The message prompt.
 * @param var The float variable to display.
 */
void debug_out(const char* prompt, float var)
{
    char txt[200];
    sprintf(txt, "%s:%f", prompt, var);
    OutputDebugString(txt);
}
#endif

/**
 * @brief Outputs a formatted debug string to the debugger.
 * @param format The format string.
 * @param ... Additional arguments for formatting.
 */
void __cdecl odprintf(const char* format, ...)
{
    char    buf[4096], * p = buf;
    va_list args;
    int     n;

    va_start(args, format);
    n = _vsnprintf(p, sizeof buf - 3, format, args); // buf-3 is room for CR/LF/NUL
    va_end(args);

    if (n < 0)
    {
        n = (int)(sizeof buf - 3);
        p[n] = '\0';
    }
    else if (n > (int)(sizeof buf - 3))
    {
        n = (int)(sizeof buf - 3);
        p[n] = '\0';
    }

    p += n;

    while (p > buf && isspace((unsigned char)p[-1]))
        *--p = '\0';

    *p++ = '\r';
    *p++ = '\n';
    *p = '\0';

    OutputDebugString(buf);
}

/**
 * @brief Entry point for the Windows application.
 * @param hInstance Handle to the current instance.
 * @param hPrevInstance Handle to the previous instance.
 * @param lpCmdLine Command line arguments.
 * @param nCmdShow Show window flag.
 * @return Application exit code.
 */
int WINAPI WinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance, LPTSTR lpCmdLine, int nCmdShow)
{
    UNREFERENCED_PARAMETER(hPrevInstance);
    UNREFERENCED_PARAMETER(lpCmdLine);

    if (theApp == NULL)
        return 0;

    if (theApp->init(hInstance))
    {
        theApp->run();
        theApp->terminate();
    }

    return 0;
}

/**
 * @brief Dispatches raw Windows mouse messages to the active gl_application instance.
 *
 * Extracts the cursor position and auxiliary mouse flags from the Windows
 * message parameters and forwards the event to the corresponding mouse
 * handler on theApp. This centralizes translation from Win32 messages
 * to the engine's input callbacks.
 *
 * @param message The Win32 mouse message identifier (e.g., WM_MOUSEMOVE, WM_LBUTTONDOWN).
 * @param wParam  Message-specific flags (button state, wheel delta, key modifiers).
 *                The lower 7 bits are masked and passed as the extra parameter
 *                to the application-level handlers.
 * @param lParam  Encodes the cursor position in client coordinates:
 *                - low word: X position
 *                - high word: Y position
 *
 * @note The following messages are handled and mapped:
 *       - WM_MOUSEMOVE       -> gl_application::onMouseMove
 *       - WM_LBUTTONDOWN     -> gl_application::onLMouseDown
 *       - WM_LBUTTONUP       -> gl_application::onLMouseUp
 *       - WM_LBUTTONDBLCLK   -> gl_application::onLDblClick
 *       - WM_RBUTTONDOWN     -> gl_application::onRMouseDown
 *       - WM_RBUTTONUP       -> gl_application::onRMouseUp
 *       - WM_RBUTTONDBLCLK   -> gl_application::onRDblClick
 *       - WM_MBUTTONDOWN     -> gl_application::onMMouseDown
 *       - WM_MBUTTONUP       -> gl_application::onMMouseUp
 *       - WM_MBUTTONDBLCLK   -> gl_application::onMDblClick
 *       - WM_MOUSEWHEEL      -> gl_application::onMouseWheel
 *         (wheel delta is extracted via GET_WHEEL_DELTA_WPARAM).
 */
static void handle_mouse_message(UINT message, WPARAM wParam, LPARAM lParam) {
    int x, y;
    x = LOWORD(lParam);
    y = HIWORD(lParam);
    unsigned __int64 extra = wParam & 0x7f;

    switch (message)
    {
    case WM_MOUSEMOVE:
        theApp->onMouseMove(x, y, extra);
        break;
    case WM_LBUTTONDOWN:
        theApp->onLMouseDown(x, y, extra);
        break;
    case WM_LBUTTONUP:
        theApp->onLMouseUp(x, y, extra);
        break;
    case WM_LBUTTONDBLCLK:
        theApp->onLDblClick(x, y, extra);
        break;
    case WM_RBUTTONDOWN:
        theApp->onRMouseDown(x, y, extra);
        break;
    case WM_RBUTTONUP:
        theApp->onRMouseUp(x, y, extra);
        break;
    case WM_RBUTTONDBLCLK:
        theApp->onRDblClick(x, y, extra);
        break;
    case WM_MBUTTONDOWN:
        theApp->onMMouseDown(x, y, extra);
        break;
    case WM_MBUTTONUP:
        theApp->onMMouseUp(x, y, extra);
        break;
    case WM_MBUTTONDBLCLK:
        theApp->onMDblClick(x, y, extra);
        break;

    case WM_MOUSEWHEEL:
        theApp->onMouseWheel(GET_WHEEL_DELTA_WPARAM(wParam), extra);
        break;
    }
}


/**
 * @brief Window procedure for handling Windows messages.
 * @param hWnd Window handle.
 * @param message Message identifier.
 * @param wParam Additional message information.
 * @param lParam Additional message information.
 * @return Result of message processing.
 */
static LRESULT CALLBACK WndProc(HWND hWnd, UINT message, WPARAM wParam, LPARAM lParam)
{
    PAINTSTRUCT ps;
    HDC hdc;
    static int mouse_set = 0;
    static int last_mouse_x;
    static int last_mouse_y;

    switch (message)
    {
    case WM_SIZE:                                                     // window size change
        switch (wParam)                                               // evaluate
        {
        case SIZE_MINIMIZED:                                          // minimized?
            theApp->windowMinimized(true);                            // set isMinimized to true
            return 0;

        case SIZE_MAXIMIZED:                                          // maximized?
            theApp->windowMaximized();                                // set isMinimized to false
            theApp->resize_window(LOWORD(lParam), HIWORD(lParam));    // recalc view - LoWord=Width, HiWord=Height
            return 0;

        case SIZE_RESTORED:                                           // restored?
            theApp->windowMinimized(false);                           // set isMinimized to false
            theApp->resize_window(LOWORD(lParam), HIWORD(lParam));    // recalc view - LoWord=Width, HiWord=Height
            return 0;
        }
        break;
    case WM_KEYDOWN:                                                  // Update Keyboard Buffers For keyDown Pressed
        if ((wParam >= 0) && (wParam <= 255)) {                       // Is Key (wParam) In A Valid Range?
            theApp->onKeyDown((int)wParam);
            return 0;
        }
        break;                                                        // Break

    case WM_KEYUP:                                                    // Update Keyboard Buffers For keyDown Released
        if ((wParam >= 0) && (wParam <= 255)) {                       // Is Key (wParam) In A Valid Range?
            theApp->onKeyUp((int)wParam);                             // Set The Selected Key (wParam) To False
            return 0;                                                 // Return
        }
        break;                                                        // Break

    case WM_MOUSEMOVE:
    case WM_LBUTTONDOWN:
    case WM_LBUTTONUP:
    case WM_LBUTTONDBLCLK:
    case WM_RBUTTONDOWN:
    case WM_RBUTTONUP:
    case WM_RBUTTONDBLCLK:
    case WM_MBUTTONDOWN:
    case WM_MBUTTONUP:
    case WM_MBUTTONDBLCLK:
    case WM_MOUSEWHEEL:
        handle_mouse_message(message, wParam, lParam);
        break;

    case WM_PAINT:
        hdc = BeginPaint(hWnd, &ps);
        EndPaint(hWnd, &ps);
        break;
    case WM_DESTROY:
        PostQuitMessage(0);
        break;
    case WM_COMMAND:
    {
        int wmId = LOWORD(wParam);
        // Parse the menu selections:
        int command_parsed = 0;
        if (theApp)
            command_parsed = theApp->onCommand(wmId);
        if (command_parsed)
            return DefWindowProc(hWnd, message, wParam, lParam);
    }
    break;
    default:
        return DefWindowProc(hWnd, message, wParam, lParam);
    }
    return 0;
}
///////////////////////////////////////////////////////////////////////////////////////
namespace base_app {

    /**
     * @brief Constructs the gl_application object and initializes window parameters.
     */
    gl_application::gl_application()
    {
        theApp = this;

        m_window.szTitle = "Physics simulation";
        m_window.prefered_width = 800;
        m_window.prefered_height = 600;
        m_window.current_width = 0;
        m_window.current_height = 0;

        szWindowClass = "Physics simulation";

        hInst = 0;
        m_window.hWnd = 0;
        m_window.hDC = 0;
        m_window.hRC = 0;
        m_window.isMinimized = false;
        bLooping = false;
    }

    /**
     * @brief Destructor for gl_application.
     */
    gl_application::~gl_application()
    {
    }

    /**
     * @brief Runs the main application loop.
     * Initializes input, OpenGL, and enters the message/render loop.
     * @return Exit code from the message loop.
     */
    int gl_application::run()
    {
        // zero the keystatus buffers
        memset(keyDown, 0, sizeof(keyDown));

        // initialize wrangler library
        glewInit();
        // initialize game
        init_application();
        // start the timer before we enter the main loop
        get_global_timer()->start();
        get_global_timer()->get_elapsed_time();
        // Main message loop:
        bLooping = true;
        MSG msg;
        while (bLooping)
        {
            // check for windows message and precess them
            if (PeekMessage(&msg, NULL, 0, 0, PM_REMOVE) != 0)
            {
                if (msg.message == WM_QUIT)
                    bLooping = false;

                TranslateMessage(&msg);
                DispatchMessage(&msg);
            }
            else  // no messages, just loop for next frame
            {
                if (m_window.isMinimized)        // if window is minimized
                {
                    WaitMessage();        // yeld back to system, do not waste processing power
                }
                else
                {
                    float fElapsed = (float)get_global_timer()->get_elapsed_time();

                    // update the scene based on the elapsed time since last update
                    step_simulation(fElapsed);
                    // render scene
                    render();
                    // Swap Buffers (Double Buffering)
                    SwapBuffers(m_window.hDC);
#ifdef _DEBUG
                    // frame counting mechanism
                    // the higher the frame rate the faster the system
                    // if frame rate drops below 30 we are in deep trouble
                    // we should either optimize the program or buy a new computer
                    static int m_nFrames = 0;                // frame Counter
                    static float tot = 0;                    // time couner
                    tot += fElapsed;                         // increment counters
                    m_nFrames++;
                    if (tot >= 1.f)                          // one second reached
                    {
                        char txt[200];
                        sprintf_s(txt, "%s, fps:%d", m_window.szTitle.c_str(), m_nFrames);
                        SetWindowText(m_window.hWnd, txt);
                        tot = 0;                             // reset counters
                        m_nFrames = 0;
                    }
#endif
                }
            }
        }
        exit_application();
        destroy_GL_window();
        // stop the clock
        get_global_timer()->stop();

        return (int)msg.wParam;
    }

    /**
     * @brief Prepares the window class for creation.
     * @param hInstance Application instance handle.
     * @return 1 on success.
     */
    int gl_application::precreate_window(HINSTANCE hInstance) {
        // create the app window, starting with the class
        ZeroMemory(&m_wcex, sizeof(WNDCLASSEX));
        m_wcex.cbSize = sizeof(WNDCLASSEX);
        m_wcex.style = CS_HREDRAW | CS_VREDRAW | CS_OWNDC | CS_DBLCLKS;
        m_wcex.lpfnWndProc = WndProc;
        m_wcex.cbClsExtra = 0;
        m_wcex.cbWndExtra = 0;
        m_wcex.hInstance = hInstance;
        m_wcex.hIcon = NULL;
        m_wcex.hCursor = LoadCursor(NULL, IDC_ARROW);
        m_wcex.hbrBackground = (HBRUSH)(COLOR_APPWORKSPACE);
        m_wcex.lpszMenuName = NULL;
        m_wcex.lpszClassName = szWindowClass.c_str();
        m_wcex.hIconSm = NULL;

        return 1;
    }

    /**
     * @brief Initializes the application and creates the main window.
     * @param hInstance Application instance handle.
     * @return 1 on success, FALSE on failure.
     */
    int gl_application::init(HINSTANCE hInstance)
    {
        // Perform standard Windows application initialization
        hInst = hInstance; // Store instance handle in our global variable
        precreate_window(hInstance);
        if (RegisterClassEx(&m_wcex) == 0)
        {
            MessageBox(HWND_DESKTOP, "RegisterClassEx Failed!", "Error", MB_OK | MB_ICONEXCLAMATION);
            return false;
        }
        // create the window
        m_window.hWnd = create_GL_window(m_window.prefered_width, m_window.prefered_height, 32, m_window.szTitle.c_str(), hInstance, m_wcex.lpszClassName, 1);
        // if anything goes wrong abort
        if (!m_window.hWnd)
            return FALSE;
        return 1;
    }

    /**
     * @brief Terminates the application.
     */
    void gl_application::terminate()
    {
    }

    /**
     * @brief Resizes the application window.
     * @param width New window width.
     * @param height New window height.
     */
    void gl_application::resize_window(int width, int height)
    {
        m_window.current_width = width;
        m_window.current_height = height;
    }

    /**
     * @brief Creates the OpenGL window and initializes rendering context.
     * @param wid Window width.
     * @param hei Window height.
     * @param bitsPerPixel Color depth.
     * @param title Window title.
     * @param hInstance Application instance handle.
     * @param classname Window class name.
     * @param stencilBuffer Stencil buffer size.
     * @return Window handle on success, 0 on failure.
     */
    HWND gl_application::create_GL_window(int wid, int hei, int bitsPerPixel, const char* title, HINSTANCE hInstance, const char* classname, int stencilBuffer)
    {
        DWORD windowStyle = WS_OVERLAPPEDWINDOW;                // define our window style
        DWORD windowExtendedStyle = WS_EX_APPWINDOW;            // define the window's extended style

        int width = wid;
        int height = hei;

        RECT windowRect = { 0, 0, width, height };              // define our window coordinates

        // adjust window, account for window borders
        AdjustWindowRectEx(&windowRect, windowStyle, 0, windowExtendedStyle);

        // create the opengl window
        HMENU hMenu = NULL;
        if (m_wcex.lpszMenuName)
            hMenu = LoadMenu(hInstance, m_wcex.lpszMenuName);
        m_window.hWnd = CreateWindowEx(windowExtendedStyle,    // extended style
            classname,                                         // class name
            title,                                             // window title
            windowStyle,                                       // window style
            0, 0,                                              // window x,y position
            windowRect.right - windowRect.left,                // window width
            windowRect.bottom - windowRect.top,                // window height
            HWND_DESKTOP,                                      // desktop is window's parent
            hMenu,                                             // no menu
            hInstance,                                         // pass the window instance
            0);

        if (m_window.hWnd == 0)                                // was window creation a success?
        {
            return m_window.hWnd;                              // if not return false
        }

        m_window.hDC = GetDC(m_window.hWnd);                   // grab a device context for this window
        if (m_window.hDC == 0)                                 // did we get a device context?
        {
            // Failed
            DestroyWindow(m_window.hWnd);                      // destroy the window
            m_window.hWnd = 0;                                 // zero the window handle
            return m_window.hWnd;                              // return false
        }

        PIXELFORMATDESCRIPTOR pfd =                            // pfd tells windows how we want things to be
        {
            sizeof(PIXELFORMATDESCRIPTOR),                     // size of this pixel format descriptor
            1,                                                 // version number
            PFD_DRAW_TO_WINDOW |                               // format must support window
            PFD_SUPPORT_OPENGL |                               // format must support opengl
            PFD_DOUBLEBUFFER,                                  // must support double buffering
            PFD_TYPE_RGBA,                                     // request an rgba format
            (BYTE)bitsPerPixel,                                // select our color depth
            0, 0, 0, 0, 0, 0,                                  // color bits ignored
            0,                                                 // no alpha buffer
            0,                                                 // shift bit ignored
            0,                                                 // no accumulation buffer
            0, 0, 0, 0,                                        // accumulation bits ignored
            16,                                                // 16bit z-buffer (depth buffer)
            (BYTE)stencilBuffer,                               // stencil buffer
            0,                                                 // no auxiliary buffer
            PFD_MAIN_PLANE,                                    // main drawing layer
            0,                                                 // reserved
            0, 0, 0                                            // layer masks ignored
        };

        GLuint PixelFormat = ChoosePixelFormat(m_window.hDC, &pfd); // find a compatible pixel format
        if (PixelFormat == 0)                                       // did we find a compatible format?
        {
            // Failed
            ReleaseDC(m_window.hWnd, m_window.hDC);                 // release our device context
            m_window.hDC = 0;                                       // zero the device context
            DestroyWindow(m_window.hWnd);                           // destroy the window
            m_window.hWnd = 0;                                      // zero the window handle
            return m_window.hWnd;                                   // return false
        }

        if (SetPixelFormat(m_window.hDC, PixelFormat, &pfd) == false)// try to set the pixel format
        {
            // Failed
            ReleaseDC(m_window.hWnd, m_window.hDC);                 // release our device context
            m_window.hDC = 0;                                       // zero the device context
            DestroyWindow(m_window.hWnd);                           // destroy the window
            m_window.hWnd = 0;                                      // zero the window handle
            return m_window.hWnd;                                   // return false
        }

        m_window.hRC = wglCreateContext(m_window.hDC);              // try to get a rendering context
        if (m_window.hRC == 0)                                      // did we get a rendering context?
        {
            // Failed
            ReleaseDC(m_window.hWnd, m_window.hDC);                 // release our device context
            m_window.hDC = 0;                                       // zero the device context
            DestroyWindow(m_window.hWnd);                           // destroy the window
            m_window.hWnd = 0;                                      // zero the window handle
            return m_window.hWnd;                                   // return false
        }

        // make the rendering context our current rendering context
        if (wglMakeCurrent(m_window.hDC, m_window.hRC) == false)    // failed
        {
            wglDeleteContext(m_window.hRC);                         // delete the rendering context
            m_window.hRC = 0;                                       // zero the rendering context
            ReleaseDC(m_window.hWnd, m_window.hDC);                 // release our device context
            m_window.hDC = 0;                                       // zero the device context
            DestroyWindow(m_window.hWnd);                           // destroy the window
            m_window.hWnd = 0;                                      // zero the window handle
            return m_window.hWnd;                                   // return false
        }
        //glfwSwapInterval(0);

        ShowWindow(m_window.hWnd, SW_NORMAL);                       // make the window visible
        m_window.isMinimized = false;                               // set isMinimized to false

        resize_window(width, height);
        UpdateWindow(m_window.hWnd);

        return m_window.hWnd;                                       // window creating was a success
    }

    /**
     * @brief Destroys the OpenGL window and releases resources.
     * @return true when resources are released.
     */
    bool gl_application::destroy_GL_window()
    {
        if (m_window.hWnd != 0)                                     // does the window have a handle?
        {
            if (m_window.hDC != 0)                                  // does the window have a device context?
            {
                wglMakeCurrent(m_window.hDC, 0);                    // set the current active rendering context to zero
                if (m_window.hRC != 0)                              // does the window have a rendering context?
                {
                    wglDeleteContext(m_window.hRC);                 // release the rendering context
                    m_window.hRC = 0;                               // zero the rendering context
                }
                ReleaseDC(m_window.hWnd, m_window.hDC);             // release the device context
                m_window.hDC = 0;                                   // zero the device context
            }
            DestroyWindow(m_window.hWnd);                           // destroy the window
            m_window.hWnd = 0;                                      // zero the window handle
        }

        return true;                                                // return true
    }

    /**
     * @brief Initializes the application, compiles shaders, and sets up rendering.
     * @return 1 on success.
     */
    int gl_application::init_application()
    {
        // geometry and color data are stored in the shader
        static const GLchar* vertex_source =
        {
        "#version 330 core"
        "out fvec4 vs_color;"
        "void main(void)"
        "{"
            "const fvec4 colors[3] = fvec4[3](fvec4(1, 0, 0, 1.0),"
            "fvec4(0, 1, 0, 1.0),"
            "fvec4(0, 0, 1, 1.0));"
            // Declare a hard-coded array of positions
            "const fvec4 vertices[3] = fvec4[3](fvec4(-0.5, -0.5, 0, 1.0),"
            "fvec4(0, 0.5, 0, 1.0),"
            "fvec4(.5, -0.5, 0, 1.0));"
            // Index into our array using gl_VertexID
            "gl_Position = vertices[gl_VertexID];"
            "vs_color = colors[gl_VertexID];"
        "}"
        };
        // Source code for fragment shader
        static const GLchar* fragment_source =
        {
        "#version 330 core"
        "in fvec4 vs_color;"
        "out fvec4 color;"
        "void main(void)"
        "{"
        " color = vs_color;"
        "}"
        };
        base_shader.reset(new gl_shader);
        base_shader->compile(vertex_source, fragment_source);
        return 1;
    }

    /**
     * @brief Advances the simulation by the elapsed time.
     * @param fElapsed Time elapsed since last update.
     */
    void gl_application::step_simulation(float fElapsed)
    {
    }

    /**
     * @brief Pauses the simulation.
     * @param fElapsed Time elapsed since last update.
     */
    void gl_application::pause_simulation(float fElapsed) {
        b_paused = true;
    }

    /**
     * @brief Resumes the simulation.
     * @param fElapsed Time elapsed since last update.
     */
    void gl_application::resume_simulation(float fElapsed) {
        b_paused = false;
    }

    /**
     * @brief Restarts the simulation.
     */
    void gl_application::restart_simulation() {
        b_paused = false;
    }

    /**
     * @brief Renders the current frame.
     */
    void gl_application::render()
    {
        glViewport(0, 0, m_window.current_width, m_window.current_height);
        // clear the screen
        glClearColor(0.2f, 0.2f, 0.2f, 1);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        // use shader
        base_shader->use();
        // draw the triangle (vertext coordinates are in the shader)
        glDrawArrays(GL_TRIANGLES, 0, 3);
        base_shader->end();
    }

    /**
     * @brief Cleans up resources and exits the application.
     */
    void gl_application::exit_application()
    {
    }

    /**
     * @brief Handles key down events.
     * @param keycode The key code pressed.
     */
    void gl_application::onKeyDown(int keycode) {
        keyDown[keycode] = true;							 // Set The Selected Key (wParam) To True
    }

    /**
     * @brief Handles key up events.
     * @param keycode The key code released.
     */
    void gl_application::onKeyUp(int keycode) {
        keyDown[keycode] = false;						     // Set The Selected Key (wParam) To False
    }

    /**
     * @brief Sets the window minimized state.
     * @param m True if minimized, false otherwise.
     */
    void gl_application::windowMinimized(bool m) {
        m_window.isMinimized = m;
    }

    /**
     * @brief Sets the window maximized state.
     */
    void gl_application::windowMaximized() {
        m_window.isMinimized = false;					    // set isMinimized to false
    }

    /**
     * @brief Handles command events from the window.
     * @param cmd The command identifier.
     * @return 0 if not handled.
     */
    int gl_application::onCommand(int cmd) {
        return 0;
    }
}