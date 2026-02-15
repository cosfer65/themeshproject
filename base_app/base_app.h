#ifndef __application_h__
#define __application_h__

#include <Windows.h>
#include <memory>
#include <string>

namespace base_opengl {
    class gl_shader;
    class gl_light;
}

/**
 * @brief Outputs a formatted debug string to the debugger output window.
 * @param format The format string.
 * @param ... Additional arguments for formatting.
 */
extern "C" void __cdecl odprintf(const char* format, ...);

namespace base_app {

#ifdef _DEBUG
    /**
     * @brief Outputs a debug message with a float variable.
     * @param prompt The message prompt.
     * @param var The float variable to output.
     */
    void debug_out(const char* prompt, float var);
#endif

    /**
     * @brief Outputs a formatted debug string to the debugger output window.
     * @param format The format string.
     * @param ... Additional arguments for formatting.
     */
    void __cdecl odprintf(const char* format, ...);

    /**
     * @brief Structure representing an OpenGL window and its properties.
     */
    struct gl_window {
        std::string szTitle;      ///< The window title.
        int prefered_width;       ///< Preferred window width.
        int prefered_height;      ///< Preferred window height.
        int current_width;        ///< Current window width.
        int current_height;       ///< Current window height.
        HWND  hWnd;               ///< Window handle.
        HDC   hDC;                ///< Device context handle.
        HGLRC hRC;                ///< OpenGL rendering context handle.
        bool isMinimized;         ///< Indicates if the window is minimized.
    };

    /**
     * @brief Base class for an OpenGL application.
     * 
     * Provides window creation, event handling, and simulation loop for OpenGL-based applications.
     */
    class gl_application {
        bool bLooping; ///< Indicates if the main loop is running.

        /**
         * @brief Creates an OpenGL window.
         * @param wid Window width.
         * @param hei Window height.
         * @param bitsPerPixel Color depth.
         * @param title Window title.
         * @param hInstance Application instance handle.
         * @param classname Window class name.
         * @param stencilBuffer Stencil buffer bits.
         * @return Window handle.
         */
        HWND create_GL_window(int wid, int hei, int bitsPerPixel, const char* title, HINSTANCE hInstance, const char* classname, int stencilBuffer);

        /**
         * @brief Destroys the OpenGL window and releases resources.
         * @return True if successful, false otherwise.
         */
        bool destroy_GL_window();

    protected:
        HINSTANCE hInst;                        ///< Application instance handle.
        std::string szWindowClass;              ///< Window class name.
        std::unique_ptr<base_opengl::gl_shader> base_shader; ///< Base OpenGL shader.
        bool b_paused = false;                  ///< Indicates if the simulation is paused.
        WNDCLASSEX m_wcex;                      ///< Window class structure.

    public:
        gl_window m_window;     ///< OpenGL window properties.

        /**
         * @brief Constructs a gl_application object.
         */
        gl_application();

        /**
         * @brief Destroys the gl_application object and releases resources.
         */
        virtual ~gl_application();

        /**
         * @brief Initializes the application.
         * @param hInstance Application instance handle.
         * @return 0 on success, non-zero on failure.
         */
        int init(HINSTANCE hInstance);

        /**
         * @brief Runs the main application loop.
         * @return Exit code.
         */
        int run();

        /**
         * @brief Called before window creation for custom setup.
         * @param hInstance Application instance handle.
         * @return 0 on success, non-zero on failure.
         */
        virtual int precreate_window(HINSTANCE hInstance);

        /**
         * @brief Terminates the application and releases resources.
         */
        virtual void terminate();

        /**
         * @brief Handles window resizing.
         * @param width New window width.
         * @param height New window height.
         */
        virtual void resize_window(int width, int height);

        /**
         * @brief Initializes application-specific resources.
         * @return 0 on success, non-zero on failure.
         */
        virtual int init_application();

        /**
         * @brief Advances the simulation by the elapsed time.
         * @param fElapsed Elapsed time in seconds.
         */
        virtual void step_simulation(float fElapsed);

        /**
         * @brief Pauses the simulation.
         * @param fElapsed Elapsed time in seconds.
         */
        virtual void pause_simulation(float fElapsed);

        /**
         * @brief Resumes the simulation.
         * @param fElapsed Elapsed time in seconds.
         */
        virtual void resume_simulation(float fElapsed);

        /**
         * @brief Restarts the simulation.
         */
        virtual void restart_simulation();

        /**
         * @brief Renders the current frame.
         */
        virtual void render();

        /**
         * @brief Exits the application.
         */
        virtual void exit_application();

        /**
         * @brief Handles key down events.
         * @param keycode Virtual key code.
         */
        virtual void onKeyDown(int keycode);

        /**
         * @brief Handles key up events.
         * @param keycode Virtual key code.
         */
        virtual void onKeyUp(int keycode);

        /**
         * @brief Handles window minimized events.
         * @param minimized True if minimized, false otherwise.
         */
        virtual void windowMinimized(bool minimized);

        /**
         * @brief Handles window maximized events.
         */
        virtual void windowMaximized();

        /**
         * @brief Handles command events.
         * @param cmd Command identifier.
         * @return 0 if handled, non-zero otherwise.
         */
        virtual int onCommand(int cmd);

        /**
         * @brief Handles mouse move events.
         *
         * Called whenever the mouse cursor moves within the client area of the window.
         *
         * @param x The X-coordinate of the mouse cursor, in client coordinates.
         * @param y The Y-coordinate of the mouse cursor, in client coordinates.
         * @param extra Additional mouse data or key/button state (implementation-defined, e.g., WPARAM/LPARAM flags).
         */
        virtual void onMouseMove(int x, int y, unsigned __int64 extra) {}

        /**
         * @brief Handles left mouse button down events.
         *
         * Called when the left mouse button is pressed within the client area.
         *
         * @param x The X-coordinate of the mouse cursor at the time of the event, in client coordinates.
         * @param y The Y-coordinate of the mouse cursor at the time of the event, in client coordinates.
         * @param extra Additional mouse data or key/button state (implementation-defined, e.g., WPARAM/LPARAM flags).
         */
        virtual void onLMouseDown(int x, int y, unsigned __int64 extra) {}

        /**
         * @brief Handles left mouse button up events.
         *
         * Called when the left mouse button is released within the client area.
         *
         * @param x The X-coordinate of the mouse cursor at the time of the event, in client coordinates.
         * @param y The Y-coordinate of the mouse cursor at the time of the event, in client coordinates.
         * @param extra Additional mouse data or key/button state (implementation-defined, e.g., WPARAM/LPARAM flags).
         */
        virtual void onLMouseUp(int x, int y, unsigned __int64 extra) {}

        /**
         * @brief Handles left mouse button double-click events.
         *
         * Called when the left mouse button is double-clicked within the client area.
         *
         * @param x The X-coordinate of the mouse cursor at the time of the event, in client coordinates.
         * @param y The Y-coordinate of the mouse cursor at the time of the event, in client coordinates.
         * @param extra Additional mouse data or key/button state (implementation-defined, e.g., WPARAM/LPARAM flags).
         */
        virtual void onLDblClick(int x, int y, unsigned __int64 extra) {}

        /**
         * @brief Handles right mouse button down events.
         *
         * Called when the right mouse button is pressed within the client area.
         *
         * @param x The X-coordinate of the mouse cursor at the time of the event, in client coordinates.
         * @param y The Y-coordinate of the mouse cursor at the time of the event, in client coordinates.
         * @param extra Additional mouse data or key/button state (implementation-defined, e.g., WPARAM/LPARAM flags).
         */
        virtual void onRMouseDown(int x, int y, unsigned __int64 extra) {}

        /**
         * @brief Handles right mouse button up events.
         *
         * Called when the right mouse button is released within the client area.
         *
         * @param x The X-coordinate of the mouse cursor at the time of the event, in client coordinates.
         * @param y The Y-coordinate of the mouse cursor at the time of the event, in client coordinates.
         * @param extra Additional mouse data or key/button state (implementation-defined, e.g., WPARAM/LPARAM flags).
         */
        virtual void onRMouseUp(int x, int y, unsigned __int64 extra) {}

        /**
         * @brief Handles right mouse button double-click events.
         *
         * Called when the right mouse button is double-clicked within the client area.
         *
         * @param x The X-coordinate of the mouse cursor at the time of the event, in client coordinates.
         * @param y The Y-coordinate of the mouse cursor at the time of the event, in client coordinates.
         * @param extra Additional mouse data or key/button state (implementation-defined, e.g., WPARAM/LPARAM flags).
         */
        virtual void onRDblClick(int x, int y, unsigned __int64 extra) {}

        /**
         * @brief Handles middle mouse button down events.
         *
         * Called when the middle mouse button is pressed within the client area.
         *
         * @param x The X-coordinate of the mouse cursor at the time of the event, in client coordinates.
         * @param y The Y-coordinate of the mouse cursor at the time of the event, in client coordinates.
         * @param extra Additional mouse data or key/button state (implementation-defined, e.g., WPARAM/LPARAM flags).
         */
        virtual void onMMouseDown(int x, int y, unsigned __int64 extra) {}

        /**
         * @brief Handles middle mouse button up events.
         *
         * Called when the middle mouse button is released within the client area.
         *
         * @param x The X-coordinate of the mouse cursor at the time of the event, in client coordinates.
         * @param y The Y-coordinate of the mouse cursor at the time of the event, in client coordinates.
         * @param extra Additional mouse data or key/button state (implementation-defined, e.g., WPARAM/LPARAM flags).
         */
        virtual void onMMouseUp(int x, int y, unsigned __int64 extra) {}

        /**
         * @brief Handles middle mouse button double-click events.
         *
         * Called when the middle mouse button is double-clicked within the client area.
         *
         * @param x The X-coordinate of the mouse cursor at the time of the event, in client coordinates.
         * @param y The Y-coordinate of the mouse cursor at the time of the event, in client coordinates.
         * @param extra Additional mouse data or key/button state (implementation-defined, e.g., WPARAM/LPARAM flags).
         */
        virtual void onMDblClick(int x, int y, unsigned __int64 extra) {}

        /**
         * @brief Handles mouse wheel scroll events.
         *
         * Called when the mouse wheel is rotated while the cursor is over the window
         * or the window has input focus.
         *
         * @param delta The scroll amount and direction. Typically corresponds to
         *              the wheel delta from `WM_MOUSEWHEEL` (positive for forward/up,
         *              negative for backward/down), often in multiples of `WHEEL_DELTA`.
         * @param extra_btn Additional mouse or keyboard state information associated
         *                  with the event (implementation-defined, e.g., `WPARAM` flags
         *                  indicating mouse buttons or modifier keys).
         */
        virtual void onMouseWheel(int delta, unsigned __int64 extra_btn) {}
    };
}

#pragma comment( lib, "base_app.lib" )


#endif // __application_h__
