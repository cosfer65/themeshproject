#include "gl_graphics.h"
#include "gl_font.h"

// works ONLY in Windows
static PFNGLWINDOWPOS2IPROC glWindowPos2i = NULL;

namespace base_opengl {

    /**
     * @brief Font map type for managing font name to font object mapping.
     */
    cg_fontMapper gl_font_manager::m_fontMap;

    /**
     * @brief Pointer to the currently selected font.
     */
    gl_font* gl_font_manager::selected = NULL;

    /**
     * @brief Constructs a new gl_font_manager object.
     */
    gl_font_manager::gl_font_manager()
    {
    }

    /**
     * @brief Destroys the gl_font_manager object.
     */
    gl_font_manager::~gl_font_manager()
    {
    }

    /**
     * @brief Adds a font to the font manager.
     * @param name The name to associate with the font.
     * @param fnt Pointer to the font object.
     */
    void gl_font_manager::add_font(const std::string name, gl_font* fnt)
    {
        cg_fontMapper::iterator it = gl_font_manager::m_fontMap.find(name);
        if (it == m_fontMap.end())
            m_fontMap.insert(std::make_pair(name, fnt));
    }

    /**
     * @brief Retrieves a font by name.
     * @param name The name of the font.
     * @return Pointer to the font object, or NULL if not found.
     */
    gl_font* gl_font_manager::get_font(const std::string  name)
    {
        cg_fontMapper::iterator it = gl_font_manager::m_fontMap.find(name);
        if (it == m_fontMap.end())
            return NULL;
        return it->second;
    }

    /**
     * @brief Selects a font by id and sets it as the currently selected font.
     * @param id The name/id of the font.
     * @return Pointer to the selected font, or NULL if not found.
     */
    gl_font* gl_font_manager::select_font(const std::string id)
    {
        selected = get_font(id);
        return selected;
    }

    /**
     * @brief Creates a new 2D font and adds it to the manager.
     * @param fontName The name to associate with the font.
     * @param sysFontName The system font name.
     * @param fontSize The font size.
     * @return Pointer to the created font object.
     */
    gl_font* gl_font_manager::create_font(const char* fontName, const char* sysFontName, int fontSize)
    {
        gl_font* gfnt = new gl_font(sysFontName, fontSize);
        add_font(fontName, gfnt);
        return gfnt;
    }

    /**
     * @brief Creates a new 3D font and adds it to the manager.
     * @param fontName The name to associate with the font.
     * @param sysFontName The system font name.
     * @param fontSize The font size.
     * @param depth The extrusion depth for 3D fonts.
     * @return Pointer to the created font object.
     */
    gl_font* gl_font_manager::create_font(const char* fontName, const char* sysFontName, int fontSize, float depth)
    {
        gl_font* menu_font = new gl_font(sysFontName, fontSize, depth);
        add_font(fontName, menu_font);
        return menu_font;
    }

    /**
     * @brief Static instance of the font manager.
     */
    static gl_font_manager fManager;

    /**
     * @brief Gets the global font manager instance.
     * @return Reference to the font manager.
     */
    gl_font_manager& get_font_manager()
    {
        return fManager;
    }

    ////////////////////////////////////////////////////////////////
    /**
     * @brief Default constructor for gl_font.
     */
    gl_font::gl_font() : m_fontName(""), m_fontSize(0), m_depth(0)
    {
        init();
    }

    /**
     * @brief Constructs a 2D font.
     * @param fontName The system font name.
     * @param fontSize The font size.
     */
    gl_font::gl_font(const char* fontName, int fontSize) : m_fontName(fontName), m_fontSize(fontSize), m_depth(0)
    {
        init();
        build_2d(m_fontName.c_str(), m_fontSize);
    }

    /**
     * @brief Constructs a 3D font.
     * @param fontName The system font name.
     * @param fontSize The font size.
     * @param depth The extrusion depth for 3D fonts.
     */
    gl_font::gl_font(const char* fontName, int fontSize, float depth) : m_fontName(fontName), m_fontSize(fontSize), m_depth(depth)
    {
        init();
        build_3d(m_fontName.c_str(), m_fontSize, m_depth);
    }

    /**
     * @brief Destructor for gl_font. Releases font resources.
     */
    gl_font::~gl_font()
    {
        clear();
    }

    /**
     * @brief Clears the font display lists.
     */
    void gl_font::clear()
    {
        glDeleteLists(listBase, 256);
    }

    /**
     * @brief Initializes font properties and OpenGL function pointers.
     */
    void gl_font::init()
    {
        screenPos = ivec2(10);

        if (!glWindowPos2i)
            glWindowPos2i = (PFNGLWINDOWPOS2IPROC)wglGetProcAddress("glWindowPos2i");
    }

    /**
     * @brief Builds a 3D font using system font and OpenGL display lists.
     * @param name The system font name.
     * @param size The font size.
     * @param depth The extrusion depth.
     */
    void gl_font::build_3d(const char* name, int size, float depth)
    {
        HFONT hFont;
        HDC hDC = wglGetCurrentDC();
        listBase = glGenLists(256);      // create storage for 96 characters

        if (strcmp(name, "symbol") == 0)
        {
            hFont = CreateFontA(-size, 0, 0, 0, FW_BOLD, FALSE, FALSE, FALSE, SYMBOL_CHARSET,
                OUT_TT_PRECIS, CLIP_DEFAULT_PRECIS, ANTIALIASED_QUALITY,
                FF_DONTCARE | DEFAULT_PITCH, name);
        }
        else
        {
            hFont = CreateFontA(size, 0, 0, 0, FW_NORMAL, FALSE, FALSE, FALSE, ANSI_CHARSET,
                OUT_TT_PRECIS, CLIP_DEFAULT_PRECIS, ANTIALIASED_QUALITY,
                FF_DONTCARE | DEFAULT_PITCH, name);
        }

        if (!hFont)
            return;

        SelectObject(hDC, hFont);
        wglUseFontOutlines(hDC, 0, 255, listBase, (FLOAT)0.0f, (FLOAT)depth, WGL_FONT_POLYGONS, gmf);
    }

    /**
     * @brief Builds a 2D font using system font and OpenGL display lists.
     * @param name The system font name.
     * @param size The font size.
     */
    void gl_font::build_2d(const char* name, int size)
    {
        HFONT hFont;
        HDC hDC = wglGetCurrentDC();
        listBase = glGenLists(256);

        if (strcmp(name, "symbol") == 0)
        {
            hFont = CreateFontA(-size, 0, 0, 0, FW_BOLD, FALSE, FALSE, FALSE, SYMBOL_CHARSET,
                OUT_TT_PRECIS, CLIP_DEFAULT_PRECIS, ANTIALIASED_QUALITY,
                FF_DONTCARE | DEFAULT_PITCH, name);
        }
        else
        {
            hFont = CreateFontA(-size, 0, 0, 0, FW_NORMAL, FALSE, FALSE, FALSE, ANSI_CHARSET,
                OUT_TT_PRECIS, CLIP_DEFAULT_PRECIS, ANTIALIASED_QUALITY,
                FF_DONTCARE | DEFAULT_PITCH, name);
        }

        if (!hFont)
            return;

        SelectObject(hDC, hFont);
        wglUseFontBitmaps(hDC, 0, 256, listBase);
    }

    /**
     * @brief Renders text using the font and a shader, with alignment.
     * @param shader The shader to use for rendering.
     * @param alignment The text alignment (ALIGN_LEFT, ALIGN_CENTER, ALIGN_RIGHT).
     * @param str The format string for the text.
     * @param ... Additional arguments for the format string.
     */
    void gl_font::render(gl_shader* shader, int alignment, const char* str, ...)
    {
        float length = 0;
        char text[1512];
        va_list args;

        if ((str == NULL))
            return;

        va_start(args, str);
        vsprintf(text, str, args);
        va_end(args);

        // center the text
        for (unsigned int loop = 0; loop < (strlen(text)); loop++)	// find length of text
        {
            length += gmf[text[loop]].gmfCellIncX;		        // increase length by character's width
        }
        fvec3 offset(0, 0, 0);
        switch (alignment)
        {
        case ALIGN_CENTER:
            offset.x() = -length / 2;
            break;
        case ALIGN_LEFT:
            break;
        case ALIGN_RIGHT:
            offset.x() = -length;
            break;
        }

        // draw the text
        glPushAttrib(GL_LIST_BIT);
        glListBase(listBase);
        for (auto i = 0; i < strlen(text); ++i)
        {
            fmat4 ttm = tmat;
            translate_matrix(ttm, offset);
            // putting translation at the end treats the text as a unit
            fmat4 ob_matrix = rmat * smat * ttm;
            ob_matrix = ob_matrix.transpose();    // convert for OpenGL!
            shader->set_mat4("model", ob_matrix);
            glCallLists(1, GL_UNSIGNED_BYTE, &text[i]);
            offset.x() += gmf[text[i]].gmfCellIncX;
        }
        glPopAttrib();
    }

    /**
     * @brief Renders 2D text at the current screen position.
     * @param str The format string for the text.
     * @param ... Additional arguments for the format string.
     */
    void gl_font::render(const char* str, ...)
    {
        char text[1512];
        va_list args;

        if (str == NULL)
            return;

        va_start(args, str);
        vsprintf(text, str, args);
        va_end(args);

        glColor4f(color2d.x(), color2d.y(), color2d.z(), color2d.w());
        glWindowPos2i(screenPos.x(), screenPos.y());
        glPushAttrib(GL_LIST_BIT);
        glListBase(listBase);
        glCallLists((GLsizei)strlen(text), GL_UNSIGNED_BYTE, text);
        glPopAttrib();
    }
}