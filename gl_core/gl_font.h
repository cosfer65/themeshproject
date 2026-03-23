/* fonts for text output */
#ifndef __font_h__
#define __font_h__
#include "gl_prim.h"
#include <map>

#define ALIGN_CENTER 0
#define ALIGN_LEFT   1
#define ALIGN_RIGHT  2

namespace base_opengl {
    class gl_shader;

    /**
     * @class gl_font
     * @brief Simple font class supporting 2D and 3D fonts.
     * 
     * Inherits from gl_prim. Provides methods for rendering text in both 2D and 3D,
     * with support for alignment, color, and positioning. Handles font creation, cleanup,
     * and display list management.
     */
    class gl_font : public gl_prim {
        unsigned int listBase;              ///< Display list base for font glyphs.
        GLYPHMETRICSFLOAT gmf[256];         ///< Holds orientation and placement metrics for each glyph.

        unsigned int texID;                 ///< Font texture ID.
        unsigned int callList;              ///< Font display list.

        fvec4 color2d;                       ///< RGBA color for 2D text rendering.
        ivec2 screenPos;                    ///< Screen position for 2D text (in pixels).
        std::string m_fontName;             ///< Name of the font.
        int m_fontSize;                     ///< Font size.
        float m_depth;                      ///< Depth for 3D fonts.

        /**
         * @brief Initializes the font object.
         * 
         * Internal helper for setting up font resources.
         */
        void init();

    public:
        /**
         * @brief Constructs a 3D font with specified depth.
         * @param fontName Name of the font.
         * @param fontSize Size of the font.
         * @param depth Depth for 3D extrusion.
         */
        gl_font(const char* fontName, int fontSize, float depth);

        /**
         * @brief Constructs a 2D font.
         * @param fontName Name of the font.
         * @param fontSize Size of the font.
         */
        gl_font(const char* fontName, int fontSize);

        /**
         * @brief Default constructor.
         */
        gl_font();

        /**
         * @brief Destructor. Cleans up font resources.
         */
        ~gl_font();

        /**
         * @brief Renders 3D text with alignment.
         * @param shader Shader to use for rendering.
         * @param alignment Text alignment (ALIGN_LEFT, ALIGN_CENTER, ALIGN_RIGHT).
         * @param str Format string (printf-style).
         * @param ... Additional arguments for formatting.
         */
        void render(gl_shader* shader, int alignment, const char* str, ...);

        /**
         * @brief Renders 2D text at the current position.
         * @param str Format string (printf-style).
         * @param ... Additional arguments for formatting.
         */
        void render(const char* str, ...);

        /**
         * @brief Builds a 3D font.
         * @param name Name of the font.
         * @param size Font size.
         * @param depth Depth for 3D extrusion.
         */
        void build_3d(const char* name, int size, float depth);

        /**
         * @brief Builds a 2D font.
         * @param name Name of the font.
         * @param size Font size.
         */
        void build_2d(const char* name, int size);

        /**
         * @brief Explicitly clears and releases font resources.
         */
        void clear();

        /**
         * @brief Sets the 2D text position in pixels.
         * @param x X coordinate.
         * @param y Y coordinate.
         */
        void set_position(int x, int y) { screenPos = ivec2(x, y); }

        /**
         * @brief Sets the color for 2D text rendering.
         * @param col RGBA color vector.
         */
        void set_color(const fvec4& col) { color2d = col; }
    };

    /**
     * @typedef cg_fontMapper
     * @brief Map of font identifiers to gl_font pointers.
     */
    typedef std::map<std::string, gl_font*> cg_fontMapper;

    /**
     * @class gl_font_manager
     * @brief Manages font objects and provides font selection and creation.
     * 
     * Handles a map of fonts, allows adding, retrieving, and selecting fonts.
     * Supports creation of 2D and 3D fonts.
     */
    class gl_font_manager
    {
        static cg_fontMapper m_fontMap;     ///< Static map of font identifiers to font objects.
        static gl_font* selected;           ///< Currently selected font.
    public:
        /**
         * @brief Constructor.
         */
        gl_font_manager();

        /**
         * @brief Destructor.
         */
        ~gl_font_manager();

        /**
         * @brief Adds a font to the manager.
         * @param id Identifier for the font.
         * @param fnt Pointer to the font object.
         */
        void add_font(const std::string id, gl_font* fnt);

        /**
         * @brief Retrieves a font by its identifier.
         * @param id Identifier for the font.
         * @return Pointer to the font object, or nullptr if not found.
         */
        gl_font* get_font(const std::string id);

        /**
         * @brief Selects a font by its identifier.
         * @param id Identifier for the font.
         * @return Pointer to the selected font, or nullptr if not found.
         */
        gl_font* select_font(const std::string id);

        /**
         * @brief Creates a 2D font and adds it to the manager.
         * @param fontName Font identifier.
         * @param sysFontName System font name.
         * @param fontSize Font size.
         * @return Pointer to the created font.
         */
        gl_font* create_font(const char* fontName, const char* sysFontName, int fontSize);

        /**
         * @brief Creates a 3D font and adds it to the manager.
         * @param fontName Font identifier.
         * @param sysFontName System font name.
         * @param fontSize Font size.
         * @param depth Depth for 3D extrusion.
         * @return Pointer to the created font.
         */
        gl_font* create_font(const char* fontName, const char* sysFontName, int fontSize, float depth);

        /**
         * @brief Returns the currently selected font.
         * @return Pointer to the selected font.
         */
        gl_font* select_font() {
            return selected;
        }
    };

    /**
     * @brief Returns a reference to the global font manager instance.
     * @return Reference to gl_font_manager.
     */
    gl_font_manager& get_font_manager();
}
#endif // __FONT_H__
