#ifndef __light_h__
#define __light_h__

#include "vector.h"
#include "gl_shaders.h"

using namespace base_math;

namespace base_opengl {

    /**
     * @class gl_light
     * @brief Represents a light source for OpenGL rendering.
     *
     * This class encapsulates the properties and behavior of a light source,
     * including its type, position or direction, and lighting components
     * (ambient, diffuse, specular). It provides methods to set and get these
     * properties, and to apply them to a shader.
     */
    class gl_light {
    protected:
        fvec3 pos_or_dir;    ///< Light position (for point/spot) or direction (for directional light).
        fvec3 ambient;       ///< Ambient color component of the light.
        fvec3 diffuse;       ///< Diffuse color component of the light.
        fvec3 specular;      ///< Specular color component of the light.
    public:
        /**
         * @enum light_type
         * @brief Enumerates the types of supported lights.
         */
        enum light_type { SPOTLIGHT = 1, DIRLIGHT };

        light_type _type;   ///< The type of the light.

        /**
         * @brief Constructs a gl_light with optional type.
         * @param _t The type of the light (default: SPOTLIGHT).
         */
        gl_light(light_type _t = SPOTLIGHT)
            : pos_or_dir(fvec3(8, 2, 8)),
              ambient(fvec3(0.4f)),
              diffuse(fvec3(0.7f)),
              specular(fvec3(0.95f)),
              _type(_t) {
        }

        /**
         * @brief Destructor for gl_light.
         */
        ~gl_light() {
        }

        /**
         * @brief Sets the ambient color component.
         * @param a The ambient color.
         */
        void set_ambient(const fvec3& a) {
            ambient = a;
        }

        /**
         * @brief Gets the ambient color component.
         * @return The ambient color.
         */
        fvec3 get_ambient() {
            return ambient;
        }

        /**
         * @brief Sets the diffuse color component.
         * @param d The diffuse color.
         */
        void set_diffuse(const fvec3& d) {
            diffuse = d;
        }

        /**
         * @brief Gets the diffuse color component.
         * @return The diffuse color.
         */
        fvec3 get_diffuse() {
            return diffuse;
        }

        /**
         * @brief Sets the specular color component.
         * @param s The specular color.
         */
        void set_specular(const fvec3& s) {
            specular = s;
        }

        /**
         * @brief Gets the specular color component.
         * @return The specular color.
         */
        fvec3 get_specular() {
            return specular;
        }

        /**
         * @brief Sets the position of the light (for point/spot lights).
         * @param p The position vector.
         */
        void set_position(const fvec3& p) {
            pos_or_dir = p;
        }

        /**
         * @brief Gets the position of the light.
         * @return The position vector.
         */
        fvec3 get_position() {
            return pos_or_dir;
        }

        /**
         * @brief Sets the direction of the light (for directional lights).
         * @param d The direction vector.
         */
        void set_direction(const fvec3& d) {
            pos_or_dir = d;
            pos_or_dir.normalize();
        }

        /**
         * @brief Gets the direction of the light.
         * @return The direction vector.
         */
        fvec3 get_direction() {
            return pos_or_dir;
        }

        /**
         * @brief Applies the light properties to the given shader.
         * @param shdr Pointer to the shader to apply the light to.
         */
        void apply(gl_shader* shdr) {
            shdr->set_int("light.type", (int)_type);
            shdr->set_vec3("light.pos_or_dir", pos_or_dir);
            shdr->set_vec3("light.ambient", ambient);
            shdr->set_vec3("light.diffuse", diffuse);
            shdr->set_vec3("light.specular", specular);
        }
    };
}
#endif // __light_h__
