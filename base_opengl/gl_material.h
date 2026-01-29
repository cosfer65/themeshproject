#ifndef __material_h__
#define __material_h__

#include "vector.h"
#include "gl_shaders.h"

using namespace base_math;

namespace base_opengl {

    /**
     * @class cg_material
     * @brief Represents a material for use in OpenGL rendering, encapsulating color and texture properties.
     *
     * The cg_material class manages ambient, diffuse, and specular color components, shininess, and optional
     * diffuse and specular texture maps. It provides methods to set and retrieve these properties, and to
     * apply them to a shader for rendering.
     */
    class cg_material {
        fvec3 ambient;        ///< Ambient color component of the material.
        fvec3 diffuse;        ///< Diffuse color component of the material.
        fvec3 specular;       ///< Specular color component of the material.
        float shine;         ///< Shininess coefficient for specular highlights.

        int diffuse_map;     ///< OpenGL texture ID for the diffuse map, or -1 if not used.
        int diffuse_index;   ///< Texture unit index for the diffuse map, or -1 if not used.
        int specular_map;    ///< OpenGL texture ID for the specular map, or -1 if not used.
        int specular_index;  ///< Texture unit index for the specular map, or -1 if not used.
    public:
        /**
         * @brief Constructs a cg_material with default color and texture values.
         */
        cg_material() : ambient(fvec3(0.2f, 0.2f, 0.2f)),
            diffuse(fvec3(.95f, .95f, .95f)), specular(fvec3(.5f, 0.5f, 0.5f)),
            shine(128.f), diffuse_map(-1), diffuse_index(-1),
            specular_map(-1), specular_index(-1) {
        }

        /**
         * @brief Destructor for cg_material.
         */
        ~cg_material() {
        }

        /**
         * @brief Gets the ambient color component.
         * @return The ambient color as a fvec3.
         */
        fvec3 get_ambient() const {
            return ambient;
        }

        /**
         * @brief Gets the diffuse color component.
         * @return The diffuse color as a fvec3.
         */
        fvec3 get_diffuse() const {
            return diffuse;
        }

        /**
         * @brief Gets the specular color component.
         * @return The specular color as a fvec3.
         */
        fvec3 get_specular() const {
            return specular;
        }

        /**
         * @brief Sets the ambient color component.
         * @param a The new ambient color.
         */
        void set_ambient(const fvec3& a) {
            ambient = a;
        }

        /**
         * @brief Sets the ambient color component using individual float values.
         * @param a0 Red component.
         * @param a1 Green component.
         * @param a2 Blue component.
         */
        void set_ambient(float a0, float a1, float a2) {
            ambient = fvec3(a0, a1, a2);
        }

        /**
         * @brief Sets the diffuse color component.
         * @param d The new diffuse color.
         */
        void set_diffuse(const fvec3& d) {
            diffuse = d;
        }

        /**
         * @brief Sets the diffuse color component using individual float values.
         * @param a0 Red component.
         * @param a1 Green component.
         * @param a2 Blue component.
         */
        void set_diffuse(float a0, float a1, float a2) {
            diffuse = fvec3(a0, a1, a2);
        }

        /**
         * @brief Sets the diffuse texture map and its texture unit index.
         * @param texture OpenGL texture ID.
         * @param active_texture Texture unit index.
         */
        void set_diffuse(int texture, int active_texture) {
            diffuse_map = texture;
            diffuse_index = active_texture;
        }

        /**
         * @brief Sets the specular color component.
         * @param s The new specular color.
         */
        void set_specular(const fvec3& s) {
            specular = s;
        }

        /**
         * @brief Sets the specular color component using individual float values.
         * @param a0 Red component.
         * @param a1 Green component.
         * @param a2 Blue component.
         */
        void set_specular(float a0, float a1, float a2) {
            specular = fvec3(a0, a1, a2);
        }

        /**
         * @brief Sets the specular texture map and its texture unit index.
         * @param texture OpenGL texture ID.
         * @param active_texture Texture unit index.
         */
        void set_specular(int texture, int active_texture) {
            specular_map = texture;
            specular_index = active_texture;
        }

        /**
         * @brief Sets the shininess coefficient.
         * @param s The new shininess value.
         */
        void set_shine(float s) {
            shine = s;
        }

        /**
         * @brief Applies the material properties to the given shader.
         *
         * Sets the shader uniforms for ambient, diffuse, specular, and shininess.
         * If texture maps are set, binds them to the appropriate texture units and sets the corresponding shader uniforms.
         *
         * @param shdr Pointer to the shader to apply the material to.
         */
        void apply(gl_shader* shdr)
        {
            shdr->set_vec3("mat.ambient", ambient);
            shdr->set_vec3("mat.diffuse", diffuse);
            shdr->set_vec3("mat.specular", specular);
            shdr->set_float("mat.shine", shine);
            if (diffuse_index >= 0)
            {
                shdr->set_int("mat.diffuse_map", diffuse_index); // GL_TEXTURE0=0...
                glActiveTexture(GL_TEXTURE0 + diffuse_index);
                glBindTexture(GL_TEXTURE_2D, diffuse_map);
            }
            if (specular_index >= 0)
            {
                shdr->set_int("mat.specular_map", specular_index); // GL_TEXTURE0=0...
                glActiveTexture(GL_TEXTURE0 + specular_index);
                glBindTexture(GL_TEXTURE_2D, specular_map);
            }
        }
    };
}

#endif // __material_h__
