#ifndef __shader_h__
#define __shader_h__

#include <string>
#include <vector>

#include "vector.h"

using namespace base_math;

namespace base_opengl {
    /**
     * @brief Represents a shader described by a file on disk.
     *
     * This struct pairs an OpenGL shader type (e.g. GL_VERTEX_SHADER)
     * with a filename that contains the shader source. The filename is
     * not opened by this struct; it is a simple value object used by
     * shader loading routines.
     */
    struct shader_file {
        /**
         * @brief Construct a shader_file mapping a shader type to a filename.
         * @param _t OpenGL shader type (GLenum). Must be one of:
         *          GL_COMPUTE_SHADER, GL_VERTEX_SHADER,
         *          GL_TESS_CONTROL_SHADER, GL_TESS_EVALUATION_SHADER,
         *          GL_GEOMETRY_SHADER, GL_FRAGMENT_SHADER.
         * @param _s Path to the shader source file.
         */
        shader_file(GLenum _t, const std::string& _s) :shaderType(_t), shader_fname(_s) {}
        /// OpenGL shader type for this file.
        GLenum shaderType;
        /// Filesystem path to the shader source.
        std::string shader_fname;
    };

    /**
     * @brief Represents shader source code provided as a string.
     *
     * This struct pairs an OpenGL shader type with an in-memory
     * string containing shader source. It is used by compilation
     * helpers that accept source strings rather than filenames.
     */
    struct shader_source {
        /**
         * @brief Construct a shader_source mapping a shader type to source text.
         * @param _t OpenGL shader type (GLenum). Must be one of:
         *          GL_COMPUTE_SHADER, GL_VERTEX_SHADER,
         *          GL_TESS_CONTROL_SHADER, GL_TESS_EVALUATION_SHADER,
         *          GL_GEOMETRY_SHADER, GL_FRAGMENT_SHADER.
         * @param _s Shader source text.
         */
        shader_source(GLenum _t, const std::string& _s) :shaderType(_t), shader_src(_s) {}
        /// OpenGL shader type for this source.
        GLenum shaderType;
        /// Shader source code as a string.
        std::string shader_src;
    };

    /**
     * @brief Simple OpenGL shader program helper.
     *
     * Manages creation, compilation and use of an OpenGL program composed
     * from one or more shaders. Provides convenience setters for common
     * uniform types and a mechanism to remember and restore the previously
     * bound program.
     */
    class gl_shader {
        /// GL program handle owned by this object (0 means none).
        GLuint program;
        /**
         * @brief Previously bound GL program used by `use()` / `end()`.
         *
         * Stored as GLint because glGetIntegerv is used to query
         * GL_CURRENT_PROGRAM.
         */
        GLint old_program;
        /// List of shader files to be loaded when `load()` is called.
        std::vector<shader_file> sh_file_list;
    public:
        /**
         * @brief Construct an empty shader helper.
         *
         * Initializes internal state; no GL resources are created until
         * `load`/`compile` are invoked.
         */
        gl_shader() :program(0) {
        }

        /**
         * @brief Destructor releases owned GL program if present.
         *
         * Calls glDeleteProgram on the stored program handle.
         */
        ~gl_shader()
        {
            if (program)
                glDeleteProgram(program);
        }

        /**
         * @brief Add a shader file to the internal list.
         * @param _t OpenGL shader type (GLenum).
         * @param _s Path to the shader file to add.
         * @return New count of entries in the internal file list.
         *
         * The actual file is not read until `load()` is called.
         */
        size_t add_file(GLenum _t, const std::string& _s) {
            sh_file_list.push_back(shader_file(_t, _s));
            return sh_file_list.size();
        }

        /**
         * @brief Load shader files previously added via `add_file` and link them into a program.
         * @return The GL program handle (or 0 on failure).
         *
         * Implementation should read each file in `sh_file_list`, compile shaders,
         * attach them to a program and link. On success `program` will be set.
         */
        GLuint load();

        /**
         * @brief Compile and link shaders from a vector of `shader_source` pointers.
         * @param code Vector of pointers to shader_source describing shader types and code.
         * @return The GL program handle (or 0 on failure).
         */
        GLuint compile(const std::vector<shader_source*>& code);

        /**
         * @brief Convenience loader for a common vertex/fragment/geometry set.
         * @param vertex_file Path to vertex shader file.
         * @param fragment_file Path to fragment shader file.
         * @param geometry_file Optional path to geometry shader file (nullable).
         * @return The GL program handle (or 0 on failure).
         */
        GLuint load(const char* vertex_file, const char* fragment_file, const char* geometry_file = nullptr);

        /**
         * @brief Compile program from in-memory source strings.
         * @param vertex_source Vertex shader source string.
         * @param fragment_source Fragment shader source string.
         * @param geom_source Optional geometry shader source string (nullable).
         * @return The GL program handle (or 0 on failure).
         */
        GLuint compile(const char* vertex_source, const char* fragment_source, const char* geom_source = nullptr);

        /**
         * @brief Set common default uniform values after binding.
         *
         * Called automatically by `use()`; sets `use_vertex_color` to 0
         * unless overridden by caller.
         */
        void set_defaults() {
            set_int("use_vertex_color", 0);
        }

        /**
         * @brief Bind this shader program for use.
         *
         * Stores the previously bound program in `old_program`, calls
         * glUseProgram with this object's program handle, then applies
         * `set_defaults`.
         */
        void use() {
            glGetIntegerv(GL_CURRENT_PROGRAM, &old_program);
            glUseProgram(program);
            set_defaults();
        }
        /**
         * @brief Restore the previously bound shader program saved by `use()`.
         */
        void end() {
            glUseProgram(old_program);
        }
        /**
         * @brief Set an integer uniform.
         * @param name Uniform name in the shader.
         * @param value Integer value to set.
         *
         * Uses glGetUniformLocation(program, name) to locate the uniform
         * and glUniform1i to set it.
         */
        void set_int(const std::string& name, int value) const {
            glUniform1i(glGetUniformLocation(program, name.c_str()), value);
        }
        /**
         * @brief Set a float uniform.
         * @param name Uniform name in the shader.
         * @param value Float value to set.
         */
        void set_float(const std::string& name, float value) const {
            glUniform1f(glGetUniformLocation(program, name.c_str()), value);
        }
        /**
         * @brief Set a 4x4 matrix uniform from a raw float pointer.
         * @param name Uniform name in the shader.
         * @param mat Pointer to 16 floats (column-major) used to set the matrix.
         */
        void set_mat4(const std::string& name, const float* mat) const {
            glUniformMatrix4fv(glGetUniformLocation(program, name.c_str()), 1, GL_FALSE, mat);
        }
        /**
         * @brief Set a 4x4 matrix uniform from a `fmat4`.
         * @param name Uniform name in the shader.
         * @param mat Reference to a `fmat4` object (must be convertible to float*).
         */
        void set_mat4(const std::string& name, fmat4& mat) const {
            glUniformMatrix4fv(glGetUniformLocation(program, name.c_str()), 1, GL_FALSE, (float*)mat);
        }
        /**
         * @brief Set a fvec4 uniform from a `fvec4`.
         * @param name Uniform name in the shader.
         * @param value Reference to a `fvec4` (must be convertible to float*).
         */
        void set_vec4(const std::string& name, fvec4& value) const {
            glUniform4fv(glGetUniformLocation(program, name.c_str()), 1, (float*)value);
        }
        /**
         * @brief Set a fvec4 uniform from a raw float pointer.
         * @param name Uniform name in the shader.
         * @param value Pointer to 4 floats.
         */
        void set_vec4(const std::string& name, const float* value) const {
            glUniform4fv(glGetUniformLocation(program, name.c_str()), 1, value);
        }
        /**
         * @brief Set a fvec3 uniform from a `fvec3`.
         * @param name Uniform name in the shader.
         * @param value Reference to a `fvec3` (must be convertible to float*).
         */
        void set_vec3(const std::string& name, fvec3& value) const {
            glUniform3fv(glGetUniformLocation(program, name.c_str()), 1, (float*)value);
        }
        /**
         * @brief Set a fvec3 uniform from a raw float pointer.
         * @param name Uniform name in the shader.
         * @param value Pointer to 3 floats.
         */
        void set_vec3(const std::string& name, const float* value) const {
            glUniform3fv(glGetUniformLocation(program, name.c_str()), 1, value);
        }
    };
}

#endif // __shader_h__
