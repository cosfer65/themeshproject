#include "gl_graphics.h"
#include "gl_shaders.h"
#include <fstream>
#include <sstream>

namespace base_opengl {
    /**
     * @brief Loads shader source files listed in sh_file_list, compiles them, and links them into a shader program.
     * 
     * This function reads each shader file specified in sh_file_list, creates a shader_source object for each,
     * compiles and links them into a single OpenGL shader program, and returns the program ID.
     * Any errors during file reading are reported (on Windows) via OutputDebugString.
     * 
     * @return GLuint The OpenGL program ID of the linked shader program.
     */
    GLuint gl_shader::load()
    {
        std::vector<shader_source*> code;

        for (auto it = sh_file_list.begin(); it != sh_file_list.end(); ++it)
        {
            std::ifstream shaderFile;
            shaderFile.exceptions(std::ifstream::failbit | std::ifstream::badbit);
            try
            {
                shaderFile.open((*it).shader_fname);
                std::stringstream shaderStream;
                shaderStream << shaderFile.rdbuf();
                shaderFile.close();
                shader_source* src = new shader_source((*it).shaderType, shaderStream.str());
                code.push_back(src);
            }
            catch (...)
            {
#ifdef WIN32
                OutputDebugString("error");
                OutputDebugString((*it).shader_fname.c_str());
#endif
            }
        }
        program = compile(code);
        for (auto i = code.begin(); i != code.end(); ++i)
        {
            delete* i;
        }

        return program;
    }

    /**
     * @brief Compiles and links a set of shader sources into an OpenGL shader program.
     * 
     * This function creates, compiles, and attaches each shader in the provided vector, then links them into a program.
     * 
     * @param code Vector of pointers to shader_source objects, each containing shader type and source code.
     * @return GLuint The OpenGL program ID of the linked shader program.
     */
    GLuint gl_shader::compile(const std::vector<shader_source*>& code)
    {
        std::vector<GLuint> ishaders;
        for (auto it = code.begin(); it != code.end(); ++it)
        {
            GLuint _shader = glCreateShader((*it)->shaderType);
            const char* src = (*it)->shader_src.c_str();
            glShaderSource(_shader, 1, &src, NULL);
            glCompileShader(_shader);
            ishaders.push_back(_shader);
        }
        program = glCreateProgram();
        for (auto it = ishaders.begin(); it != ishaders.end(); ++it)
            glAttachShader(program, *it);
        glLinkProgram(program);

        return program;
    }

    /**
     * @brief Loads and compiles vertex, fragment, and optionally geometry shader files into a shader program.
     * 
     * This function adds the specified shader files to the file list, then loads and compiles them into a program.
     * 
     * @param vertex_file Path to the vertex shader file.
     * @param fragment_file Path to the fragment shader file.
     * @param geometry_file Optional path to the geometry shader file (default is nullptr).
     * @return GLuint The OpenGL program ID of the linked shader program.
     */
    GLuint gl_shader::load(const char* vertex_file, const char* fragment_file, const char* geometry_file /*= nullptr*/)
    {
        add_file(GL_VERTEX_SHADER, vertex_file);
        add_file(GL_FRAGMENT_SHADER, fragment_file);
        if (geometry_file)
            add_file(GL_GEOMETRY_SHADER, geometry_file);
        program = load();
        return program;
    }

    /**
     * @brief Compiles and links provided shader source strings into an OpenGL shader program.
     * 
     * This function creates shader_source objects from the provided source strings, compiles and links them.
     * 
     * @param vertex_source Source code for the vertex shader.
     * @param fragment_source Source code for the fragment shader.
     * @param geom_source Optional source code for the geometry shader (default is NULL).
     * @return GLuint The OpenGL program ID of the linked shader program.
     */
    GLuint gl_shader::compile(const char* vertex_source, const char* fragment_source, const char* geom_source/* = NULL*/)
    {
        std::vector<shader_source*> code;
        if (vertex_source)
            code.push_back(new shader_source(GL_VERTEX_SHADER, vertex_source));
        code.push_back(new shader_source(GL_FRAGMENT_SHADER, fragment_source));
        if (geom_source)
            code.push_back(new shader_source(GL_GEOMETRY_SHADER, geom_source));
        program = compile(code);
        for (auto i = code.begin(); i != code.end(); ++i)
        {
            delete* i;
        }
        return program;
    }
}