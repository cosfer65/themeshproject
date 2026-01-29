/*
    Vertex Shader for Mesh Tools

    This GLSL vertex shader transforms mesh vertices from model space to clip space and assigns colors based on either a uniform object color or per-vertex color. It also computes the transformed normal vector for lighting calculations.

    Inputs:
        - aPos (location = 0): Vertex position in model space.
        - aNormal (location = 1): Vertex normal in model space.
        - aVertexCol (location = 2): Per-vertex color (used if object_or_vertex_color == 1).

    Uniforms:
        - object_or_vertex_color: Determines color source (0 = object color, 1 = per-vertex color).
        - object_color: Uniform color for the entire object (used if object_or_vertex_color == 0).
        - camera: Camera transformation matrix (model-view-projection).
        - model: Model transformation matrix.

    Outputs:
        - vertex_color: Final color for the vertex, passed to the fragment shader.
        - Normal: Transformed normal vector in world space.
        - pos: Transformed vertex position in world space.

    Main Function:
        - Transforms the vertex position by the model matrix and then by the camera matrix.
        - Selects the vertex color based on the object_or_vertex_color uniform.
        - Calculates the world-space position and normal for use in lighting and shading.
*/

#version 330 core

layout(location = 0) in vec3 aPos;
layout(location = 1) in vec3 aNormal;

uniform int object_or_vertex_color; // 0= object color, 1= vertex color
// used when per vertex color is used
layout(location = 2) in vec3 aVertexCol;
// used when uniform object color is selected
uniform vec4 object_color;    // object color

uniform mat4 camera;
uniform mat4 model;

out vec4 vertex_color;
out vec3 Normal;
out vec3 pos;

void main()
{
    vec4 p = model * vec4(aPos, 1);
    gl_Position = camera * p;
    
    if (object_or_vertex_color == 0)  // use object color
        vertex_color = object_color;
    else                           // use per vertex color
        vertex_color = vec4(aVertexCol,1);

    pos = vec3(p);
    Normal = mat3(transpose(inverse(model))) * aNormal;
}
