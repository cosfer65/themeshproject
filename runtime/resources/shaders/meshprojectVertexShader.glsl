#version 330 core

layout(location = 0) in vec3 aPos;
layout(location = 1) in vec3 aNormal;
layout(location = 2) in vec3 aVertexColor;

uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;

uniform int useVertexColor; // 1 if using vertex color, 0 if using object color
uniform vec3 object_color;

out vec3 FragPos;
out vec3 Normal;
out vec3 VertexColor;

void main()
{
    gl_Position = projection * view * model * vec4(aPos, 1.0);
    FragPos = vec3(model * vec4(aPos, 1.0));
    // The normal vector should be transformed by the inverse transpose of the model matrix
    Normal = mat3(transpose(inverse(model))) * aNormal;
    if (useVertexColor == 1) {
        VertexColor = aVertexColor;
    } else {
        VertexColor = object_color; // use object color
    }
}
