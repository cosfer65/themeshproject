#version 330 core

in vec3 FragPos;
in vec3 Normal;
in vec3 VertexColor;

uniform vec3 light_position;
uniform vec3 light_direction;   // must be normalized
uniform vec3 light_color;
uniform vec3 cameraPos;

uniform vec3 object_color;
uniform int useVertexColor; // 1 if using vertex color, 0 if using object color

out vec4 fragColor;

void main()
{
    vec3 colorToUse = useVertexColor == 1 ? VertexColor : object_color;
    float ambient_value = 0.65;
    vec3 ambient = ambient_value * light_color;
    
    vec3 viewPoint = vec3(0,0,0);//FragPos
    vec3 norm = normalize(Normal);
    vec3 lightDir = normalize(light_position - viewPoint);
    float diff = max(dot(norm, lightDir), 0.0) * 0.65;
    vec3 diffuse = diff * light_color;

    vec3 result = (ambient + diffuse) * colorToUse;
    
    fragColor = vec4(result, 1.0);
}
