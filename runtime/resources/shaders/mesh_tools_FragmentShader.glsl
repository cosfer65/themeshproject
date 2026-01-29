/*
    Fragment Shader for Mesh Tools

    This shader implements per-fragment lighting using the Phong reflection model. It supports two types of light sources: spotlight and directional light, as indicated by the `type` field in the `lightsource` struct.

    Uniforms:
    - `light`: Contains the properties of the active light source, including type, position/direction, and ambient/diffuse/specular components.
    - `cameraPos`: The position of the camera/viewer in world space.

    Inputs:
    - `vertex_color`: The color passed from the vertex shader, representing the object's base color.
    - `Normal`: The interpolated surface normal at the fragment.
    - `pos`: The world-space position of the fragment.

    Output:
    - `color`: The final color output for the fragment.

    Main Lighting Calculations:
    - The normal vector is normalized for accurate lighting.
    - The direction to the light is computed and normalized.
    - The view direction (towards the camera) is calculated.
    - Specular reflection is computed using the Phong model, with a shininess exponent of 64.
    - Diffuse reflection is calculated as the dot product of the normal and light direction, clamped to zero.
    - Ambient lighting is taken directly from the light source.
    - The final color is the sum of ambient, diffuse, and specular components, each modulated by the vertex color.

    This shader is suitable for rendering lit 3D objects with basic material properties and a single light source.
*/
#version 330 core

struct lightsource {
    int type;             // SPOTLIGHT=1, DIRLIGHT=2
    vec3 pos_or_dir;      // light location
    vec3 ambient;         // the ambience property of the light
    vec3 diffuse;         // the diffuse property of the light
    vec3 specular;        // the specular property of the light
};
uniform lightsource light;

uniform vec3 cameraPos;   // viewer location

in vec4 vertex_color;     // vertex/object coloring (determined in vertex shader)
in vec3 Normal;           // surface normal
in vec3 pos;              // drawing position

out vec4 color;           // resulting drawing color

void main(){
    vec3 norm = normalize(Normal);
    vec3 lightDir = normalize(light.pos_or_dir - pos);
    
    vec3 viewDir = normalize(cameraPos - pos);
    vec3 reflectDir = reflect(lightDir, norm);
    float spec = pow(max(dot(viewDir, reflectDir), 0.0), 64);
    vec3 specular = spec * light.specular;

    float diff = max(dot(norm, lightDir), 0.0);
    vec3 diffuse = diff * light.diffuse;

    vec3 ambient = light.ambient;

	vec4 a_color = vec4((ambient),1) * vertex_color;
	vec4 d_color = vec4((diffuse),1) * vertex_color;
	vec4 s_color = vec4((specular),1) * vertex_color;
    color = a_color + d_color + s_color;
}
