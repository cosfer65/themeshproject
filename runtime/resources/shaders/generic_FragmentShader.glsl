#version 330 core

uniform vec3 object_color;

out vec4 color; // drawing color for OpenGL to use

void main() {
	color = vec4(object_color, 1);
}