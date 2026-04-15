#version 330 core

uniform vec4 object_color; // object color

// drawing color for OpenGL to use
out vec4 color;
void main() {
	color = object_color;
}