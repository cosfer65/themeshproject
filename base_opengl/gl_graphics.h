#ifndef __graphics_h__
#define __graphics_h__

#define _CRT_SECURE_NO_WARNINGS

// OpenGL
// use wrangler headers because gl headers are not updated in windows
#define GLEW_STATIC
#include "glew.h"
#include "wglew.h"

#pragma comment( lib, "opengl32.lib" )
#pragma comment( lib, "glu32.lib" )

#pragma comment( lib, "base_opengl.lib" )

#endif // __graphics_h__
