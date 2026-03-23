#include "gl_graphics.h"

#include "app_core.h"
#include "glew.h"
#include "gl_timer.h"
#include "gl_shaders.h"

using namespace base_opengl;

void debug_out(const char* prompt, float var)
{
    char txt[200];
    sprintf(txt, "%s:%f", prompt, var);
    OutputDebugString(txt);
}

void __cdecl odprintf(const char* format, ...)
{
    char    buf[4096], * p = buf;
    va_list args;
    int     n;

    va_start(args, format);
    n = _vsnprintf(p, sizeof buf - 3, format, args); // buf-3 is room for CR/LF/NUL
    va_end(args);

    if (n < 0)
    {
        n = (int)(sizeof buf - 3);
        p[n] = '\0';
    }
    else if (n > (int)(sizeof buf - 3))
    {
        n = (int)(sizeof buf - 3);
        p[n] = '\0';
    }

    p += n;

    while (p > buf && isspace((unsigned char)p[-1]))
        *--p = '\0';

    *p++ = '\r';
    *p++ = '\n';
    *p = '\0';

    OutputDebugString(buf);
}