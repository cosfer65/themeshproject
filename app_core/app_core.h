#ifndef __app_core_h__
#define __app_core_h__

void debug_out(const char* prompt, float var);
void __cdecl odprintf(const char* format, ...);

#pragma comment( lib, "app_core.lib" )

#endif // __app_core_h__
