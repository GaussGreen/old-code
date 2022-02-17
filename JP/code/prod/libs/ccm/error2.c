#include <stdio.h>
#include <stdarg.h>

typedef void (*ErrFuncPtr)(const char* , va_list);
ErrFuncPtr errCallBack = NULL;

void DR_Error(const char* format, ...)
{
	va_list args;
	va_start(args, format);
    if (errCallBack)
    {
        (*errCallBack)(format, args);
    }
    else
    {
        vfprintf(stdout,format,args);
	    fprintf(stdout,"\n");
	    fflush(stdout);
    }
	va_end(args);
}


