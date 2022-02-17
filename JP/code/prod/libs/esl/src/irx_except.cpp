/*
***************************************************************************
** SOURCE FILE: irx_except.c
**
** C++ error handling.
**
***************************************************************************
*/

#include "irx/error.h"
#include "irx/irxexcept.hpp"


/*f
***************************************************************************
** C++ exception handler.
** calling irxErrorV internally.
***************************************************************************
*/
irxException::irxException(const char *format, ...)
{
    va_list args;

    va_start(args, format);
    irxErrorV (format, args);
    va_end(args);
}
