#include <stdlib.h>
#include <string.h>

#include "drutils.h"

char 
KDccToDRDcc(KDcc dcc)
{
    switch (dcc)
    {
        case DCC_30_360:    return '3'; break;
        case DCC_ACT_360:   return '0'; break;
        case DCC_ACT_365:   return '5'; break;
        case DCC_ACT_ACT:   return 'A'; break;
    }
    return '?';
}

void* 
ArrayNewCopy(void const* src,
             size_t      size)
{
    void* ptr = malloc(size);
    if (ptr)
        memcpy(ptr, src, size);
    return ptr;
}

void 
SafeFree(void* ptr)
{
    if (ptr)
        free(ptr);
}

