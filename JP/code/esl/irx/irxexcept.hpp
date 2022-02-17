/*
***************************************************************************
** HEADER FILE: irxexcept.h
**
** C++ Error handling for the credit & rates exotics.
***************************************************************************
*/

#ifndef IRX_EXCEPT_H
#define IRX_EXCEPT_H

#include <stdarg.h>          /* defines va_list */
#include "error.h"
#include <exception>

class irxException : public std::exception
{
public:
        irxException(){}
        irxException(const char *fmt, ...); 

        virtual ~irxException() throw () {}
};




/*
 * Useful macros.
 */

#define IF_FAILED_THROW(statement)  \
    {if (statement != SUCCESS) throw irxException("Condition `%s' failed.\n", \
     #statement);}

#define ASSERT_OR_THROW(statement)  \
    {if (!(statement)) throw irxException("Assertion `%s' failed.\n",\
     #statement);}

#define READ_DATA_CPP(type,ptr,str)        \
    { if (irxFScanVType(fp, type, (void*) ptr) != SUCCESS) \
          { throw irxException("%s: can't read %s.\n", routine, str);}}

#endif

