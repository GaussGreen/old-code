/*
***************************************************************************
** SOURCE FILE: error.c
**
** Error handling for the CX product set.
**
** $Header$
***************************************************************************
*/

#include "irx/error.h"

#include <stdio.h>
#include <stdarg.h>

#include "irx/platdep.h"
#include "irx/macros.h"

typedef struct 
{
    IrxTErrorCallBackV errorCallBackV;
    IrxTErrorCallBack  errorCallBack;
} ERROR_GLOBAL;

static ERROR_GLOBAL global =
{
    NULL,
    NULL
};

#define ERRBUFFSIZE    1024L           /* size of error message      */
static char gErrMsg[ERRBUFFSIZE]="";   /* error buffer for call back */


/*f
***************************************************************************
** Error handler.
**
** By default this prints to stderr. However by setting a call back routine
** you can intercept the messages with your own message handler.
**
** The message printed should be terminated with a '\n' as this function
** does not add it.
***************************************************************************
*/
void irxError(const char *format, ...)
{
    va_list args;
    va_start(args, format);
    irxErrorV (format, args);
    va_end(args);
    return;
}

/*f
***************************************************************************
** Error handler using va_list instead of ...
**
** By default this prints to stderr. However by setting a call back routine
** you can intercept the messages with your own message handler.
**
** The message printed should be terminated with a '\n' as this function
** does not add it.
***************************************************************************
*/
void irxErrorV(const char *format, va_list args)
{
    int doStdMsg = TRUE;

    if (global.errorCallBack != NULL) 
    {
        char errorBuff[ERRBUFFSIZE];
        vsnprintf(errorBuff, ERRBUFFSIZE, format, args);
        if (global.errorCallBack(errorBuff))
            doStdMsg = FALSE;
    }

    if (global.errorCallBackV != NULL) 
    {
        if (global.errorCallBackV(format, args))
            doStdMsg = FALSE;
    }

    if (doStdMsg)
    {
        vfprintf(stderr, format, args);
        fflush(stderr);
    }
}

/*f
***************************************************************************
** Handles the message that a routine failed. Designed so that this message
** is consistent throughout the library.
***************************************************************************
*/
int irxErrorFailure(const char* routine)
{
    irxError ("%s: Failed.\n", routine);
    return FAILURE;
}

/*f
***************************************************************************
** Sets an error callback routine.
**
** This callback routine takes the formatted string.
** If the callback returns FALSE, then the message is no longer printed
** to stderr as well.
**
** See also irxSetErrorCallBackV
***************************************************************************
*/
void irxSetErrorCallBack(IrxTErrorCallBack callBack)
{
    global.errorCallBack = callBack;
}

/*f
***************************************************************************
** Sets an error callback routine.
**
** This callback routine takes the format string and va_list arguments.
** If the callback returns FALSE, then the message is no longer printed
** to stderr as well.
**
** See also irxSetErrorCallBack
***************************************************************************
*/
void irxSetErrorCallBackV(IrxTErrorCallBackV callBackV)
{
    global.errorCallBackV = callBackV;
}



/*f
***************************************************************************
**
**       Store the error message in a string buffer
**
***************************************************************************
*/
void irxErrorCallBack(const char* newErrMsg)
{
    long    oldErrSize, newErrSizeMax;

    /* Existing error msg size */
    oldErrSize = strlen(gErrMsg);

    /* Max allowed size of new error msg */
    newErrSizeMax = ERRBUFFSIZE - oldErrSize;

    /* Append the new error msg to the existing error msg in the buffer, 
     * up to the full buffer including \0.  Do this only if existing 
     * error bufffer is not full ((including \0).
     */
    if (newErrSizeMax > 1)
        strncat(gErrMsg, newErrMsg, newErrSizeMax );

    return;
}



/*f
***************************************************************************
**
**       Retrieve the error message in a string buffer
**       Store the error message in a string buffer
**
***************************************************************************
*/
const char* irxErrorRetrieve()
{
    /* In case of truncation, separated by new line */
    if (gErrMsg[strlen(gErrMsg)-1] != '\n')
        gErrMsg[strlen(gErrMsg)-1] = '\n';

    return &gErrMsg[0];
}

