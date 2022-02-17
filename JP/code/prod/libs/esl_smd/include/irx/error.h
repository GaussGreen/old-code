/*
***************************************************************************
** HEADER FILE: error.h
**
** Error handling for the credit & rates exotics.
***************************************************************************
*/

#ifndef IRX_ERROR_H
#define IRX_ERROR_H

#include <stdarg.h>          /* defines va_list */

/**
***************************************************************************
** Defines a call back routine for a message handler. This version takes
** the formatted message.
**
** Routines of this type can be registered with the error handling
** mechanism, and are called before the standard error handling. In such
** cases, if the call back routine returns TRUE, then the message is no
** longer sent to stderr (the default action of irxError).
***************************************************************************
*/
typedef int (*IrxTErrorCallBack)(const char *);

/*t
***************************************************************************
** Defines a call back routine for a message handler. This version takes
** the unformatted message and the variable argument list.
**
** Routines of this type can be registered with the error handling
** mechanism, and are called before the standard error handling. In such
** cases, if the call back routine returns TRUE, then the message is no
** longer sent to stderr (the default action of irxError).
***************************************************************************
*/
typedef int (*IrxTErrorCallBackV)(const char *, va_list);

/**
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
#if __GNUC__ > 1
__attribute__ ((format(printf,1,2)))
#endif
    ;

/**
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
void irxErrorV(const char *format, va_list args);

/**
***************************************************************************
** Handles the message that a routine failed. Designed so that this message
** is consistent throughout the library. Returns FAILURE, so that you can
** end a routine with:
\begin{verbatim}
    return irxErrorFailure(routine);
\end{verbatim}
***************************************************************************
*/
int irxErrorFailure(const char* routine);

/**
***************************************************************************
** Sets an error callback routine.
**
** This callback routine takes the formatted string.
** If the callback returns TRUE, then the message is no longer printed
** to stderr as well.
**
** See also irxSetErrorCallBackV
***************************************************************************
*/
void irxSetErrorCallBack(IrxTErrorCallBack callBack);

/**
***************************************************************************
** Sets error callback routines.
**
** You provide two types of routines.
**
** The first one takes the format string and va_list arguments. 
** If the callback returns TRUE, then the message is no longer printed
** to stderr as well. This is used when irxError or irxErrorV is called.
**
** The second one takes a pre-formatted message without arguments.
** If the callback returns FALSE, then the message is no longer printed
** to stderr as well. This is used when irxErrorWrite is called.
**
** See also irxSetErrorCallBack
***************************************************************************
*/
void irxSetErrorCallBackV(IrxTErrorCallBackV callBackV);


/**
***************************************************************************
**
**       Store the error message in a string buffer
**
***************************************************************************
*/
void irxErrorCallBack(const char* newErrMsg);


/**
***************************************************************************
**
**       Retrieve the error message in a string buffer
**       Store the error message in a string buffer
**
***************************************************************************
*/
const char* irxErrorRetrieve();


#endif

