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
#include <memory.h>          /* for size_t */
#include <irx/cgeneral.h>    /* for IrxTBool */

#ifdef  __cplusplus
extern "C" {
#endif

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


typedef int (*IrxTTraceCallbackV)(const char*, va_list);

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

void irxSetTraceCallBackV(IrxTTraceCallbackV callBackV);


/**
**************************************************************************
**
** If space is available, stores message in internal buffer for later
** retrieval by irxErrorRetrieve.  Message may be truncated in buffer.
** Returns FAILURE if message was truncated and/or buffer is full.
**
***************************************************************************
*/
int irxBufferErrorMsg(const char* message);

/**
***************************************************************************
**
** Sets size of buffer (in bytes) for storing error messages for later
** retrieval.  Any existing messages in the buffer will be lost!
**
***************************************************************************
*/
void irxSetErrorBufferSize(size_t bufSize);

/**
***************************************************************************
**
**       Retrieves any error message that was stored by irxErrorCallBack.
**
***************************************************************************
*/
const char* irxErrorRetrieve();

/**
***************************************************************************
**
** Clears the error buffer.
**
***************************************************************************
*/
void irxClearErrorBuffer();

/**
***************************************************************************
**
** Utility error callback that simply buffers error messages via
** irxBufferErrorMsg (see above).
**
***************************************************************************
*/
int irxErrorCallBack(const char* newErrMsg);

/**
***************************************************************************
*
*       Print an trace message.
*
***************************************************************************
*/
void DR_Trace(char* format, ...);

/**
***************************************************************************
*
*       Enable/disable error messages.  Returns previous setting.
*
***************************************************************************
*/
IrxTBool irxEnableErrorMsgs(IrxTBool isEnabled);

/**
***************************************************************************
*
*       Enable/disable trace messages.  Returns previous setting.
*
***************************************************************************
*/
IrxTBool irxEnableTraceMsgs(IrxTBool isEnabled);

/**
***************************************************************************
*
*       Enable/disable flushing of error and trace messages each time they
*       are output .  Returns previous setting.
*
***************************************************************************
*/
IrxTBool irxSetFlushOnWriteMsg(IrxTBool isEnabled);

/**
***************************************************************************
*
*       Manually flushes any OS-buffered error/trace messages, respectively.
*
***************************************************************************
*/
void irxFlushErrorMsgs();
void irxFlushTraceMsgs();

/**
***************************************************************************
*
*       Open/close error file.  If opened, it will be used instead of stdout
*       as the destination for error messages (if enabled) that are not
*       completely handled by any error callback function.  Returns
*       SUCCESS or FAILURE.
*
***************************************************************************
*/
int irxOpenErrorFile(const char* filename, int append);
int irxCloseErrorFile();

#ifdef  __cplusplus
}
#endif

#endif
