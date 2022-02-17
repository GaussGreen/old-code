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

#include "DRDiagnosticMsgCB.h"
#include "irx/platdep.h"
#include "irx/macros.h"

typedef struct 
{
    IrxTErrorCallBackV errorCallBackV;
    IrxTErrorCallBack  errorCallBack;
    IrxTTraceCallbackV traceCallBack;

} ERROR_GLOBAL;

static ERROR_GLOBAL global =
{
    NULL,
    NULL,
    NULL
};

#define DEFAULT_ERR_BUFF_SIZE 1024
static size_t gErrBuffSize = 0;
static char* gErrBuff = NULL;
static IrxTBool gAreErrorsEnabled = TRUE;
static IrxTBool gIsTracingEnabled = TRUE;
static FILE* gErrorFile = NULL;
static IrxTBool gFlushOnWrite = TRUE;

/*
** Some bogus code to force linking with esl_error.o and esl_warning.o  Otherwise,
** these callback registration functions may not be referenced anywhere
** and thus excluded by the linker.
*/
#ifdef __cplusplus
extern "C" {
#endif
typedef void (*DRDiagnosticMsgSetCB)(DRDiagnosticMsgCB);
#ifdef __cplusplus
}
#endif
void linkDRDiagnosticMsgCB()
{
  static DRDiagnosticMsgSetCB linkDR_ErrorSetCallBack = NULL;
  static DRDiagnosticMsgSetCB linkDR_WarningSetCallBack = NULL;
  static DRDiagnosticMsgSetCB linkDR_TraceSetCallBack = NULL;
  linkDR_ErrorSetCallBack = DR_ErrorSetCallBack;
  linkDR_WarningSetCallBack = DR_WarningSetCallBack;
  linkDR_TraceSetCallBack = DR_TraceSetCallBack;
}


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
** By default this prints to either an open error log file or stdout. However
** by setting a callback routine you can intercept the messages with your own
** message handler.
**
** The message printed should be terminated with a '\n' as this function
** does not add it.
***************************************************************************
*/
void irxErrorV(const char *format, va_list args)
{
    int doStdMsg = TRUE;
    FILE* errorDest = gErrorFile != NULL ? gErrorFile : stdout;

    if (!gAreErrorsEnabled)
        return;

    if (global.errorCallBack != NULL) 
    {
        int bufSize = gErrBuffSize > 0 ? gErrBuffSize : DEFAULT_ERR_BUFF_SIZE;
        char* errorBuff = (char*)malloc(bufSize);
        if (errorBuff != NULL) {
            vsnprintf(errorBuff, bufSize, format, args);
            if (global.errorCallBack(errorBuff))
                doStdMsg = FALSE;
            free(errorBuff);
        }
        else {
            fprintf(errorDest, "Failed to allocate buffer for error message.");
            /* And fall through to other error handlers... */
        }
    }

    if (global.errorCallBackV != NULL) 
    {
        if (global.errorCallBackV(format, args))
            doStdMsg = FALSE;
    }

    if (doStdMsg)
    {
        vfprintf(errorDest, format, args);
        if (gFlushOnWrite)
            fflush(errorDest);
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

void irxSetTraceCallBackV(IrxTTraceCallbackV callBackV)
{
    global.traceCallBack = callBackV;
}


/*------------------------------------------------------------------------
** E R R O R   B U F F E R
------------------------------------------------------------------------*/
/*f
***************************************************************************
**
** Sets size of buffer (in bytes) for storing error messages for later
** retrieval.  Any existing messages in the buffer will be lost!
**
***************************************************************************
*/
void irxSetErrorBufferSize(size_t bufSize) {
    if (gErrBuff != NULL)
    {
        free(gErrBuff);
        gErrBuff = NULL;
    }
    gErrBuffSize = bufSize;
}


static void allocErrorBuffer() {
    if (gErrBuffSize <= 0)
        gErrBuffSize = DEFAULT_ERR_BUFF_SIZE;
    gErrBuff = (char*)malloc(gErrBuffSize);
    if (gErrBuff != NULL)
        gErrBuff[0] = '\0';
}


/*f
***************************************************************************
**
** Stores the error message in a string buffer for later retrieval.
** Returns FAILURE if the message was truncated and/or the buffer is full.
**
***************************************************************************
*/
int irxBufferErrorMsg(const char* message) {
    static const char* truncMsg = "\nMESSAGE TRUNCATED\n";

    int rc = SUCCESS;
    size_t remaining; 

    if (gErrBuff == NULL)
        allocErrorBuffer();

    if (gErrBuff != NULL && message != NULL) {
        remaining = gErrBuffSize - strlen(gErrBuff) - 1;

        if (remaining > 0) {
            if (strlen(message) <= remaining)
                strcat(gErrBuff, message);
            else {
                rc = FAILURE;
                remaining -= strlen(truncMsg);
                if (remaining > 0)
                    strncat(gErrBuff, message, remaining);
                if (remaining >= 0)
                    strcat(gErrBuff, truncMsg);
            }
        }
    }

    return rc;
}


/*f
***************************************************************************
**
** Retrieves the buffered error message.
**
***************************************************************************
*/
const char* irxErrorRetrieve()
{
    return gErrBuff != NULL ? gErrBuff : "";
}


/*f
***************************************************************************
**
** Clears the error buffer.
**
***************************************************************************
*/
void irxClearErrorBuffer() {
    if (gErrBuff != NULL && gErrBuffSize > 0)
        gErrBuff[0] = '\0';
}


/*f
***************************************************************************
**
** Utility error callback that simply buffers error messages via
** irxBufferErrorMsg (see above).
**
***************************************************************************
*/
int irxErrorCallBack(const char* newErrMsg) {
    irxBufferErrorMsg(newErrMsg);
    return 1; /* Continue processing error message as usual */
}


/********************************************************************************
** Flags to control whether error and trace messages are handled or ignored and
** whether to flush messages after each write.
********************************************************************************/
IrxTBool irxEnableErrorMsgs(IrxTBool isEnabled) {
    IrxTBool prevSetting = gAreErrorsEnabled;
    gAreErrorsEnabled = isEnabled;
    return prevSetting;
}

IrxTBool irxEnableTraceMsgs(IrxTBool isEnabled) {
    IrxTBool prevSetting = gIsTracingEnabled;
    gIsTracingEnabled = isEnabled;
    return prevSetting;
}

IrxTBool irxSetFlushOnWriteMsg(IrxTBool isEnabled) {
    IrxTBool prevSetting = gFlushOnWrite;
    gFlushOnWrite = isEnabled;
    return prevSetting;
}

/********************************************************************************
** Manual flushing of OS-buffered error and trace message.
********************************************************************************/
void irxFlushErrorMsgs() {
    if (gErrorFile != NULL)
        fflush(gErrorFile);
    else
        fflush(stdout);
}

void irxFlushTraceMsgs() {
    fflush(stdout);
}

/********************************************************************************
** Error file management.
********************************************************************************/
int irxOpenErrorFile(const char* filename, int append) {
    int rc = FAILURE;

    irxCloseErrorFile();

    if (filename != NULL) {
        if ((gErrorFile = fopen(filename, append ? "a" : "w")) != NULL)
            rc = SUCCESS;
        else
            irxError("Could not open error file: %s\n", filename);
    }
    else
        irxError("NULL filename passed to irxOpenErrorFile.\n");

    return rc;
}

int irxCloseErrorFile() {
    if (gErrorFile != NULL)
    {
        fclose(gErrorFile);
        gErrorFile = NULL;
    }

    return SUCCESS;
}


/*****  DR_Trace  ***********************************************************/
/*
*       Print an trace message.
*/
void DR_Trace(char* format, ...)
{
	int handled = 0;
	va_list args;

    if (!gIsTracingEnabled)
        return;

	va_start(args, format);

	if (global.traceCallBack != NULL)
    {
		handled = (*global.traceCallBack)(format, args);
		(*global.traceCallBack)("\n", NULL);
	}

	if (!handled)
    {
		vfprintf(stdout,format,args);
		fprintf(stdout,"\n");
        if (gFlushOnWrite)
		    fflush(stdout);
	}
	va_end(args);
}
