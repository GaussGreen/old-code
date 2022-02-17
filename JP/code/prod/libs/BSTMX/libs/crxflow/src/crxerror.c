/*********************************************************************************
 * CRXERROR.C
 * error utils
 *
 ********************************************************************************/

#include <stdio.h>
#include <string.h>
#include <math.h>

#include "crxerror.h"

#define ERRBUFFSIZE    1024L        /* size of error message */
static char gErrMsg[ERRBUFFSIZE];


#if(defined(WIN32)||defined(WINNT)||defined(_WIN32))
#define vsnprintf   _vsnprintf
#endif


static DRDiagnosticMsgCB ErrorCallbackFnp = NULL;
static DRDiagnosticMsgCB WarningCallbackFnp = NULL;
static DRDiagnosticMsgCB TraceCallbackFnp = NULL;

/**
 * set the error callback function as porposed by Jerry.
 */
void DR_ErrorSetCallBack(DRDiagnosticMsgCB errorCallbackFnp)
{
        /* Initialize error msg */
        strcpy (gErrMsg, "");
        ErrorCallbackFnp = errorCallbackFnp;
}

/**
 * place holder for other additional functions proposed by Jerry.
 */
void DR_WarningSetCallBack(DRDiagnosticMsgCB f)
{
        WarningCallbackFnp=f;
}
void DR_TraceSetCallBack(DRDiagnosticMsgCB f)
{
        TraceCallbackFnp=f;
}



/*****  DR_Error  ***********************************************************/
/*
*       Print an error message.
*/
void DR_Error(const char* format, ...)
{
        int handled = 0;
        va_list args;

        va_start(args, format);
        if (ErrorCallbackFnp != NULL)
        {
                handled = (*ErrorCallbackFnp)(format, args);
                (*ErrorCallbackFnp)("\n", NULL);
        }

        if (!handled)
        {
                vfprintf(stdout,format,args);
                fprintf(stdout,"\n");
                fflush(stdout);
        }
        va_end(args);
}



/*****  DR_ErrorCallback  ***************************************************/
/*
*       Store the error message in a string buffer
*/
int DR_ErrorCallback(const char* format, va_list ap)
{
    char    newErr[ERRBUFFSIZE];
    long    oldErrSize, newErrSizeMax;

    /* Existing error msg size */
    oldErrSize = strlen(gErrMsg);

    /* Max size of new error msg */
    newErrSizeMax = ERRBUFFSIZE - oldErrSize;

    /* Error bufffer full (including \0), no action */
    if (newErrSizeMax <= 1) 
        return TRUE;



    /* Print new error msg in the buffer, up to the full buffer including \0 
     */
    vsnprintf(newErr, newErrSizeMax , format, ap);

    /* Append to the existing error msg */
    strncat(gErrMsg, newErr, newErrSizeMax );

    return TRUE;

}


/*****  DR_ErrorRetrieve  ***************************************************/
/*
*       Retrieve the error message in a string buffer
*/
const char* DR_ErrorRetrieve()
{
    /* In case of truncation, separated by new line */
    if (gErrMsg[strlen(gErrMsg)-1] != '\n')
        gErrMsg[strlen(gErrMsg)-1] = '\n';

    return &gErrMsg[0];
}

