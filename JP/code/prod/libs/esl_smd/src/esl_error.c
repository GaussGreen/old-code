/****************************************************************************/
/*      Error message handling.                                             */
/****************************************************************************/
/*      ERROR.c                                                             */
/****************************************************************************/
/*
$Header$
*/
#include <stdio.h>
#include <stdlib.h>  
#include <math.h>
#include <string.h>

#include <DRDiagnosticMsgCB.h>

static DRDiagnosticMsgCB ErrorCallbackFnp = NULL;
static DRDiagnosticMsgCB WarningCallbackFnp = NULL;
static DRDiagnosticMsgCB TraceCallbackFnp = NULL;

/* Callback set by the ESL Log facility */
DRDiagnosticMsgCB EslLogCallbackFnp = NULL;

/**
 * set the error callback function as porposed by Jerry.
 */
void DR_ErrorSetCallBack(DRDiagnosticMsgCB errorCallbackFnp)
{
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
/**
*       Print an error message.
*/

void DR_Error(char const* format, ...)
{
	int handled = 0;
	va_list args;

	va_start(args, format);

    /* if ESL Log enabled print only there */
    if (EslLogCallbackFnp)
    {
		(*EslLogCallbackFnp)(format, args);
    }
    else
    {
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
	}
	va_end(args);
}


